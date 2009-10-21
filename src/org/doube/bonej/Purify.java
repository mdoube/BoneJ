package org.doube.bonej;

/**
 *  
 * Purify_ plugin for ImageJ
 * Copyright 2009 Michael Doube 
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

import ij.*;
import ij.plugin.PlugIn;
import ij.measure.ResultsTable;
import ij.gui.GenericDialog;
import java.util.Arrays;

/**
 * <p>
 * Purify_ plugin for ImageJ
 * </p>
 * 
 * <p>
 * Prepare binary stack for connectivity analysis by reducing number of
 * reference phase (foreground) particles to 1, filling cavities within the
 * single reference phase particle and ensuring there is only 1 particle in the
 * background phase.
 * </p>
 * 
 * <p>
 * Foreground is 26-connected and background is 6-connected.
 * </p>
 * 
 * <p>
 * This plugin is based on Object_Counter3D by Fabrice P Cordelires and Jonathan
 * Jackson, but with significant speed increases through reduction of recursion
 * and multi-threading. Thanks to Robert Barbour for the suggestion to 'chunk'
 * the stack. Chunking works as follows:
 * </p>
 * <ol>
 * <li>Perform initial labelling on the whole stack in a single thread</li>
 * <li>for <i>n</i> discrete, contiguous chunks within the labelling array,
 * connectStructures()
 * <ol type="a">
 * <li>connectStructures() can run in a separate thread for each chunk</li>
 * <li>chunks are approximately equal-sized sets of slices</li>
 * </ol>
 * <li>stitchChunks() for the pixels on the first slice of each chunk, except
 * for the first chunk, restricting replaceLabels() to the current and all
 * previous chunks.
 * <ol type="a">
 * <li>stitchChunks() iterates through the slice being stitched in a single
 * thread</li>
 * </ol>
 * </li>
 * 
 * </ol>
 * <p>
 * The performance improvement should be in the region of a factor of <i>n</i>
 * if run linearly, and if multithreaded over <i>c</i> processors, speed
 * increase should be in the region of <i>n</i> * <i>c</i>, minus overhead.
 * </p>
 * 
 * @author Michael Doube
 * @version 1.0
 * @see <p>
 *      Odgaard A, Gundersen HJG (1993) Quantification of connectivity in
 *      cancellous bone, with special emphasis on 3-D reconstructions. Bone 14:
 *      173-182. <a
 *      href="http://dx.doi.org/10.1016/8756-3282(93)90245-6">doi:10.1016
 *      /8756-3282(93)90245-6</a>
 *      </p>
 *      <p>
 *      Object counter 3D <a
 *      href="http://rsbweb.nih.gov/ij/plugins/track/objects.html">home page at
 *      the NIH</a>
 *      </p>
 * 
 */
public class Purify implements PlugIn {

	private final static int foreground = -1, background = 0;

	private int width, height, nSlices, nThreads, nChunks, sliceSize;

	private int ID;

	private String chunkString = "";

	private String sPhase = "";

	public void run(String arg) {
		ImagePlus imp = IJ.getImage();
		ImageCheck ic = new ImageCheck();
		if (!ic.isBinary(imp)){
			IJ.error("Purify requires a binary image");
			return;
		}
		GenericDialog gd = new GenericDialog("Setup");
		gd.addNumericField("Chunk Size", 4, 0, 4, "slices");
		gd.addCheckbox("Performance Log", false);
		gd.addCheckbox("Make copy", true);
		gd.showDialog();
		int slicesPerChunk = (int) Math.floor(gd.getNextNumber());
		boolean showPerformance = gd.getNextBoolean();
		boolean doCopy = gd.getNextBoolean();
		if (gd.wasCanceled()) {
			return;
		}

		Object[] result = purify(imp, slicesPerChunk, doCopy, showPerformance);
		IJ.freeMemory();
		if (null != result) {
			ImagePlus purified = (ImagePlus) result[1];

			if (doCopy){
				purified.show();
				if (!purified.isInvertedLut())
					IJ.run("Invert LUT");
			}
			else {
				ImageStack stack2 = purified.getStack();
				imp.setStack(imp.getTitle(), stack2);
				imp.show();
				if (!imp.isInvertedLut())
					IJ.run("Invert LUT");
			}
		}
		System.gc();
		return;
	}

	public Object[] purify(ImagePlus imp, int slicesPerChunk, boolean doCopy,
			boolean showPerformance) {
		this.width = imp.getWidth();
		this.height = imp.getHeight();
		this.nSlices = imp.getStackSize();
		this.sliceSize = this.width * this.height;
		nThreads = Runtime.getRuntime().availableProcessors();
		nChunks = (int) Math.floor((double) nSlices / (double) slicesPerChunk);
		if (nChunks == 0)
			nChunks = 1;

		int remainder = nSlices - nChunks * slicesPerChunk;

		if (remainder > 0) {
			nChunks += 1 + (int) Math.floor((double) remainder
					/ (double) slicesPerChunk);
		}

		// set up the chunks
		final int[][] chunkRanges = getChunkRanges(nChunks, slicesPerChunk);
		final int[][] stitchRanges = getStitchRanges(nChunks, slicesPerChunk);

		long startTime = System.currentTimeMillis();

		byte[][] workArray = makeWorkArray(imp);

		sPhase = "foreground";
		int[][] particleLabels = firstIDAttribution(workArray, foreground);

		// connect foreground particles within chunks
		ConnectStructuresThread[] cptf = new ConnectStructuresThread[nThreads];
		for (int thread = 0; thread < nThreads; thread++) {
			cptf[thread] = new ConnectStructuresThread(thread, nThreads,
					workArray, particleLabels, foreground, nChunks, chunkRanges);
			cptf[thread].start();
		}
		try {
			for (int thread = 0; thread < nThreads; thread++) {
				cptf[thread].join();
			}
		} catch (InterruptedException ie) {
			IJ.error("A thread was interrupted.");
		}

		// connect foreground particles between chunks
		if (nChunks > 1) {
			chunkString = ": stitching...";
			connectStructures(workArray, particleLabels, foreground,
					stitchRanges);
		}

		long[] particleSizes = getParticleSizes(workArray, particleLabels,
				foreground);
		removeSmallParticles(workArray, particleLabels, particleSizes,
				foreground);

		Arrays.fill(particleSizes, 0);
		for (int i = 0; i < particleLabels.length; i++) {
			Arrays.fill(particleLabels[i], 0);
		}
		IJ.freeMemory();

		sPhase = "background";
		particleLabels = firstIDAttribution(workArray, background);

		// connect background particles within chunks
		ConnectStructuresThread[] cptb = new ConnectStructuresThread[nThreads];
		for (int thread = 0; thread < nThreads; thread++) {
			cptb[thread] = new ConnectStructuresThread(thread, nThreads,
					workArray, particleLabels, background, nChunks, chunkRanges);
			cptb[thread].start();
		}
		try {
			for (int thread = 0; thread < nThreads; thread++) {
				cptb[thread].join();
			}
		} catch (InterruptedException ie) {
			IJ.error("A thread was interrupted.");
		}

		// connect background particles between chunks
		if (nChunks > 1) {
			chunkString = ": stitching...";
			connectStructures(workArray, particleLabels, background,
					stitchRanges);
		}
		particleSizes = getParticleSizes(workArray, particleLabels, background);
		touchEdges(workArray, particleLabels, particleSizes, background);
		particleSizes = getParticleSizes(workArray, particleLabels, background);
		removeSmallParticles(workArray, particleLabels, particleSizes,
				background);

		double duration = ((double) System.currentTimeMillis() - (double) startTime)
				/ (double) 1000;
		if (showPerformance)
			showResults(chunkRanges, duration, imp, slicesPerChunk);

		IJ.showStatus("Image Purified");

		ImageStack stack = new ImageStack(this.width, this.height);
		for (int z = 0; z < workArray.length; z++) {
			stack.addSlice("", workArray[z]);
		}
		ImagePlus purified = new ImagePlus("Purified", stack);
		purified.setCalibration(imp.getCalibration());
		Object[] result = { duration, purified };
		IJ.freeMemory();
		System.gc();
		return result;
	}

	/**
	 * Create a work array
	 * 
	 * @return byte[] work array
	 */
	private byte[][] makeWorkArray(ImagePlus imp) {
		int s = imp.getStackSize();
		int p = imp.getWidth() * imp.getHeight();
		byte[][] workArray = new byte[s][p];
		ImageStack stack = imp.getStack();
		for (int z = 0; z < s; z++) {
			byte[] slicePixels = (byte[]) stack.getPixels(z + 1);
			System.arraycopy(slicePixels, 0, workArray[z], 0, p);
		}
		return workArray;
	}

	/**
	 * Get a 2 d array that defines the z-slices to scan within while connecting
	 * particles within chunkified stacks.
	 * 
	 * @param nC
	 *            number of chunks
	 * @return scanRanges int[][] containing 4 limits: int[0][] - start of outer
	 *         for; int[1][] end of outer for; int[3][] start of inner for;
	 *         int[4] end of inner 4. Second dimension is chunk number.
	 */
	private int[][] getChunkRanges(int nC, int slicesPerChunk) {
		int[][] scanRanges = new int[4][nC];
		scanRanges[0][0] = 0; // the first chunk starts at the first (zeroth)
		// slice
		scanRanges[2][0] = 0; // and that is what replaceLabel() will work on
		// first

		if (nC == 1) {
			scanRanges[1][0] = nSlices;
			scanRanges[3][0] = nSlices;
		} else if (nC > 1) {
			scanRanges[1][0] = slicesPerChunk;
			scanRanges[3][0] = slicesPerChunk;

			for (int c = 1; c < nC; c++) {
				for (int i = 0; i < 4; i++) {
					scanRanges[i][c] = scanRanges[i][c - 1] + slicesPerChunk;
				}
			}
			// reduce the last chunk to nSlices
			scanRanges[1][nC - 1] = nSlices;
			scanRanges[3][nC - 1] = nSlices;
		}
		return scanRanges;
	}

	/**
	 * Return scan ranges for stitching. The first 2 values for each chunk are
	 * the first slice of the next chunk and the last 2 values are the range
	 * through which to replaceLabels()
	 * 
	 * Running replace labels over incrementally increasing volumes as chunks
	 * are added is OK (for 1st interface connect chunks 0 & 1, for 2nd connect
	 * chunks 0, 1, 2, etc.)
	 * 
	 * @param nC
	 *            number of chunks
	 * @return scanRanges list of scan limits for connectStructures() to stitch
	 *         chunks back together
	 */
	private int[][] getStitchRanges(int nC, int slicesPerChunk) {
		if (nC < 2) {
			return null;
		}
		int[][] scanRanges = new int[4][3 * (nC - 1)]; // there are nC - 1
		// interfaces

		for (int c = 0; c < nC - 1; c++) {
			scanRanges[0][c] = (c + 1) * slicesPerChunk;
			scanRanges[1][c] = (c + 1) * slicesPerChunk + 1;
			scanRanges[2][c] = c * slicesPerChunk; // forward and reverse
			// algorithm
			// scanRanges[2][c] = 0; //cumulative algorithm - reliable but O^2
			// hard
			scanRanges[3][c] = (c + 2) * slicesPerChunk;
		}
		// stitch back
		for (int c = nC - 1; c < 2 * (nC - 1); c++) {
			scanRanges[0][c] = (2 * nC - c - 2) * slicesPerChunk - 1;
			scanRanges[1][c] = (2 * nC - c - 2) * slicesPerChunk;
			scanRanges[2][c] = (2 * nC - c - 3) * slicesPerChunk;
			scanRanges[3][c] = (2 * nC - c - 1) * slicesPerChunk;
		}
		// stitch forwards (paranoid third pass)
		for (int c = 2 * (nC - 1); c < 3 * (nC - 1); c++) {
			scanRanges[0][c] = (-2 * nC + c + 3) * slicesPerChunk;
			scanRanges[1][c] = (-2 * nC + c + 3) * slicesPerChunk + 1;
			scanRanges[2][c] = (-2 * nC + c + 2) * slicesPerChunk;
			scanRanges[3][c] = (-2 * nC + c + 4) * slicesPerChunk;
		}
		for (int i = 0; i < scanRanges.length; i++) {
			for (int c = 0; c < scanRanges[i].length; c++) {
				if (scanRanges[i][c] > nSlices) {
					scanRanges[i][c] = nSlices;
				}
			}
		}
		scanRanges[3][nC - 2] = nSlices;
		return scanRanges;
	}

	/**
	 * Go through all pixels and assign initial particle label
	 * 
	 * @param workArray
	 *            byte[] array containing pixel values
	 * @param phase
	 *            foreground = -1, background = 0
	 * @return particleLabels int[] array containing label associating every
	 *         pixel with a particle
	 */
	private int[][] firstIDAttribution(byte[][] workArray, int phase) {
		IJ.showStatus("Finding " + sPhase + " structures");
		int[][] particleLabels = new int[nSlices][sliceSize];
		ID = 1;

		if (phase == foreground) {
			for (int z = 0; z < nSlices; z++) {
				for (int y = 0; y < height; y++) {
					int rowIndex = y * width;
					for (int x = 0; x < width; x++) {
						int arrayIndex = rowIndex + x;
						if (workArray[z][arrayIndex] == phase) {
							particleLabels[z][arrayIndex] = ID;
							int minTag = ID;
							// Find the minimum particleLabel in the
							// neighbouring pixels
							for (int vZ = z - 1; vZ <= z + 1; vZ++) {
								for (int vY = y - 1; vY <= y + 1; vY++) {
									for (int vX = x - 1; vX <= x + 1; vX++) {
										if (withinBounds(vX, vY, vZ, 0, nSlices)) {
											int offset = getOffset(vX, vY);
											if (workArray[vZ][offset] == phase) {
												int tagv = particleLabels[vZ][offset];
												if (tagv != 0 && tagv < minTag) {
													minTag = tagv;
												}
											}
										}
									}
								}
							}
							// assign the smallest particle label from the
							// neighbours to the pixel
							particleLabels[z][arrayIndex] = minTag;
							// increment the particle label
							if (minTag == ID) {
								ID++;
							}
						}
					}
				}
				IJ.showProgress(z, nSlices);
			}
			ID++;
		} else if (phase == background) {
			for (int z = 0; z < nSlices; z++) {
				for (int y = 0; y < height; y++) {
					int rowIndex = y * width;
					for (int x = 0; x < width; x++) {
						int arrayIndex = rowIndex + x;
						if (workArray[z][arrayIndex] == phase) {
							particleLabels[z][arrayIndex] = ID;
							int minTag = ID;
							// Find the minimum particleLabel in the
							// neighbouring pixels
							int nX = x, nY = y, nZ = z;
							for (int n = 0; n < 7; n++) {
								switch (n) {
								case 0:
									break;
								case 1:
									nX = x - 1;
									break;
								case 2:
									nX = x + 1;
									break;
								case 3:
									nY = y - 1;
									nX = x;
									break;
								case 4:
									nY = y + 1;
									break;
								case 5:
									nZ = z - 1;
									nY = y;
									break;
								case 6:
									nZ = z + 1;
									break;
								}
								if (withinBounds(nX, nY, nZ, 0, nSlices)) {
									int offset = getOffset(nX, nY);
									if (workArray[nZ][offset] == phase) {
										int tagv = particleLabels[nZ][offset];
										if (tagv != 0 && tagv < minTag) {
											minTag = tagv;
										}
									}
								}
							}
							// assign the smallest particle label from the
							// neighbours to the pixel
							particleLabels[z][arrayIndex] = minTag;
							// increment the particle label
							if (minTag == ID) {
								ID++;
							}
						}
					}
				}
				IJ.showProgress(z, nSlices);
			}
			ID++;
		}
		return particleLabels;
	}

	/**
	 * Connect structures = minimisation of IDs
	 * 
	 * @param workArray
	 * @param particleLabels
	 * @param phase
	 *            foreground or background
	 * @param scanRanges
	 *            int[][] listing ranges to run connectStructures on
	 * @return particleLabels with all particles connected
	 */
	private void connectStructures(byte[][] workArray, int[][] particleLabels,
			int phase, int[][] scanRanges) {
		IJ.showStatus("Connecting " + sPhase + " structures" + chunkString);
		for (int c = 0; c < scanRanges[0].length; c++) {
			int sR0 = scanRanges[0][c];
			int sR1 = scanRanges[1][c];
			int sR2 = scanRanges[2][c];
			int sR3 = scanRanges[3][c];
			if (phase == foreground) {
				for (int z = sR0; z < sR1; z++) {
					for (int y = 0; y < height; y++) {
						int rowIndex = y * width;
						for (int x = 0; x < width; x++) {
							int arrayIndex = rowIndex + x;
							if (workArray[z][arrayIndex] == phase
									&& particleLabels[z][arrayIndex] > 1) {
								int minTag = particleLabels[z][arrayIndex];
								// Find the minimum particleLabel in the
								// neighbours' pixels
								for (int vZ = z - 1; vZ <= z + 1; vZ++) {
									for (int vY = y - 1; vY <= y + 1; vY++) {
										for (int vX = x - 1; vX <= x + 1; vX++) {
											if (withinBounds(vX, vY, vZ, sR2,
													sR3)) {
												int offset = getOffset(vX, vY);
												if (workArray[vZ][offset] == phase) {
													int tagv = particleLabels[vZ][offset];
													if (tagv != 0
															&& tagv < minTag) {
														minTag = tagv;
													}
												}
											}
										}
									}
								}
								// Replacing particleLabel by the minimum
								// particleLabel found
								for (int vZ = z - 1; vZ <= z + 1; vZ++) {
									for (int vY = y - 1; vY <= y + 1; vY++) {
										for (int vX = x - 1; vX <= x + 1; vX++) {
											if (withinBounds(vX, vY, vZ, sR2,
													sR3)) {
												int offset = getOffset(vX, vY);
												if (workArray[vZ][offset] == phase) {
													int tagv = particleLabels[vZ][offset];
													if (tagv != 0
															&& tagv != minTag) {
														replaceLabel(
																particleLabels,
																tagv, minTag,
																sR2, sR3);
													}
												}
											}
										}
									}
								}
							}
						}
					}
					IJ.showStatus("Connecting foreground structures"
							+ chunkString);
					IJ.showProgress(z, nSlices);
				}
			} else if (phase == background) {
				for (int z = sR0; z < sR1; z++) {
					for (int y = 0; y < height; y++) {
						int rowIndex = y * width;
						for (int x = 0; x < width; x++) {
							int arrayIndex = rowIndex + x;
							if (workArray[z][arrayIndex] == phase) {
								int minTag = particleLabels[z][arrayIndex];
								// Find the minimum particleLabel in the
								// neighbours' pixels
								int nX = x, nY = y, nZ = z;
								for (int n = 0; n < 7; n++) {
									switch (n) {
									case 0:
										break;
									case 1:
										nX = x - 1;
										break;
									case 2:
										nX = x + 1;
										break;
									case 3:
										nY = y - 1;
										nX = x;
										break;
									case 4:
										nY = y + 1;
										break;
									case 5:
										nZ = z - 1;
										nY = y;
										break;
									case 6:
										nZ = z + 1;
										break;
									}
									if (withinBounds(nX, nY, nZ, sR2, sR3)) {
										int offset = getOffset(nX, nY);
										if (workArray[nZ][offset] == phase) {
											int tagv = particleLabels[nZ][offset];
											if (tagv != 0 && tagv < minTag) {
												minTag = tagv;
											}
										}
									}
								}
								// Replacing particleLabel by the minimum
								// particleLabel found
								for (int n = 0; n < 7; n++) {
									switch (n) {
									case 0:
										nZ = z;
										break; // last switch block left nZ = z
									// + 1;
									case 1:
										nX = x - 1;
										break;
									case 2:
										nX = x + 1;
										break;
									case 3:
										nY = y - 1;
										nX = x;
										break;
									case 4:
										nY = y + 1;
										break;
									case 5:
										nZ = z - 1;
										nY = y;
										break;
									case 6:
										nZ = z + 1;
										break;
									}
									if (withinBounds(nX, nY, nZ, sR2, sR3)) {
										int offset = getOffset(nX, nY);
										if (workArray[nZ][offset] == phase) {
											int tagv = particleLabels[nZ][offset];
											if (tagv != 0 && tagv != minTag) {
												replaceLabel(particleLabels,
														tagv, minTag, sR2, sR3);
											}
										}
									}
								}
							}
						}
					}
					IJ.showStatus("Connecting background structures"
							+ chunkString);
					IJ.showProgress(z, nSlices + 1);
				}
			}
		}
		return;
	}

	class ConnectStructuresThread extends Thread {
		final int thread, nThreads, nChunks, phase;// w,h,d,nR;

		final byte[][] workArray;

		final int[][] particleLabels;

		final int[][] chunkRanges;

		public ConnectStructuresThread(int thread, int nThreads,
				byte[][] workArray, int[][] particleLabels, int phase,
				int nChunks, int[][] chunkRanges) {
			this.thread = thread;
			this.nThreads = nThreads;
			this.workArray = workArray;
			this.particleLabels = particleLabels;
			this.phase = phase;
			this.nChunks = nChunks;
			this.chunkRanges = chunkRanges;
		}

		public void run() {
			for (int k = this.thread; k < this.nChunks; k += this.nThreads) {
				// assign singleChunkRange for chunk k from chunkRanges
				int[][] singleChunkRange = new int[4][1];
				for (int i = 0; i < 4; i++) {
					singleChunkRange[i][0] = this.chunkRanges[i][k];
				}
				chunkString = ": chunk " + (k + 1) + "/" + nChunks;
				connectStructures(this.workArray, this.particleLabels,
						this.phase, singleChunkRange);
			}
		}
	}// ConnectStructuresThread

	/**
	 * Get the sizes of all the particles
	 * 
	 * @param workArray
	 * @param particleLabels
	 * @param phase
	 *            foreground = -1, background = 0
	 * @return particleSizes
	 */
	private long[] getParticleSizes(byte[][] workArray, int[][] particleLabels,
			int phase) {
		IJ.showStatus("Getting " + sPhase + " particle sizes");
		long[] particleSizes = new long[ID];
		for (int z = 0; z < nSlices; z++) {
			int arrayIndex = 0;
			for (int y = 0; y < height; y++) {
				for (int x = 0; x < width; x++) {
					if (workArray[z][arrayIndex] == phase) {
						particleSizes[particleLabels[z][arrayIndex]]++;
					}
					arrayIndex++;
				}
			}
			IJ.showProgress(z, nSlices);
		}
		return particleSizes;
	}

	/**
	 * <p>
	 * Find particles of phase that touch the stack sides and assign them the ID
	 * of the biggest particle of phase. Euler number calculation assumes that
	 * the background phase is connected outside the image stack, so apparently
	 * isolated background particles touching the sides should be assigned to
	 * the single background particle.
	 * </p>
	 * 
	 * @param workArray
	 * @param particleLabels
	 * @param particleSizes
	 * @param phase
	 * @return particleLabels
	 */
	private void touchEdges(byte[][] workArray, int[][] particleLabels,
			long[] particleSizes, int phase) {
		String status = "Background particles touching ";
		// find the label associated with the biggest
		// particle in phase
		long maxVoxCount = 0;
		int biggestParticle = 0;
		for (int i = 0; i < particleSizes.length; i++) {
			if (particleSizes[i] > maxVoxCount) {
				maxVoxCount = particleSizes[i];
				biggestParticle = i;
			}
		}
		// check each face of the stack for pixels that are touching edges and
		// replace that particle's label in particleLabels with
		// the label of the biggest particle
		int x, y, z;

		// up
		z = 0;
		for (y = 0; y < height; y++) {
			IJ.showStatus(status + "top");
			IJ.showProgress(y, height);
			int rowOffset = y * width;
			for (x = 0; x < width; x++) {
				int offset = rowOffset + x;
				if (workArray[z][offset] == phase
						&& particleLabels[z][offset] != biggestParticle) {
					replaceLabel(particleLabels, particleLabels[z][offset],
							biggestParticle, 0, nSlices);
				}
			}
		}

		// down
		z = nSlices - 1;
		for (y = 0; y < height; y++) {
			IJ.showStatus(status + "bottom");
			IJ.showProgress(y, height);
			int rowOffset = y * width;
			for (x = 0; x < width; x++) {
				int offset = rowOffset + x;
				if (workArray[z][offset] == phase
						&& particleLabels[z][offset] != biggestParticle) {
					replaceLabel(particleLabels, particleLabels[z][offset],
							biggestParticle, 0, nSlices);
				}
			}
		}

		// left
		x = 0;
		for (z = 0; z < nSlices; z++) {
			IJ.showStatus(status + "left");
			IJ.showProgress(z, nSlices);
			for (y = 0; y < height; y++) {
				int offset = y * width;
				if (workArray[z][offset] == phase
						&& particleLabels[z][offset] != biggestParticle) {
					replaceLabel(particleLabels, particleLabels[z][offset],
							biggestParticle, 0, nSlices);
				}
			}
		}

		// right
		x = width - 1;
		for (z = 0; z < nSlices; z++) {
			IJ.showStatus(status + "right");
			IJ.showProgress(z, nSlices);
			for (y = 0; y < height; y++) {
				int offset = y * width + x;
				if (workArray[z][offset] == phase
						&& particleLabels[z][offset] != biggestParticle) {
					replaceLabel(particleLabels, particleLabels[z][offset],
							biggestParticle, 0, nSlices);
				}
			}
		}

		// front
		y = height - 1;
		int rowOffset = y * width;
		for (z = 0; z < nSlices; z++) {
			IJ.showStatus(status + "front");
			IJ.showProgress(z, nSlices);
			for (x = 0; x < width; x++) {
				int offset = rowOffset + x;
				if (workArray[z][offset] == phase
						&& particleLabels[z][offset] != biggestParticle) {
					replaceLabel(particleLabels, particleLabels[z][offset],
							biggestParticle, 0, nSlices);
				}
			}
		}

		// back
		y = 0;
		for (z = 0; z < nSlices; z++) {
			IJ.showStatus(status + "back");
			IJ.showProgress(z, nSlices);
			for (x = 0; x < width; x++) {
				int offset = x;
				if (workArray[z][offset] == phase
						&& particleLabels[z][offset] != biggestParticle) {
					replaceLabel(particleLabels, particleLabels[z][offset],
							biggestParticle, 0, nSlices);
				}
			}
		}
		return;
	}

	/**
	 * Remove all but the largest phase particle from workArray
	 * 
	 * @param workArray
	 * @param particleLabels
	 * @param particleSizes
	 * @param phase
	 * @return workArray
	 */
	private void removeSmallParticles(byte[][] workArray,
			int[][] particleLabels, long[] particleSizes, int phase) {
		long maxVoxCount = 0;
		for (int i = 0; i < particleSizes.length; i++) {
			if (particleSizes[i] > maxVoxCount) {
				maxVoxCount = particleSizes[i];
			}
		}
		if (phase == foreground) {
			// go through work array and turn all
			// smaller foreground particles into background (0)
			for (int z = 0; z < nSlices; z++) {
				for (int i = 0; i < sliceSize; i++) {
					if (workArray[z][i] == foreground) {
						if (particleSizes[particleLabels[z][i]] < maxVoxCount) {
							workArray[z][i] = background;
						}
					}
				}
				IJ.showStatus("Removing foreground particles");
				IJ.showProgress(z, nSlices);
			}
		} else if (phase == background) {
			// go through work array and turn all
			// smaller background particles into foreground
			for (int z = 0; z < nSlices; z++) {
				for (int i = 0; i < sliceSize; i++) {
					if (workArray[z][i] == background) {
						if (particleSizes[particleLabels[z][i]] < maxVoxCount) {
							workArray[z][i] = foreground;
						}
					}
				}
				IJ.showStatus("Removing background particles");
				IJ.showProgress(z, nSlices);
			}
		}
		return;
	}

	/**
	 * Check to see if the pixel at (m,n,o) is within the bounds of the current
	 * stack
	 * 
	 * @param m
	 *            x co-ordinate
	 * @param n
	 *            y co-ordinate
	 * @param o
	 *            z co-ordinate
	 * @param startZ
	 *            first Z coordinate to use
	 * 
	 * @param endZ
	 *            last Z coordinate to use
	 * 
	 * @return True if the pixel is within the bounds of the current stack
	 */
	private boolean withinBounds(int m, int n, int o, int startZ, int endZ) {
		return (m >= 0 && m < width && n >= 0 && n < height && o >= startZ && o < endZ);
	}

	/**
	 * Find the offset within a 1D array given 2 (x, y) offset values
	 * 
	 * @param m
	 *            x difference
	 * @param n
	 *            y difference
	 * 
	 * @return Integer offset for looking up pixel in work array
	 */
	private int getOffset(int m, int n) {
		return m + n * width;
	}

	/**
	 * Check whole array replacing m with n
	 * 
	 * @param m
	 *            value to be replaced
	 * @param n
	 *            new value
	 * @param startZ
	 *            first z coordinate to check
	 * @param endZ
	 *            last+1 z coordinate to check
	 */
	private void replaceLabel(int[][] particleLabels, int m, int n,
			int startZ, int endZ) {
		for (int z = startZ; z < endZ; z++) {
			for (int i = 0; i < sliceSize; i++)
				if (particleLabels[z][i] == m) {
					particleLabels[z][i] = n;
				}
		}
	}

	/**
	 * Show a Results table containing some performance information
	 * 
	 * @param chunkRanges
	 * @param duration
	 */
	private void showResults(int[][] chunkRanges, double duration,
			ImagePlus imp, int slicesPerChunk) {
		ResultsTable rt = ResultsTable.getResultsTable();
		rt.incrementCounter();
		rt.addLabel(imp.getTitle());
		rt.addValue("Threads", nThreads);
		rt.addValue("Slices", nSlices);
		rt.addValue("Chunks", nChunks);
		rt.addValue("Chunk size", slicesPerChunk);
		rt.addValue("Last chunk size", chunkRanges[1][nChunks - 1]
				- chunkRanges[0][nChunks - 1]);
		rt.addValue("Duration (s)", duration);
		rt.show("Results");
		return;
	}
}
