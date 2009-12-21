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

	public final static int FORE = -1, BACK = 0;

	private int width, height;
	private int nSlices, nThreads;
	private int nChunks, sliceSize;

	private int ID;

	private String chunkString = "";

	private String sPhase = "";

	public void run(String arg) {
		if (!ImageCheck.checkIJVersion())
			return;
		ImagePlus imp = IJ.getImage();
		ImageCheck ic = new ImageCheck();
		if (!ic.isBinary(imp)) {
			IJ.error("Purify requires a binary image");
			return;
		}
		GenericDialog gd = new GenericDialog("Setup");
		gd.addNumericField("Chunk Size", 4, 0, 4, "slices");
		gd.addCheckbox("Performance Log", false);
		gd.addCheckbox("Make_copy", true);
		gd.addCheckbox("Show_particles", false);
		gd.showDialog();
		int slicesPerChunk = (int) Math.floor(gd.getNextNumber());
		if (gd.wasCanceled()) {
			return;
		}
		boolean showPerformance = gd.getNextBoolean();
		boolean doCopy = gd.getNextBoolean();
		boolean doParticleImage = gd.getNextBoolean();
		Object[] result = purify(imp, slicesPerChunk, showPerformance);
		if (null != result) {
			ImagePlus purified = (ImagePlus) result[1];

			if (doCopy) {
				purified.show();
				if (!purified.isInvertedLut())
					IJ.run("Invert LUT");
			} else {
				ImageStack stack2 = purified.getStack();
				imp.setStack(null, stack2);
				if (!imp.isInvertedLut())
					IJ.run("Invert LUT");
			}
			if (doParticleImage) {
				displayParticleLabels((int[][]) result[2], imp).show();
			}
		}
		return;
	}

	public Object[] purify(ImagePlus imp, int slicesPerChunk,
			boolean showPerformance) {

		long startTime = System.currentTimeMillis();

		sPhase = "foreground";
		Object[] foregroundParticles = getParticles(imp, slicesPerChunk, FORE);
		byte[][] workArray = (byte[][]) foregroundParticles[0];
		int[][] particleLabels = (int[][]) foregroundParticles[1];
		long[] particleSizes = (long[]) foregroundParticles[2];
		removeSmallParticles(workArray, particleLabels, particleSizes, FORE);

		sPhase = "background";
		Object[] backgroundParticles = getParticles(imp, workArray,
				slicesPerChunk, BACK);
		particleLabels = (int[][]) backgroundParticles[1];
		particleSizes = (long[]) backgroundParticles[2];
		touchEdges(workArray, particleLabels, particleSizes, BACK);
		particleSizes = getParticleSizes(workArray, particleLabels, BACK);
		removeSmallParticles(workArray, particleLabels, particleSizes, BACK);

		double duration = ((double) System.currentTimeMillis() - (double) startTime)
				/ (double) 1000;
		if (showPerformance)
			showResults(duration, imp, slicesPerChunk);

		IJ.showStatus("Image Purified");

		ImageStack stack = new ImageStack(this.width, this.height);
		final int nSlices = workArray.length;
		for (int z = 0; z < nSlices; z++) {
			stack.addSlice(imp.getStack().getSliceLabel(z + 1), workArray[z]);
		}
		ImagePlus purified = new ImagePlus("Purified", stack);
		purified.setCalibration(imp.getCalibration());
		Object[] result = { duration, purified, particleLabels };
		return result;
	}

	/**
	 * Create a work array
	 * 
	 * @return byte[] work array
	 */
	private byte[][] makeWorkArray(ImagePlus imp) {
		final int s = imp.getStackSize();
		final int p = imp.getWidth() * imp.getHeight();
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
	private int[][] firstIDAttribution(final byte[][] workArray, final int phase) {
		final int w = width;
		final int h = height;
		final int d = nSlices;
		IJ.showStatus("Finding " + sPhase + " structures");
		int[][] particleLabels = new int[nSlices][sliceSize];
		ID = 1;

		if (phase == FORE) {
			for (int z = 0; z < d; z++) {
				for (int y = 0; y < h; y++) {
					final int rowIndex = y * w;
					for (int x = 0; x < w; x++) {
						final int arrayIndex = rowIndex + x;
						if (workArray[z][arrayIndex] == phase) {
							particleLabels[z][arrayIndex] = ID;
							int minTag = ID;
							// Find the minimum particleLabel in the
							// neighbouring pixels
							for (int vZ = z - 1; vZ <= z + 1; vZ++) {
								for (int vY = y - 1; vY <= y + 1; vY++) {
									for (int vX = x - 1; vX <= x + 1; vX++) {
										if (withinBounds(vX, vY, vZ, w, h, 0, d)) {
											final int offset = getOffset(vX,
													vY, w);
											if (workArray[vZ][offset] == phase) {
												final int tagv = particleLabels[vZ][offset];
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
				IJ.showProgress(z, d);
			}
			ID++;
		} else if (phase == BACK) {
			for (int z = 0; z < d; z++) {
				for (int y = 0; y < h; y++) {
					final int rowIndex = y * w;
					for (int x = 0; x < w; x++) {
						final int arrayIndex = rowIndex + x;
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
								if (withinBounds(nX, nY, nZ, w, h, 0, d)) {
									final int offset = getOffset(nX, nY, w);
									if (workArray[nZ][offset] == phase) {
										final int tagv = particleLabels[nZ][offset];
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
				IJ.showProgress(z, d);
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
	private void connectStructures(final byte[][] workArray,
			int[][] particleLabels, final int phase, final int[][] scanRanges) {
		IJ.showStatus("Connecting " + sPhase + " structures" + chunkString);
		final int w = width;
		final int h = height;
		final int d = nSlices;
		for (int c = 0; c < scanRanges[0].length; c++) {
			final int sR0 = scanRanges[0][c];
			final int sR1 = scanRanges[1][c];
			final int sR2 = scanRanges[2][c];
			final int sR3 = scanRanges[3][c];
			if (phase == FORE) {
				for (int z = sR0; z < sR1; z++) {
					for (int y = 0; y < h; y++) {
						final int rowIndex = y * w;
						for (int x = 0; x < w; x++) {
							final int arrayIndex = rowIndex + x;
							if (workArray[z][arrayIndex] == phase
									&& particleLabels[z][arrayIndex] > 1) {
								int minTag = particleLabels[z][arrayIndex];
								// Find the minimum particleLabel in the
								// neighbours' pixels
								for (int vZ = z - 1; vZ <= z + 1; vZ++) {
									for (int vY = y - 1; vY <= y + 1; vY++) {
										for (int vX = x - 1; vX <= x + 1; vX++) {
											if (withinBounds(vX, vY, vZ, w, h,
													sR2, sR3)) {
												final int offset = getOffset(
														vX, vY, w);
												if (workArray[vZ][offset] == phase) {
													final int tagv = particleLabels[vZ][offset];
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
											if (withinBounds(vX, vY, vZ, w, h,
													sR2, sR3)) {
												final int offset = getOffset(
														vX, vY, w);
												if (workArray[vZ][offset] == phase) {
													final int tagv = particleLabels[vZ][offset];
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
					IJ.showProgress(z, d);
				}
			} else if (phase == BACK) {
				for (int z = sR0; z < sR1; z++) {
					for (int y = 0; y < h; y++) {
						final int rowIndex = y * w;
						for (int x = 0; x < w; x++) {
							final int arrayIndex = rowIndex + x;
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
									if (withinBounds(nX, nY, nZ, w, h, sR2, sR3)) {
										final int offset = getOffset(nX, nY, w);
										if (workArray[nZ][offset] == phase) {
											final int tagv = particleLabels[nZ][offset];
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
									if (withinBounds(nX, nY, nZ, w, h, sR2, sR3)) {
										final int offset = getOffset(nX, nY, w);
										if (workArray[nZ][offset] == phase) {
											final int tagv = particleLabels[nZ][offset];
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
					IJ.showProgress(z, d + 1);
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
				byte[][] workArray, int[][] particleLabels, final int phase,
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
	private long[] getParticleSizes(final byte[][] workArray,
			final int[][] particleLabels, final int phase) {
		IJ.showStatus("Getting " + sPhase + " particle sizes");
		final int w = width;
		final int h = height;
		final int d = nSlices;
		final int wh = w * h;
		long[] particleSizes = new long[ID];
		for (int z = 0; z < d; z++) {
			for (int arrayIndex = 0; arrayIndex < wh; arrayIndex++) {
				if (workArray[z][arrayIndex] == phase) {
					particleSizes[particleLabels[z][arrayIndex]]++;
				}
			}
			IJ.showProgress(z, d);
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
	private void touchEdges(final byte[][] workArray, int[][] particleLabels,
			final long[] particleSizes, final int phase) {
		String status = "Background particles touching ";
		final int w = width;
		final int h = height;
		final int d = nSlices;
		// find the label associated with the biggest
		// particle in phase
		long maxVoxCount = 0;
		int bigP = 0;
		final int nPartSizes = particleSizes.length;
		for (int i = 0; i < nPartSizes; i++) {
			if (particleSizes[i] > maxVoxCount) {
				maxVoxCount = particleSizes[i];
				bigP = i;
			}
		}
		final int biggestParticle = bigP;
		// check each face of the stack for pixels that are touching edges and
		// replace that particle's label in particleLabels with
		// the label of the biggest particle
		int x, y, z;

		// up
		z = 0;
		for (y = 0; y < h; y++) {
			IJ.showStatus(status + "top");
			IJ.showProgress(y, h);
			final int rowOffset = y * w;
			for (x = 0; x < w; x++) {
				final int offset = rowOffset + x;
				if (workArray[z][offset] == phase
						&& particleLabels[z][offset] != biggestParticle) {
					replaceLabel(particleLabels, particleLabels[z][offset],
							biggestParticle, 0, d);
				}
			}
		}

		// down
		z = d - 1;
		for (y = 0; y < h; y++) {
			IJ.showStatus(status + "bottom");
			IJ.showProgress(y, h);
			final int rowOffset = y * w;
			for (x = 0; x < w; x++) {
				final int offset = rowOffset + x;
				if (workArray[z][offset] == phase
						&& particleLabels[z][offset] != biggestParticle) {
					replaceLabel(particleLabels, particleLabels[z][offset],
							biggestParticle, 0, d);
				}
			}
		}

		// left
		x = 0;
		for (z = 0; z < d; z++) {
			IJ.showStatus(status + "left");
			IJ.showProgress(z, d);
			for (y = 0; y < h; y++) {
				final int offset = y * w;
				if (workArray[z][offset] == phase
						&& particleLabels[z][offset] != biggestParticle) {
					replaceLabel(particleLabels, particleLabels[z][offset],
							biggestParticle, 0, d);
				}
			}
		}

		// right
		x = w - 1;
		for (z = 0; z < d; z++) {
			IJ.showStatus(status + "right");
			IJ.showProgress(z, d);
			for (y = 0; y < h; y++) {
				final int offset = y * w + x;
				if (workArray[z][offset] == phase
						&& particleLabels[z][offset] != biggestParticle) {
					replaceLabel(particleLabels, particleLabels[z][offset],
							biggestParticle, 0, d);
				}
			}
		}

		// front
		y = h - 1;
		int rowOffset = y * w;
		for (z = 0; z < d; z++) {
			IJ.showStatus(status + "front");
			IJ.showProgress(z, d);
			for (x = 0; x < w; x++) {
				final int offset = rowOffset + x;
				if (workArray[z][offset] == phase
						&& particleLabels[z][offset] != biggestParticle) {
					replaceLabel(particleLabels, particleLabels[z][offset],
							biggestParticle, 0, d);
				}
			}
		}

		// back
		y = 0;
		for (z = 0; z < d; z++) {
			IJ.showStatus(status + "back");
			IJ.showProgress(z, d);
			for (x = 0; x < w; x++) {
				final int offset = x;
				if (workArray[z][offset] == phase
						&& particleLabels[z][offset] != biggestParticle) {
					replaceLabel(particleLabels, particleLabels[z][offset],
							biggestParticle, 0, d);
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
			final int[][] particleLabels, final long[] particleSizes,
			final int phase) {
		final int d = nSlices;
		final int s = sliceSize;
		long maxVC = 0;
		final int nPartSizes = particleSizes.length;
		for (int i = 0; i < nPartSizes; i++) {
			if (particleSizes[i] > maxVC) {
				maxVC = particleSizes[i];
			}
		}
		final long maxVoxCount = maxVC;
		if (phase == FORE) {
			// go through work array and turn all
			// smaller foreground particles into background (0)
			for (int z = 0; z < d; z++) {
				for (int i = 0; i < s; i++) {
					if (workArray[z][i] == FORE) {
						if (particleSizes[particleLabels[z][i]] < maxVoxCount) {
							workArray[z][i] = BACK;
						}
					}
				}
				IJ.showStatus("Removing foreground particles");
				IJ.showProgress(z, d);
			}
		} else if (phase == BACK) {
			// go through work array and turn all
			// smaller background particles into foreground
			for (int z = 0; z < d; z++) {
				for (int i = 0; i < s; i++) {
					if (workArray[z][i] == BACK) {
						if (particleSizes[particleLabels[z][i]] < maxVoxCount) {
							workArray[z][i] = FORE;
						}
					}
				}
				IJ.showStatus("Removing background particles");
				IJ.showProgress(z, d);
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
	private boolean withinBounds(int m, int n, int o, int w, int h, int startZ,
			int endZ) {
		return (m >= 0 && m < w && n >= 0 && n < h && o >= startZ && o < endZ);
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
	private int getOffset(int m, int n, int w) {
		return m + n * w;
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
	private void replaceLabel(int[][] particleLabels, final int m, int n,
			int startZ, final int endZ) {
		final int s = sliceSize;
		for (int z = startZ; z < endZ; z++) {
			for (int i = 0; i < s; i++)
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
	private void showResults(double duration, ImagePlus imp, int slicesPerChunk) {
		int[][] chunkRanges = getChunkRanges(nChunks, slicesPerChunk);
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

	/**
	 * Display the particle labels as an ImagePlus
	 * 
	 * @param particleLabels
	 * @param imp
	 *            original image, used for image dimensions, calibration and
	 *            titles
	 */
	private ImagePlus displayParticleLabels(int[][] particleLabels,
			ImagePlus imp) {
		final int w = imp.getWidth();
		final int h = imp.getHeight();
		final int d = imp.getImageStackSize();
		ImageStack stackParticles = new ImageStack(w, h);
		for (int z = 1; z <= d; z++) {
			stackParticles.addSlice(imp.getImageStack().getSliceLabel(z),
					particleLabels[z - 1]);
		}
		ImagePlus impParticles = new ImagePlus(imp.getShortTitle() + "_parts",
				stackParticles);
		impParticles.setCalibration(imp.getCalibration());
		return impParticles;
	}

	/**
	 * Get particles, particle labels and particle sizes from a 3D ImagePlus
	 * 
	 * @param imp
	 *            Binary input image
	 * @param slicesPerChunk
	 *            number of slices per chunk. 2 is generally good.
	 * @param phase
	 *            foreground or background (FORE or BACK)
	 * @return Object[] {byte[][], int[][], long[]} containing a binary
	 *         workArray, particle labels and particle sizes indexed by particle
	 *         label
	 */
	public Object[] getParticles(ImagePlus imp, int slicesPerChunk, int phase) {
		byte[][] workArray = makeWorkArray(imp);
		return getParticles(imp, workArray, slicesPerChunk, phase);
	}

	/**
	 * Get particles, particle labels and sizes from a workArray
	 * 
	 * @param imp
	 *            input binary image
	 * @param slicesPerChunk
	 *            number of slices to use for each chunk
	 * @param phase
	 *            one of either Purify.foreground or Purify.background
	 * @return Object[] array containing a binary workArray, particle labels and
	 *         particle sizes
	 */
	private Object[] getParticles(ImagePlus imp, byte[][] workArray,
			int slicesPerChunk, int phase) {
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

		int[][] particleLabels = firstIDAttribution(workArray, phase);

		// connect foreground particles within chunks
		ConnectStructuresThread[] cptf = new ConnectStructuresThread[nThreads];
		for (int thread = 0; thread < nThreads; thread++) {
			cptf[thread] = new ConnectStructuresThread(thread, nThreads,
					workArray, particleLabels, phase, nChunks, chunkRanges);
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
			connectStructures(workArray, particleLabels, phase, stitchRanges);
		}

		long[] particleSizes = getParticleSizes(workArray, particleLabels,
				phase);

		Object[] result = { workArray, particleLabels, particleSizes };
		return result;
	}
}
