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
		gd.showDialog();
		int slicesPerChunk = (int) Math.floor(gd.getNextNumber());
		if (gd.wasCanceled()) {
			return;
		}
		boolean showPerformance = gd.getNextBoolean();
		boolean doCopy = gd.getNextBoolean();
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
		}
		return;
	}

	/**
	 * 
	 * @param imp
	 * @param slicesPerChunk
	 * @param showPerformance
	 * @return
	 */
	public Object[] purify(ImagePlus imp, int slicesPerChunk,
			boolean showPerformance) {

		long startTime = System.currentTimeMillis();
		ParticleCounter pc = new ParticleCounter();

		final int fg = ParticleCounter.FORE;
		Object[] foregroundParticles = pc.getParticles(imp, slicesPerChunk, fg);
		byte[][] workArray = (byte[][]) foregroundParticles[0];
		int[][] particleLabels = (int[][]) foregroundParticles[1];
		long[] particleSizes = (long[]) foregroundParticles[2];
		removeSmallParticles(workArray, particleLabels, particleSizes, fg);
		
		final int bg = ParticleCounter.BACK;
		Object[] backgroundParticles = pc.getParticles(imp, workArray,
				slicesPerChunk, bg);
		particleLabels = (int[][]) backgroundParticles[1];
		particleSizes = (long[]) backgroundParticles[2];
		touchEdges(imp, workArray, particleLabels, particleSizes, bg);
		particleSizes = pc.getParticleSizes(particleLabels);
		removeSmallParticles(workArray, particleLabels, particleSizes, bg);

		double duration = ((double) System.currentTimeMillis() - (double) startTime)
				/ (double) 1000;

		if (showPerformance)
			showResults(duration, imp, slicesPerChunk);

		IJ.showStatus("Image Purified");

		ImageStack stack = new ImageStack(imp.getWidth(), imp.getHeight());
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
	private void touchEdges(ImagePlus imp, final byte[][] workArray,
			int[][] particleLabels, final long[] particleSizes, final int phase) {
		String status = "Background particles touching ";
		final int w = imp.getWidth();
		final int h = imp.getHeight();
		final int d = imp.getImageStackSize();
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

		ParticleCounter pc = new ParticleCounter();
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
					pc.replaceLabel(particleLabels, particleLabels[z][offset],
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
					pc.replaceLabel(particleLabels, particleLabels[z][offset],
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
					pc.replaceLabel(particleLabels, particleLabels[z][offset],
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
					pc.replaceLabel(particleLabels, particleLabels[z][offset],
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
					pc.replaceLabel(particleLabels, particleLabels[z][offset],
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
					pc.replaceLabel(particleLabels, particleLabels[z][offset],
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
		final int d = workArray.length;
		final int wh = workArray[0].length;
		final int fg = ParticleCounter.FORE;
		final int bg = ParticleCounter.BACK;
		long maxVC = 0;
		final int nPartSizes = particleSizes.length;
		for (int i = 0; i < nPartSizes; i++) {
			if (particleSizes[i] > maxVC) {
				maxVC = particleSizes[i];
			}
		}
		final long maxVoxCount = maxVC;
		if (phase == fg) {
			// go through work array and turn all
			// smaller foreground particles into background (0)
			for (int z = 0; z < d; z++) {
				for (int i = 0; i < wh; i++) {
					if (workArray[z][i] == fg) {
						if (particleSizes[particleLabels[z][i]] < maxVoxCount) {
							workArray[z][i] = bg;
						}
					}
				}
				IJ.showStatus("Removing foreground particles");
				IJ.showProgress(z, d);
			}
		} else if (phase == bg) {
			// go through work array and turn all
			// smaller background particles into foreground
			for (int z = 0; z < d; z++) {
				for (int i = 0; i < wh; i++) {
					if (workArray[z][i] == bg) {
						if (particleSizes[particleLabels[z][i]] < maxVoxCount) {
							workArray[z][i] = fg;
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
	 * Show a Results table containing some performance information
	 * 
	 * @param chunkRanges
	 * @param duration
	 */
	private void showResults(double duration, ImagePlus imp, int slicesPerChunk) {
		ParticleCounter pc = new ParticleCounter();
		final int nChunks = pc.getNChunks(imp, slicesPerChunk);
		int[][] chunkRanges = pc.getChunkRanges(imp, nChunks, slicesPerChunk);
		ResultsTable rt = ResultsTable.getResultsTable();
		rt.incrementCounter();
		rt.addLabel(imp.getTitle());
		rt.addValue("Threads", Runtime.getRuntime().availableProcessors());
		rt.addValue("Slices", imp.getImageStackSize());
		rt.addValue("Chunks", nChunks);
		rt.addValue("Chunk size", slicesPerChunk);
		rt.addValue("Last chunk size", chunkRanges[1][nChunks - 1]
				- chunkRanges[0][nChunks - 1]);
		rt.addValue("Duration (s)", duration);
		rt.show("Results");
		return;
	}
}
