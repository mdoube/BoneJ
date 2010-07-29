package org.doube.bonej;

/**
 *  
 * Purify plugin for ImageJ
 * Copyright 2009 2010 Michael Doube 
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

import java.awt.AWTEvent;
import java.awt.Choice;
import java.awt.TextField;
import java.util.Vector;
import java.util.concurrent.atomic.AtomicInteger;

import org.doube.util.DialogModifier;
import org.doube.util.ImageCheck;
import org.doube.util.Multithreader;

import ij.*;
import ij.plugin.PlugIn;
import ij.measure.ResultsTable;
import ij.gui.DialogListener;
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
 * 
 */
public class Purify implements PlugIn, DialogListener {

	public void run(String arg) {
		if (!ImageCheck.checkEnvironment())
			return;
		ImagePlus imp = IJ.getImage();
		ImageCheck ic = new ImageCheck();
		if (!ic.isBinary(imp)) {
			IJ.error("Purify requires a binary image");
			return;
		}
		GenericDialog gd = new GenericDialog("Setup");
		String[] items = { "Multithreaded", "Linear" };
		gd.addChoice("Labelling algorithm", items, items[0]);
		gd.addNumericField("Chunk Size", 4, 0, 4, "slices");
		gd.addCheckbox("Performance Log", false);
		gd.addCheckbox("Make_copy", true);
		gd.addHelp("http://bonej.org/purify");
		gd.addDialogListener(this);
		gd.showDialog();
		if (gd.wasCanceled())
			return;
		final String choice = gd.getNextChoice();
		int labelMethod;
		if (choice.equals(items[0])) {
			labelMethod = ParticleCounter.MULTI;
		} else {
			labelMethod = ParticleCounter.LINEAR;
		}
		int slicesPerChunk = (int) Math.floor(gd.getNextNumber());
		boolean showPerformance = gd.getNextBoolean();
		boolean doCopy = gd.getNextBoolean();
		long startTime = System.currentTimeMillis();
		ImagePlus purified = purify(imp, slicesPerChunk, labelMethod);
		if (null != purified) {
			if (doCopy) {
				purified.show();
				if (imp.isInvertedLut() && !purified.isInvertedLut())
					IJ.run("Invert LUT");
			} else {
				imp.setStack(null, purified.getStack());
				if (!imp.isInvertedLut())
					IJ.run("Invert LUT");
			}
		}
		double duration = ((double) System.currentTimeMillis() - (double) startTime)
				/ (double) 1000;

		if (showPerformance)
			showResults(duration, imp, slicesPerChunk, labelMethod);
		return;
	}

	/**
	 * Find all foreground and particles in an image and remove all but the
	 * largest. Foreground is 26-connected and background is 8-connected.
	 * 
	 * @param imp
	 *            input image
	 * @param slicesPerChunk
	 *            number of slices to send to each CPU core as a chunk
	 * @return purified image
	 */
	public ImagePlus purify(ImagePlus imp, int slicesPerChunk, int labelMethod) {

		ParticleCounter pc = new ParticleCounter();
		pc.setLabelMethod(labelMethod);

		final int fg = ParticleCounter.FORE;
		Object[] foregroundParticles = pc.getParticles(imp, slicesPerChunk, 0,
				Double.POSITIVE_INFINITY, fg);
		byte[][] workArray = (byte[][]) foregroundParticles[0];
		int[][] particleLabels = (int[][]) foregroundParticles[1];
		// index 0 is background particle's size...
		long[] particleSizes = pc.getParticleSizes(particleLabels);
		removeSmallParticles(workArray, particleLabels, particleSizes, fg);

		final int bg = ParticleCounter.BACK;
		Object[] backgroundParticles = pc.getParticles(imp, workArray,
				slicesPerChunk, 0, Double.POSITIVE_INFINITY, bg);
		particleLabels = (int[][]) backgroundParticles[1];
		particleSizes = pc.getParticleSizes(particleLabels);
		touchEdges(imp, workArray, particleLabels, particleSizes, bg);
		particleSizes = pc.getParticleSizes(particleLabels);
		removeSmallParticles(workArray, particleLabels, particleSizes, bg);

		IJ.showStatus("Image Purified");

		ImageStack stack = new ImageStack(imp.getWidth(), imp.getHeight());
		final int nSlices = workArray.length;
		for (int z = 0; z < nSlices; z++) {
			stack.addSlice(imp.getStack().getSliceLabel(z + 1), workArray[z]);
		}
		ImagePlus purified = new ImagePlus("Purified", stack);
		purified.setCalibration(imp.getCalibration());
		return purified;
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
							biggestParticle, 0, d, true);
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
							biggestParticle, 0, d, true);
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
							biggestParticle, 0, d, true);
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
							biggestParticle, 0, d, true);
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
							biggestParticle, 0, d, true);
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
							biggestParticle, 0, d, true);
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
	private void removeSmallParticles(final byte[][] workArray,
			final int[][] particleLabels, final long[] particleSizes,
			final int phase) {
		final int d = workArray.length;
		final int wh = workArray[0].length;
		final int fg = ParticleCounter.FORE;
		final int bg = ParticleCounter.BACK;
		long maxVC = 0;
		final int nPartSizes = particleSizes.length;
		for (int i = 1; i < nPartSizes; i++) {
			if (particleSizes[i] > maxVC) {
				maxVC = particleSizes[i];
			}
		}
		final long maxVoxCount = maxVC;
		final AtomicInteger ai = new AtomicInteger(0);
		Thread[] threads = Multithreader.newThreads();
		for (int thread = 0; thread < threads.length; thread++) {
			threads[thread] = new Thread(new Runnable() {
				public void run() {
					if (phase == fg) {
						// go through work array and turn all
						// smaller foreground particles into background (0)
						for (int z = ai.getAndIncrement(); z < d; z = ai
								.getAndIncrement()) {
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
						for (int z = ai.getAndIncrement(); z < d; z = ai
								.getAndIncrement()) {
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
				}
			});
		}
		Multithreader.startAndJoin(threads);
		return;
	}

	/**
	 * Show a Results table containing some performance information
	 * 
	 * @param chunkRanges
	 * @param duration
	 */
	private void showResults(double duration, ImagePlus imp,
			int slicesPerChunk, int labelMethod) {
		if (labelMethod == ParticleCounter.LINEAR)
			slicesPerChunk = imp.getImageStackSize();
		ParticleCounter pc = new ParticleCounter();
		final int nChunks = pc.getNChunks(imp, slicesPerChunk);
		int[][] chunkRanges = pc.getChunkRanges(imp, nChunks, slicesPerChunk);
		ResultsTable rt = ResultsTable.getResultsTable();
		rt.incrementCounter();
		rt.addLabel(imp.getTitle());
		rt.addValue("Algorithm", labelMethod);
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

	@Override
	public boolean dialogItemChanged(GenericDialog gd, AWTEvent e) {
		Vector<?> choices = gd.getChoices();
		Vector<?> numbers = gd.getNumericFields();
		Choice choice = (Choice) choices.get(0);
		TextField num = (TextField) numbers.get(0);
		if (choice.getSelectedItem().contentEquals("Multithreaded")) {
			num.setEnabled(true);
		} else {
			num.setEnabled(false);
		}
		DialogModifier.registerMacroValues(gd, gd.getComponents());
		return true;
	}
}
