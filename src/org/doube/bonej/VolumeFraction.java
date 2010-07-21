package org.doube.bonej;

/**
 * VolumeFraction plugin for ImageJ
 * Copyright 2009 2010 Michael Doube
 * 
 *This program is free software: you can redistribute it and/or modify
 *it under the terms of the GNU General Public License as published by
 *the Free Software Foundation, either version 3 of the License, or
 *(at your option) any later version.
 *
 *This program is distributed in the hope that it will be useful,
 *but WITHOUT ANY WARRANTY; without even the implied warranty of
 *MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *GNU General Public License for more details.
 *
 *You should have received a copy of the GNU General Public License
 *along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

import java.awt.AWTEvent;
import java.awt.Choice;
import java.awt.Rectangle;
import java.awt.TextField;
import java.util.List;
import java.util.Vector;

import javax.vecmath.Color3f;
import javax.vecmath.Point3f;

import marchingcubes.MCTriangulator;

import org.doube.util.DialogModifier;
import org.doube.util.ImageCheck;
import org.doube.util.ResultInserter;

import customnode.CustomTriangleMesh;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.measure.Calibration;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;
import ij.plugin.PlugIn;

import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.gui.WaitForUserDialog;

public class VolumeFraction implements PlugIn, DialogListener {

	public void run(String arg) {
		if (!ImageCheck.checkEnvironment())
			return;
		final ImagePlus imp = IJ.getImage();
		if (null == imp) {
			IJ.noImage();
			return;
		}
		if (imp.getBitDepth() == 32 || imp.getBitDepth() == 24) {
			IJ
					.error("Volume Fraction requires a binary, 8-bit or 16-bit image");
			return;
		}

		GenericDialog gd = new GenericDialog("Volume");
		String[] types = { "Voxel", "Surface" };
		gd.addChoice("Algorithm", types, types[0]);
		gd.addNumericField("Surface resampling", 6, 0);
		gd.addHelp("http://bonej.org/volumefraction");
		gd.addDialogListener(this);
		gd.showDialog();
		if (gd.wasCanceled()) {
			return;
		}
		String type = gd.getNextChoice();
		final int resampling = (int) Math.floor(gd.getNextNumber());

		final double[] thresholds = setThreshold(imp);
		final double minT = thresholds[0];
		final double maxT = thresholds[1];

		double[] volumes = new double[2];
		if (type.equals(types[0])) {
			volumes = getVolumes(imp, minT, maxT);
		} else if (type.equals(types[1])) {
			try {
				volumes = getSurfaceVolume(imp, minT, maxT, resampling);
			} catch (Exception e) {
				IJ.log(e.getMessage());
			}
		}
		double volBone = volumes[0];
		double volTotal = volumes[1];
		double p = volBone / volTotal;
		Calibration cal = imp.getCalibration();

		ResultInserter ri = ResultInserter.getInstance();
		ri.setResultInRow(imp, "BV (" + cal.getUnits() + "³)", volBone);
		ri.setResultInRow(imp, "TV (" + cal.getUnits() + "³)", volTotal);
		ri.setResultInRow(imp, "BV/TV", p);
		ri.updateTable();
		return;
	}

	/**
	 * Get the total and thresholded volumes of a masked area
	 * 
	 * @param imp
	 *            Image
	 * @param minT
	 *            minimum threshold (inclusive)
	 * @param maxT
	 *            maximum threshold (inclusive)
	 * @return double[2] containing the foreground and total volumes
	 * 
	 */
	public double[] getVolumes(ImagePlus imp, double minT, double maxT) {
		final ImageStack stack = imp.getImageStack();
		final ImageProcessor mask = imp.getProcessor().getMask();
		final Rectangle r = imp.getProcessor().getRoi();

		int nThreads = Runtime.getRuntime().availableProcessors();
		long[] volTotalT = new long[nThreads];
		long[] volBoneT = new long[nThreads];
		VolumesThread[] vt = new VolumesThread[nThreads];
		for (int thread = 0; thread < nThreads; thread++) {
			vt[thread] = new VolumesThread(thread, nThreads, stack, minT, maxT,
					mask, r, volTotalT, volBoneT);
			vt[thread].start();
		}
		try {
			for (int thread = 0; thread < nThreads; thread++) {
				vt[thread].join();
			}
		} catch (InterruptedException ie) {
			IJ.error("A thread was interrupted.");
		}

		long volTotal = 0;
		long volBone = 0;
		for (int i = 0; i < nThreads; i++){
			volTotal += volTotalT[i];
			volBone += volBoneT[i];
		}
		Calibration cal = imp.getCalibration();
		double voxelVol = cal.pixelWidth * cal.pixelHeight * cal.pixelDepth;
		double[] volumes = { volBone * voxelVol, volTotal * voxelVol };
		return volumes;
	}

	class VolumesThread extends Thread {
		int thread, nThreads;
		ImageStack stack;
		double minT, maxT;
		long[] volTotalT, volBoneT;
		ImageProcessor mask;
		Rectangle r;

		public VolumesThread(int thread, int nThreads, ImageStack stack,
				double minT, double maxT, ImageProcessor mask, Rectangle r,
				long[] volTotalT, long[] volBoneT) {
			this.thread = thread;
			this.nThreads = nThreads;
			this.stack = stack;
			this.minT = minT;
			this.maxT = maxT;
			this.volTotalT = volTotalT;
			this.volBoneT = volBoneT;
			this.mask = mask;
			this.r = r;
		}

		public void run() {
			final int nSlices = stack.getSize();
			final int rLeft = r.x;
			final int rTop = r.y;
			final int rRight = rLeft + r.width;
			final int rBottom = rTop + r.height;
			final boolean hasMask = (mask != null);
			for (int s = thread + 1; s <= nSlices; s += nThreads) {
				ImageProcessor ipSlice = stack.getProcessor(s);
				for (int v = rTop; v < rBottom; v++) {
					final int vrTop = v - rTop;
					for (int u = rLeft; u < rRight; u++) {
						if (!hasMask || mask.get(u - rLeft, vrTop) > 0) {
							volTotalT[thread]++;
							final double pixel = ipSlice.get(u, v);
							if (pixel >= minT && pixel <= maxT) {
								volBoneT[thread]++;
							}
						}
					}
				}
			}
		}
	}

	/**
	 * Calculate the foreground (bone) and total volumes, using surface meshes.
	 * 
	 * @param imp
	 *            Input ImagePlus
	 * @param minT
	 *            threshold minimum
	 * @param maxT
	 *            threshold maximum
	 * @param resampling
	 *            voxel resampling for mesh creation; higher values result in
	 *            simpler meshes
	 * @return
	 */
	@SuppressWarnings("unchecked")
	public double[] getSurfaceVolume(ImagePlus imp, double minT, double maxT,
			int resampling) {
		ImageProcessor ip = imp.getProcessor();
		final ImageStack stack = imp.getImageStack();
		final ImageProcessor mask = ip.getMask();
		final boolean hasMask = (mask != null);
		Rectangle r = ip.getRoi();
		final int rLeft = r.x;
		final int rTop = r.y;
		final int rRight = rLeft + r.width;
		final int rBottom = rTop + r.height;
		final int nSlices = imp.getStackSize();
		ImageStack outStack = new ImageStack(r.width, r.height);
		ImageStack maskStack = new ImageStack(r.width, r.height);
		for (int s = 1; s <= nSlices; s++) {
			ImageProcessor ipSlice = stack.getProcessor(s);
			ByteProcessor ipOut = new ByteProcessor(r.width, r.height);
			ByteProcessor ipMask = new ByteProcessor(r.width, r.height);
			for (int v = rTop; v < rBottom; v++) {
				final int vrTop = v - rTop;
				for (int u = rLeft; u < rRight; u++) {
					if (!hasMask || mask.get(u - rLeft, vrTop) > 0) {
						ipMask.set(u - rLeft, v - rTop, (byte) 255);
						final double pixel = ipSlice.get(u, v);
						if (pixel >= minT && pixel <= maxT) {
							ipOut.set(u - rLeft, v - rTop, (byte) 255);
						} else {
							ipOut.set(u - rLeft, v - rTop, (byte) 0);
						}
					}
				}
			}
			outStack.addSlice(stack.getSliceLabel(s), ipOut);
			maskStack.addSlice(stack.getSliceLabel(s), ipMask);
		}
		ImagePlus impOut = new ImagePlus();
		impOut.setStack("Out", outStack);
		impOut.setCalibration(imp.getCalibration());
		final Color3f colour = new Color3f(0.0f, 0.0f, 0.0f);
		boolean[] channels = { true, false, false };
		MCTriangulator mct = new MCTriangulator();
		List<Point3f> points = mct.getTriangles(impOut, 128, channels,
				resampling);
		CustomTriangleMesh surface = new CustomTriangleMesh(points, colour,
				0.0f);
		double boneVolume = surface.getVolume();
		ImagePlus maskImp = new ImagePlus("Mask", maskStack);
		maskImp.setCalibration(imp.getCalibration());
		points = mct.getTriangles(maskImp, 128, channels, resampling);
		surface = new CustomTriangleMesh(points, colour, 0.0f);
		double totalVolume = surface.getVolume();
		double[] volumes = { boneVolume, totalVolume };

		return volumes;
	}

	private double[] setThreshold(ImagePlus imp) {
		double[] thresholds = new double[2];
		ImageCheck ic = new ImageCheck();
		if (ic.isBinary(imp)) {
			thresholds[0] = 128;
			thresholds[1] = 255;
		} else {
			IJ.run("Threshold...");
			new WaitForUserDialog("Set the threshold, then click OK.").show();
			thresholds[0] = imp.getProcessor().getMinThreshold();
			thresholds[1] = imp.getProcessor().getMaxThreshold();
		}
		return thresholds;
	}

	public boolean dialogItemChanged(GenericDialog gd, AWTEvent e) {
		Vector<?> choices = gd.getChoices();
		Choice choice = (Choice) choices.get(0);
		Vector<?> numbers = gd.getNumericFields();
		TextField num = (TextField) numbers.get(0);
		if (choice.getSelectedIndex() == 1) {
			num.setEnabled(true);
		} else {
			num.setEnabled(false);
		}
		DialogModifier.registerMacroValues(gd, gd.getComponents());
		return true;
	}
}
