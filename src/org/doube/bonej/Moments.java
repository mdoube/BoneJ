package org.doube.bonej;

/** Moments ImageJ plugin Copyright 2008 2009 2010 Michael Doube 
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
 *
 */

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ImageProcessor;
import ij.process.StackConverter;
import ij.plugin.Duplicator;
import ij.plugin.PlugIn;
import ij.measure.Calibration;
import ij.gui.*;
import ij3d.Content;
import ij3d.Image3DUniverse;

import java.awt.AWTEvent;
//import java.awt.Checkbox;
import java.awt.Choice;
import java.awt.Rectangle;
import java.awt.TextField;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Vector;

import javax.vecmath.Color3f;
import javax.vecmath.Point3f;

import org.doube.jama.Matrix;
import org.doube.jama.EigenvalueDecomposition;
import org.doube.util.DialogModifier;
import org.doube.util.ImageCheck;
import org.doube.util.ResultInserter;
import org.doube.util.ThresholdGuesser;

import customnode.CustomPointMesh;

/**
 * Calculate centroid and principal axes of a thresholded stack; originally
 * designed for 16-bit CT scans of a bone in air so default thresholds are 0 and
 * 4000 HU, but most greyscale images should be handled
 * 
 * @author Michael Doube
 * 
 */
public class Moments implements PlugIn, DialogListener {

	private boolean fieldUpdated = false;
	private Calibration cal;
	private String pixUnits = "grey";
	/** Use pixel values for density */
	private static final int PIXEL = 0;
	/** Use calibration table for density */
	private static final int TABLE = 1;
	/** Calculate density from pixel value and user-defined linear function */
	private static final int LINEAR_PIXEL = 2;
	/** Density from calibration table and user-defined linear function */
	private static final int LINEAR_CAL = 3;
	/** Don't do any density calibration */
	private static final int NONE = 4;
	private int calMode = PIXEL;
	
	public void run(String arg) {
		if (!ImageCheck.checkEnvironment())
			return;
		final ImagePlus imp = IJ.getImage();
		if (null == imp) {
			IJ.noImage();
			return;
		}
		cal = imp.getCalibration();
		double[] thresholds = ThresholdGuesser.setDefaultThreshold(imp);
		double min = thresholds[0];
		double max = thresholds[1];
		// if (ImageCheck.huCalibrated(imp)) {
		// pixUnits = "HU";
		// fieldUpdated = true;
		// } else
		// pixUnits = "grey";
		GenericDialog gd = new GenericDialog("Setup");
		gd.addNumericField("Start Slice:", 1, 0);
		gd.addNumericField("End Slice:", imp.getStackSize(), 0);
		String[] items = { "Pixel Values", "Calibration Table",
				"Linear Fit to Pixels", "Linear Fit to Calibration", "None" };
		gd.addChoice("Density from...", items, items[calMode]);

		// gd.addCheckbox("HU Calibrated", ImageCheck.huCalibrated(imp));
		gd.addNumericField("Bone_Min:", min, 1, 6, pixUnits + " ");
		gd.addNumericField("Bone_Max:", max, 1, 6, pixUnits + " ");
		gd
				.addMessage("Only pixels >= bone min\n"
						+ "and <= bone max are used.");
		gd.addMessage("Density calibration coefficients");
		gd.addNumericField("Slope", 0, 4, 6, "g.cm^-3 / " + pixUnits + " ");
		gd.addNumericField("Y_Intercept", 1.8, 4, 6, "g.cm^-3");
		gd.addCheckbox("Align result", true);
		gd.addCheckbox("Show axes (2D)", false);
		gd.addCheckbox("Show axes (3D)", true);
		gd.addHelp("http://bonej.org/moments");
		gd.addDialogListener(this);
		gd.showDialog();
		if (gd.wasCanceled()) {
			return;
		}
		final int startSlice = (int) gd.getNextNumber();
		final int endSlice = (int) gd.getNextNumber();
		calMode = gd.getNextChoiceIndex();
		// boolean isHUCalibrated = gd.getNextBoolean();
		min = gd.getNextNumber();
		max = gd.getNextNumber();
		// if (isHUCalibrated) {
		// min = cal.getRawValue(min);
		// max = cal.getRawValue(max);
		// }
		if (calMode == TABLE || calMode == LINEAR_CAL){
			min = cal.getRawValue(min);
			max = cal.getRawValue(max);
		}
		final double m = gd.getNextNumber();
		final double c = gd.getNextNumber();
		final boolean doAlign = gd.getNextBoolean();
		final boolean doAxes = gd.getNextBoolean();
		final boolean doAxes3D = gd.getNextBoolean();

		double[] centroid = getCentroid3D(imp, startSlice, endSlice, min, max,
				m, c);
		if (centroid[0] < 0) {
			IJ.error("Empty Stack",
					"No voxels are available for calculation.\n"
							+ "Check your ROI and threshold.");
			return;
		}
		Object[] momentResults = calculateMoments(imp, startSlice, endSlice,
				centroid, min, max, m, c);

		EigenvalueDecomposition E = (EigenvalueDecomposition) momentResults[0];
		double[] moments = (double[]) momentResults[1];

		String units = imp.getCalibration().getUnits();
		ResultInserter ri = ResultInserter.getInstance();
		ri.setResultInRow(imp, "Xc (" + units + ")", centroid[0]);
		ri.setResultInRow(imp, "Yc (" + units + ")", centroid[1]);
		ri.setResultInRow(imp, "Zc (" + units + ")", centroid[2]);
		ri.setResultInRow(imp, "Vol (" + units + "³)", moments[0]);
		ri.setResultInRow(imp, "Mass (g)", moments[1]);
		ri.setResultInRow(imp, "Icxx (kg.m²)", moments[2]);
		ri.setResultInRow(imp, "Icyy (kg.m²)", moments[3]);
		ri.setResultInRow(imp, "Iczz (kg.m²)", moments[4]);
		ri.setResultInRow(imp, "Icxy (kg.m²)", moments[5]);
		ri.setResultInRow(imp, "Icxz (kg.m²)", moments[6]);
		ri.setResultInRow(imp, "Icyz (kg.m²)", moments[7]);
		ri.setResultInRow(imp, "I1 (kg.m²)", E.getD().get(2, 2));
		ri.setResultInRow(imp, "I2 (kg.m²)", E.getD().get(1, 1));
		ri.setResultInRow(imp, "I3 (kg.m²)", E.getD().get(0, 0));
		ri.updateTable();

		if (doAlign)
			alignToPrincipalAxes(imp, E.getV(), centroid, startSlice, endSlice,
					min, max, doAxes).show();

		if (doAxes3D)
			show3DAxes(imp, E.getV(), centroid, startSlice, endSlice, min, max);

		return;
	}

	// cortical bone apparent density (material density * volume fraction) from
	// Mow & Huiskes (2005) p.140
	// using 1.8 g.cm^-3: 1mm³ = 0.0018 g = 0.0000018 kg = 1.8*10^-6 kg; 1mm²
	// = 10^-6 m²
	// conversion coefficient from mm^5 to kg.m² = 1.8*10^-12
	// double cc = 1.8*Math.pow(10, -12);

	/**
	 * Convert a pixel value <i>x</i> to a voxel density <i>y</i> given
	 * calibration constants <i>m</i>, <i>c</i> for the equation <i>y</i> =
	 * <i>mx</i> + <i>c</i>
	 * 
	 * @param pixelValue
	 *            raw pixel value
	 * @param m
	 *            slope of regression line
	 * @param c
	 *            y intercept of regression line
	 * 
	 * @return voxelDensity
	 */
	private double voxelDensity(double pixelValue, double m, double c,
			double factor) {
		double voxelDensity = (m * pixelValue + c) / factor;
		if (voxelDensity < 0)
			voxelDensity = 0;
		return voxelDensity;
	}/* end voxelDensity */

	/**
	 * Get a scale factor because density is in g / cm³ but our units are mm so
	 * density is 1000* too high
	 * 
	 * @param imp
	 * @return
	 */
	private double getDensityFactor(ImagePlus imp) {
		String units = imp.getCalibration().getUnits();
		double factor = 1;
		if (units.contains("mm")) {
			factor = 1000;
		} else {
			factor = 1;
		}
		return factor;
	}

	/**
	 * Return an empty pixel array of the type appropriate for the bit depth
	 * required. Returns an Object, which can be used when adding an empty slice
	 * to a stack
	 * 
	 * @param w
	 * @param h
	 * @param bitDepth
	 * @return Object containing an array of the type needed for an image with
	 *         bitDepth
	 */
	public static Object getEmptyPixels(int w, int h, int bitDepth) {

		Object emptyPixels = new Object();
		if (bitDepth == 8) {
			byte[] bytePixels = new byte[w * h];
			emptyPixels = bytePixels;
		} else if (bitDepth == 16) {
			short[] shortPixels = new short[w * h];
			emptyPixels = shortPixels;
		} else if (bitDepth == 32) {
			float[] floatPixels = new float[w * h];
			emptyPixels = floatPixels;
		}
		return emptyPixels;
	}

	/**
	 * Calculate a density-weighted centroid in an image using z-clip planes,
	 * threshold clipping and density = m * pixel value + c density equation
	 * 
	 * @param imp
	 *            ImagePlus
	 * @param startSlice
	 *            first slice to use
	 * @param endSlice
	 *            last slice to use
	 * @param min
	 *            minimum threshold value
	 * @param max
	 *            maximum threshold value
	 * @param m
	 *            slope of density equation (set to 0 if constant density)
	 * @param c
	 *            constant in density equation
	 * @return double[] containing (x,y,z) centroid in scaled units
	 */
	public double[] getCentroid3D(ImagePlus imp, int startSlice, int endSlice,
			final double min, final double max, final double m, final double c) {
		final ImageStack stack = imp.getImageStack();
		final Rectangle r = imp.getProcessor().getRoi();
		final int rW = r.x + r.width;
		final int rH = r.y + r.height;
		final int rX = r.x;
		final int rY = r.y;
		Calibration cal = imp.getCalibration();
		final double vW = cal.pixelWidth;
		final double vH = cal.pixelHeight;
		final double vD = cal.pixelDepth;
		final double voxVol = vW * vH * vD;
		final double factor = getDensityFactor(imp);
		double sumx = 0;
		double sumy = 0;
		double sumz = 0;
		double sumMass = 0;
		for (int z = startSlice; z <= endSlice; z++) {
			IJ.showStatus("Calculating centroid...");
			IJ.showProgress(z - startSlice, endSlice - startSlice);
			ImageProcessor ip = stack.getProcessor(z);
			for (int y = rY; y < rH; y++) {
				for (int x = rX; x < rW; x++) {
					final double testPixel = ip.get(x, y);
					if (testPixel >= min && testPixel <= max) {
						final double voxelMass = voxelDensity(testPixel, m, c,
								factor)
								* voxVol;
						sumMass += voxelMass;
						sumx += x * voxelMass;
						sumy += y * voxelMass;
						sumz += z * voxelMass;
					}
				}
			}
		}
		// centroid in pixels
		double centX = sumx / sumMass;
		double centY = sumy / sumMass;
		double centZ = sumz / sumMass;
		if (sumMass == 0) {
			centX = -1;
			centY = -1;
			centZ = -1;
		}
		// centroid in real units
		double[] centroid = { centX * vW, centY * vH, centZ * vD };
		return centroid;
	}/* end findCentroid3D */

	public Object[] calculateMoments(ImagePlus imp, int startSlice,
			int endSlice, double[] centroid, final double min,
			final double max, double m, double c) {
		// START OF 3D MOMENT CALCULATIONS
		Calibration cal = imp.getCalibration();
		final double vW = cal.pixelWidth;
		final double vH = cal.pixelHeight;
		final double vD = cal.pixelDepth;
		final ImageStack stack = imp.getImageStack();
		final Rectangle r = imp.getProcessor().getRoi();
		final int rW = r.x + r.width;
		final int rH = r.y + r.height;
		final int rX = r.x;
		final int rY = r.y;
		final double cX = centroid[0];
		final double cY = centroid[1];
		final double cZ = centroid[2];
		final double factor = getDensityFactor(imp);
		final double voxVol = vW * vH * vD;
		final double voxVhVd = (vH * vH + vD * vD) / 12;
		final double voxVwVd = (vW * vW + vD * vD) / 12;
		final double voxVhVw = (vH * vH + vW * vW) / 12;
		double sumVoxVol = 0;
		double sumVoxMass = 0;
		double Icxx = 0;
		double Icyy = 0;
		double Iczz = 0;
		double Icxy = 0;
		double Icxz = 0;
		double Icyz = 0;
		for (int z = startSlice; z <= endSlice; z++) {
			IJ.showStatus("Calculating inertia tensor...");
			IJ.showProgress(z - startSlice, endSlice - startSlice);
			ImageProcessor ip = stack.getProcessor(z);
			for (int y = rY; y < rH; y++) {
				for (int x = rX; x < rW; x++) {
					final double testPixel = (double) ip.get(x, y);
					if (testPixel < min || testPixel > max) {
						continue;
					} else {
						sumVoxVol += voxVol;
						final double voxMass = voxelDensity(testPixel, m, c,
								factor)
								* voxVol;
						sumVoxMass += voxMass;
						final double xvWcX = x * vW - cX;
						final double yvHcY = y * vH - cY;
						final double zvDcZ = z * vD - cZ;
						Icxx += (yvHcY * yvHcY + zvDcZ * zvDcZ + voxVhVd)
								* voxMass;
						Icyy += (xvWcX * xvWcX + zvDcZ * zvDcZ + voxVwVd)
								* voxMass;
						Iczz += (yvHcY * yvHcY + xvWcX * xvWcX + voxVhVw)
								* voxMass;
						Icxy += xvWcX * yvHcY * voxMass;
						Icxz += xvWcX * zvDcZ * voxMass;
						Icyz += yvHcY * zvDcZ * voxMass;
					}
				}
			}
		}
		// create the inertia tensor matrix
		double[][] inertiaTensor = new double[3][3];
		inertiaTensor[0][0] = Icxx;
		inertiaTensor[1][1] = Icyy;
		inertiaTensor[2][2] = Iczz;
		inertiaTensor[0][1] = -Icxy;
		inertiaTensor[0][2] = -Icxz;
		inertiaTensor[1][0] = -Icxy;
		inertiaTensor[1][2] = -Icyz;
		inertiaTensor[2][0] = -Icxz;
		inertiaTensor[2][1] = -Icyz;
		Matrix inertiaTensorMatrix = new Matrix(inertiaTensor);

		// do the Eigenvalue decomposition
		EigenvalueDecomposition E = new EigenvalueDecomposition(
				inertiaTensorMatrix);
		E.getD().printToIJLog("Eigenvalues");
		E.getV().printToIJLog("Eigenvectors");

		double[] moments = { sumVoxVol, sumVoxMass, Icxx, Icyy, Iczz, Icxy,
				Icxz, Icyz };
		Object[] result = { E, moments };

		return result;
	}

	/**
	 * Draw a copy of the original image aligned to its principal axes
	 * 
	 * @param imp
	 *            Input image
	 * @param E
	 *            Rotation matrix
	 * @param centroid
	 *            3-element array containing centroid coordinates, {x,y,z}
	 * @param startSlice
	 *            first slice to copy
	 * @param endSlice
	 *            final slice to copy
	 * @param doAxes
	 *            if true, draw axes on the aligned copy
	 * @return ImagePlus copy of the input image
	 */
	public ImagePlus alignToPrincipalAxes(ImagePlus imp, Matrix E,
			double[] centroid, int startSlice, int endSlice, double min,
			double max, boolean doAxes) {
		final ImageStack sourceStack = imp.getImageStack();
		Calibration cal = imp.getCalibration();
		final double vW = cal.pixelWidth;
		final double vH = cal.pixelHeight;
		final double vD = cal.pixelDepth;
		final double vS = Math.min(vW, Math.min(vH, vD));
		final int d = sourceStack.getSize();
		int[] sides = getRotatedSize(E, imp, centroid, startSlice, endSlice,
				min, max);

		// Rotation matrix to rotate data 90 deg around x axis
		double[][] rotX = new double[3][3];
		rotX[0][0] = 1;
		rotX[1][2] = -1;
		rotX[2][1] = 1;
		Matrix RotX = new Matrix(rotX);

		// Rotation matrix to rotate data 90 deg around y axis
		double[][] rotY = new double[3][3];
		rotY[0][2] = 1;
		rotY[1][1] = 1;
		rotY[2][0] = -1;
		Matrix RotY = new Matrix(rotY);

		// Rotation matrix to rotate data 90 deg around z axis
		double[][] rotZ = new double[3][3];
		rotZ[0][1] = -1;
		rotZ[1][0] = 1;
		rotZ[2][2] = 1;
		Matrix RotZ = new Matrix(rotZ);

		final int wi = sides[0];
		final int hi = sides[1];
		final int di = sides[2];

		E.printToIJLog("Original Rotation Matrix (Source -> Target)");
		Matrix rotation = new Matrix(new double[3][3]);

		// put long axis in z, middle axis in y and short axis in x
		if (wi <= hi && hi <= di) {
			IJ.log("Case 0");
			rotation = E;
		} else if (wi <= di && di <= hi) {
			IJ.log("Case 1");
			rotation = E.times(RotX);
		} else if (hi <= wi && wi <= di) {
			IJ.log("Case 2");
			rotation = E.times(RotZ);
		} else if (di <= hi && hi <= wi) {
			IJ.log("Case 3");
			rotation = E.times(RotY);
		} else if (hi <= di && di <= wi) {
			IJ.log("Case 4");
			rotation = E.times(RotY).times(RotZ);
		} else if (di <= wi && wi <= hi) {
			IJ.log("Case 5");
			rotation = E.times(RotX).times(RotZ);
		} else {
			IJ.log("Case 6");
			rotation = E;
		}
		// Check for z-flipping
		if (rotation.isZFlipped()) {
			rotation = rotation.times(RotX).times(RotX);
			IJ.log("Corrected Z-flipping");
		}

		// Check for reflection
		if (!rotation.isRightHanded()) {
			double[][] reflectY = new double[3][3];
			reflectY[0][0] = -1;
			reflectY[1][1] = 1;
			reflectY[2][2] = 1;
			Matrix RefY = new Matrix(reflectY);
			rotation = rotation.times(RefY);
			IJ.log("Reflected the rotation matrix");
		}
		rotation.printToIJLog("Rotation Matrix (Source -> Target)");

		Matrix eVecInv = rotation.inverse();
		eVecInv.printToIJLog("Inverse Rotation Matrix (Target -> Source)");
		final double[][] eigenVecInv = eVecInv.getArrayCopy();

		// create the target stack
		Arrays.sort(sides);
		final int wT = sides[0];
		final int hT = sides[1];
		final int dT = sides[2];
		final double xTc = wT * vS / 2;
		final double yTc = hT * vS / 2;
		final double zTc = dT * vS / 2;

		// for each voxel in the target stack,
		// find the corresponding source voxel

		// Cache the sourceStack's processors
		ImageProcessor[] sliceProcessors = new ImageProcessor[d + 1];
		for (int z = 1; z <= d; z++) {
			sliceProcessors[z] = sourceStack.getProcessor(z);
		}
		// Initialise an empty stack and tartgetStack's processor array
		ImageStack targetStack = new ImageStack(wT, hT, dT);
		ImageProcessor[] targetProcessors = new ImageProcessor[dT + 1];
		for (int z = 1; z <= dT; z++) {
			targetStack.setPixels(getEmptyPixels(wT, hT, imp.getBitDepth()), z);
			targetProcessors[z] = targetStack.getProcessor(z);
		}

		// Multithread start
		int nThreads = Runtime.getRuntime().availableProcessors();
		AlignThread[] alignThread = new AlignThread[nThreads];
		for (int thread = 0; thread < nThreads; thread++) {
			alignThread[thread] = new AlignThread(thread, nThreads, imp,
					sliceProcessors, targetProcessors, eigenVecInv, centroid,
					wT, hT, dT, startSlice, endSlice);
			alignThread[thread].start();
		}
		try {
			for (int thread = 0; thread < nThreads; thread++) {
				alignThread[thread].join();
			}
		} catch (InterruptedException ie) {
			IJ.error("A thread was interrupted.");
		}
		// end multithreading
		if (doAxes) {
			// draw axes on stack
			final int xCent = (int) Math.floor(xTc / vS);
			final int yCent = (int) Math.floor(yTc / vS);
			final int zCent = (int) Math.floor(zTc / vS);
			final int axisColour = Integer.MAX_VALUE;
			for (int z = 1; z <= dT; z++) {
				// z axis
				try {
					ImageProcessor axisIP = targetStack.getProcessor(z);
					axisIP.set(xCent, yCent, axisColour);
				} catch (NullPointerException npe) {
					IJ.handleException(npe);
					break;
				}
			}
			ImageProcessor axisIP = targetStack.getProcessor(zCent);
			axisIP.setColor(Integer.MAX_VALUE);
			// x axis
			axisIP.drawLine(0, yCent, wT, yCent);

			// y axis
			axisIP.drawLine(xCent, 0, xCent, hT);
		}
		ImagePlus impTarget = new ImagePlus("Aligned_" + imp.getTitle(),
				targetStack);
		impTarget.setCalibration(imp.getCalibration());
		Calibration targetCal = impTarget.getCalibration();
		targetCal.pixelDepth = vS;
		targetCal.pixelHeight = vS;
		targetCal.pixelWidth = vS;
		impTarget.setDisplayRange(imp.getDisplayRangeMin(), imp
				.getDisplayRangeMax());
		return impTarget;
	}

	/**
	 * Multithreading class to look up aligned voxel values, processing each
	 * slice in its own thread
	 * 
	 * @author Michael Doube
	 * 
	 */
	class AlignThread extends Thread {
		final int thread, nThreads, wT, hT, dT, startSlice, endSlice;
		final ImagePlus impT;
		final ImageStack stackT;
		final ImageProcessor[] sliceProcessors, targetProcessors;
		final double[][] eigenVecInv;
		final double[] centroid;

		public AlignThread(int thread, int nThreads, ImagePlus imp,
				ImageProcessor[] sliceProcessors,
				ImageProcessor[] targetProcessors, double[][] eigenVecInv,
				double[] centroid, int wT, int hT, int dT, int startSlice,
				int endSlice) {
			this.impT = imp;
			this.stackT = this.impT.getStack();
			this.thread = thread;
			this.nThreads = nThreads;
			this.sliceProcessors = sliceProcessors;
			this.targetProcessors = targetProcessors;
			this.eigenVecInv = eigenVecInv;
			this.centroid = centroid;
			this.wT = wT;
			this.hT = hT;
			this.dT = dT;
			this.startSlice = startSlice;
			this.endSlice = endSlice;
		}

		public void run() {
			final Rectangle r = this.impT.getProcessor().getRoi();
			final int rW = r.x + r.width;
			final int rH = r.y + r.height;
			final int rX = r.x;
			final int rY = r.y;
			Calibration cal = this.impT.getCalibration();
			final double vW = cal.pixelWidth;
			final double vH = cal.pixelHeight;
			final double vD = cal.pixelDepth;
			final double vS = Math.min(vW, Math.min(vH, vD));
			final double xC = centroid[0];
			final double yC = centroid[1];
			final double zC = centroid[2];
			final double xTc = wT * vS / 2;
			final double yTc = hT * vS / 2;
			final double zTc = dT * vS / 2;
			final double dXc = xC - xTc;
			final double dYc = yC - yTc;
			final double dZc = zC - zTc;
			final double eVI00 = eigenVecInv[0][0];
			final double eVI10 = eigenVecInv[1][0];
			final double eVI20 = eigenVecInv[2][0];
			final double eVI01 = eigenVecInv[0][1];
			final double eVI11 = eigenVecInv[1][1];
			final double eVI21 = eigenVecInv[2][1];
			final double eVI02 = eigenVecInv[0][2];
			final double eVI12 = eigenVecInv[1][2];
			final double eVI22 = eigenVecInv[2][2];
			for (int z = this.thread + 1; z <= this.dT; z += this.nThreads) {
				IJ.showStatus("Aligning image stack...");
				IJ.showProgress(z, this.dT);
				// this.targetStack.setPixels(getEmptyPixels(this.wT, this.hT,
				// this.impT.getBitDepth()), z);
				ImageProcessor targetIP = targetProcessors[z];
				final double zD = z * vS - zTc;
				final double zDeVI00 = zD * eVI20;
				final double zDeVI01 = zD * eVI21;
				final double zDeVI02 = zD * eVI22;
				for (int y = 0; y < this.hT; y++) {
					final double yD = y * vS - yTc;
					final double yDeVI10 = yD * eVI10;
					final double yDeVI11 = yD * eVI11;
					final double yDeVI12 = yD * eVI12;
					for (int x = 0; x < this.wT; x++) {
						final double xD = x * vS - xTc;
						final double xAlign = xD * eVI00 + yDeVI10 + zDeVI00
								+ xTc;
						final double yAlign = xD * eVI01 + yDeVI11 + zDeVI01
								+ yTc;
						final double zAlign = xD * eVI02 + yDeVI12 + zDeVI02
								+ zTc;
						// possibility to do some voxel interpolation instead
						// of just rounding in next 3 lines
						final int xA = (int) Math.floor((xAlign + dXc) / vW);
						final int yA = (int) Math.floor((yAlign + dYc) / vH);
						final int zA = (int) Math.floor((zAlign + dZc) / vD);

						if (xA < rX || xA >= rW || yA < rY || yA >= rH
								|| zA < this.startSlice || zA > this.endSlice) {
							continue;
						} else {
							targetIP.set(x, y, this.sliceProcessors[zA].get(xA,
									yA));
						}
					}
				}
			}
		}
	}

	/**
	 * Find side lengths in pixels of the smallest stack to fit the aligned
	 * image
	 * 
	 * @param E
	 *            Rotation matrix
	 * @param imp
	 *            Source image
	 * @param centroid
	 *            3D centroid in 3-element array {x,y,z}
	 * @param startSlice
	 *            first slice of source image
	 * @param endSlice
	 *            last slice of source image
	 * @param min
	 *            minimum threshold
	 * @param max
	 *            maximum threshold
	 * @return Width, height and depth of a stack that will 'just fit' the
	 *         aligned image
	 */
	private int[] getRotatedSize(Matrix E, ImagePlus imp, double[] centroid,
			int startSlice, int endSlice, double min, double max) {
		final ImageStack stack = imp.getImageStack();
		final Calibration cal = imp.getCalibration();
		final double xC = centroid[0];
		final double yC = centroid[1];
		final double zC = centroid[2];
		final Rectangle r = imp.getProcessor().getRoi();
		final int rW = r.x + r.width;
		final int rH = r.y + r.height;
		final int rX = r.x;
		final int rY = r.y;

		final double vW = cal.pixelWidth;
		final double vH = cal.pixelHeight;
		final double vD = cal.pixelDepth;

		final double[][] v = E.getArrayCopy();
		final double v00 = v[0][0];
		final double v10 = v[1][0];
		final double v20 = v[2][0];
		final double v01 = v[0][1];
		final double v11 = v[1][1];
		final double v21 = v[2][1];
		final double v02 = v[0][2];
		final double v12 = v[1][2];
		final double v22 = v[2][2];

		double xTmax = 0;
		double yTmax = 0;
		double zTmax = 0;

		for (int z = startSlice; z <= endSlice; z++) {
			IJ.showStatus("Getting aligned stack dimensions...");
			ImageProcessor ip = stack.getProcessor(z);
			final double zCz = z * vD - zC;
			final double zCzv20 = zCz * v20;
			final double zCzv21 = zCz * v21;
			final double zCzv22 = zCz * v22;
			for (int y = rY; y < rH; y++) {
				final double yCy = y * vH - yC;
				final double yCyv10 = yCy * v10;
				final double yCyv11 = yCy * v11;
				final double yCyv12 = yCy * v12;
				for (int x = rX; x < rW; x++) {
					final double pixel = ip.get(x, y);
					if (pixel < min || pixel > max)
						continue;
					else {
						// distance from centroid in
						// original coordinate system
						// xCx, yCx, zCx
						final double xCx = x * vW - xC;

						// now transform each coordinate
						// transformed coordinate is dot product of original
						// coordinates
						// and eigenvectors
						final double xT = xCx * v00 + yCyv10 + zCzv20;
						final double yT = xCx * v01 + yCyv11 + zCzv21;
						final double zT = xCx * v02 + yCyv12 + zCzv22;

						// keep the biggest value to find the greatest distance
						// in x, y and z
						xTmax = Math.max(xTmax, Math.abs(xT));
						yTmax = Math.max(yTmax, Math.abs(yT));
						zTmax = Math.max(zTmax, Math.abs(zT));
					}
				}
			}
		}

		// use the smallest input voxel dimension as the voxel size
		double vS = Math.min(vW, Math.min(vH, vD));

		int tW = (int) Math.floor(2 * xTmax / vS) + 5;
		int tH = (int) Math.floor(2 * yTmax / vS) + 5;
		int tD = (int) Math.floor(2 * zTmax / vS) + 5;

		int[] size = { tW, tH, tD };
		return size;
	}

	/**
	 * Create a copy of the original image aligned to the tensor defined by a
	 * 3x3 Eigenvector matrix
	 * 
	 * @param imp
	 *            input ImagePlus stack
	 * @param E
	 *            Eigenvector rotation matrix
	 * @param doAxes
	 *            if true, draws axes on the aligned image
	 * @param startSlice
	 *            first slice to use
	 * @param endSlice
	 *            last slice to use
	 * @param min
	 *            minimum threshold
	 * @param max
	 *            maximum threshold
	 * @param m
	 *            slope of pixel to density equation, d = m * p + c
	 * @param c
	 *            intercept of density equation, d = m * p + c
	 * @return aligned ImagePlus
	 */
	public ImagePlus alignImage(ImagePlus imp, Matrix E, boolean doAxes,
			int startSlice, int endSlice, double min, double max, double m,
			double c) {
		final double[] centroid = getCentroid3D(imp, startSlice, endSlice, min,
				max, m, c);
		ImagePlus alignedImp = alignToPrincipalAxes(imp, E, centroid,
				startSlice, endSlice, min, max, doAxes);
		return alignedImp;
	}

	/**
	 * Display principal axes on a 3D rendered version of the image
	 * 
	 * @param imp
	 *            Original image
	 * @param E
	 *            eigenvectors of the principal axes
	 * @param centroid
	 *            in real units
	 * @param startSlice
	 *            first slice
	 * @param endSlice
	 *            last slice
	 * @param min
	 *            lower threshold
	 * @param max
	 *            upper threshold
	 */
	private void show3DAxes(ImagePlus imp, Matrix E, double[] centroid,
			int startSlice, int endSlice, double min, double max) {
		Calibration cal = imp.getCalibration();
		// copy the data from inside the ROI and convert it to 8-bit
		Duplicator d = new Duplicator();
		ImagePlus roiImp = d.run(imp, startSlice, endSlice);
		if (imp.getRoi() != null) {
			Rectangle roi = imp.getRoi().getBounds();
			centroid[0] -= roi.getX() * cal.pixelWidth;
			centroid[1] -= roi.getY() * cal.pixelHeight;
		}
		centroid[2] -= startSlice * cal.pixelDepth;
		final double cX = centroid[0];
		final double cY = centroid[1];
		final double cZ = centroid[2];

		Point3f cent = new Point3f();
		cent.x = (float) cX;
		cent.y = (float) cY;
		cent.z = (float) cZ;

		// initialise and show the 3D universe
		Image3DUniverse univ = new Image3DUniverse();
		univ.show();

		// show the centroid
		ArrayList<Point3f> point = new ArrayList<Point3f>();
		point.add(cent);
		CustomPointMesh mesh = new CustomPointMesh(point);
		mesh.setPointSize(5.0f);
		float red = 0.0f;
		float green = 0.5f;
		float blue = 1.0f;
		Color3f cColour = new Color3f(red, green, blue);
		mesh.setColor(cColour);
		try {
			univ.addCustomMesh(mesh, "Centroid").setLocked(true);
		} catch (NullPointerException npe) {
			IJ.log("3D Viewer was closed before rendering completed.");
			return;
		}

		// show the axes
		int[] sideLengths = getRotatedSize(E, roiImp, centroid, 1, endSlice
				- startSlice + 1, min, max);
		final double vS = Math.min(cal.pixelWidth, Math.min(cal.pixelHeight,
				cal.pixelDepth));
		final double l1 = sideLengths[0] * vS;
		final double l2 = sideLengths[1] * vS;
		final double l3 = sideLengths[2] * vS;
		List<Point3f> axes = new ArrayList<Point3f>();
		Point3f start1 = new Point3f();
		start1.x = (float) (cX - E.get(0, 0) * l1);
		start1.y = (float) (cY - E.get(1, 0) * l1);
		start1.z = (float) (cZ - E.get(2, 0) * l1);
		axes.add(start1);

		Point3f end1 = new Point3f();
		end1.x = (float) (cX + E.get(0, 0) * l1);
		end1.y = (float) (cY + E.get(1, 0) * l1);
		end1.z = (float) (cZ + E.get(2, 0) * l1);
		axes.add(end1);

		Point3f start2 = new Point3f();
		start2.x = (float) (cX - E.get(0, 1) * l2);
		start2.y = (float) (cY - E.get(1, 1) * l2);
		start2.z = (float) (cZ - E.get(2, 1) * l2);
		axes.add(start2);

		Point3f end2 = new Point3f();
		end2.x = (float) (cX + E.get(0, 1) * l2);
		end2.y = (float) (cY + E.get(1, 1) * l2);
		end2.z = (float) (cZ + E.get(2, 1) * l2);
		axes.add(end2);

		Point3f start3 = new Point3f();
		start3.x = (float) (cX - E.get(0, 2) * l3);
		start3.y = (float) (cY - E.get(1, 2) * l3);
		start3.z = (float) (cZ - E.get(2, 2) * l3);
		axes.add(start3);

		Point3f end3 = new Point3f();
		end3.x = (float) (cX + E.get(0, 2) * l3);
		end3.y = (float) (cY + E.get(1, 2) * l3);
		end3.z = (float) (cZ + E.get(2, 2) * l3);
		axes.add(end3);

		Color3f aColour = new Color3f(red, green, blue);
		try {
			univ.addLineMesh(axes, aColour, "Principal axes", false).setLocked(
					true);
		} catch (NullPointerException npe) {
			IJ.log("3D Viewer was closed before rendering completed.");
			return;
		}

		// show the stack
		try {
			new StackConverter(roiImp).convertToGray8();
			Content c = univ.addVoltex(roiImp);
			c.setLocked(true);
		} catch (NullPointerException npe) {
			IJ.log("3D Viewer was closed before rendering completed.");
			return;
		}
		return;
	}

	public boolean dialogItemChanged(GenericDialog gd, AWTEvent e) {
		if (!DialogModifier.allNumbersValid(gd.getNumericFields()))
			return false;
		// Vector<?> checkboxes = gd.getCheckboxes();
		Vector<?> nFields = gd.getNumericFields();
		Vector<?> choices = gd.getChoices();
		// Checkbox box0 = (Checkbox) checkboxes.get(0);
		// boolean isHUCalibrated = box0.getState();
		TextField minT = (TextField) nFields.get(2);
		TextField maxT = (TextField) nFields.get(3);
		TextField slope = (TextField) nFields.get(4);
		TextField intercept = (TextField) nFields.get(5);
		double min = Double.parseDouble(minT.getText());
		double max = Double.parseDouble(maxT.getText());
		if (calMode == TABLE || calMode == LINEAR_CAL){
			min = cal.getRawValue(min);
			max = cal.getRawValue(max);
		}
		Choice calChoice = (Choice) choices.get(0);
		if (e.getSource().equals(calChoice)) {
			fieldUpdated = false;
			calMode = calChoice.getSelectedIndex();
			switch (calMode) {
			case PIXEL:
				updateUnits(gd, "grey");
				slope.setEnabled(false);
				intercept.setEnabled(false);
				if (!fieldUpdated) {
					minT.setText("" + min);
					maxT.setText("" + max);
					fieldUpdated = true;
				}
				break;
			case TABLE:
				updateUnits(gd, cal.getValueUnit());
				slope.setEnabled(false);
				intercept.setEnabled(false);
				if (!fieldUpdated) {
					minT.setText("" + cal.getCValue(min));
					maxT.setText("" + cal.getCValue(max));
					fieldUpdated = true;
				}
				break;
			case LINEAR_CAL:
				updateUnits(gd, cal.getValueUnit());
				slope.setEnabled(true);
				intercept.setEnabled(true);
				if (!fieldUpdated) {
					minT.setText("" + cal.getCValue(min));
					maxT.setText("" + cal.getCValue(max));
					fieldUpdated = true;
				}
				break;
			case LINEAR_PIXEL:
				updateUnits(gd, "grey");
				slope.setEnabled(true);
				intercept.setEnabled(true);
				if (!fieldUpdated) {
					minT.setText("" + min);
					maxT.setText("" + max);
					fieldUpdated = true;
				}
				break;
			case NONE:
				updateUnits(gd, "grey");
				slope.setEnabled(false);
				intercept.setEnabled(false);
				if (!fieldUpdated) {
					minT.setText("" + min);
					maxT.setText("" + max);
					fieldUpdated = true;
				}
				break;
			}
		}
		// if (isHUCalibrated && !fieldUpdated) {
		// minT.setText("" + cal.getCValue(min));
		// maxT.setText("" + cal.getCValue(max));
		// fieldUpdated = true;
		// }
		// // if (!isHUCalibrated && fieldUpdated) {
		// // minT.setText("" + cal.getRawValue(min));
		// // maxT.setText("" + cal.getRawValue(max));
		// // fieldUpdated = false;
		// // }
		// if (isHUCalibrated)
		// updateUnits(gd, "HU");
		DialogModifier.registerMacroValues(gd, gd.getComponents());
		return true;
	}

	private void updateUnits(GenericDialog gd, String newUnits) {
		DialogModifier.replaceUnitString(gd, pixUnits, newUnits);
		pixUnits = newUnits;
	}

}
