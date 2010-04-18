package org.doube.bonej;

/** Moments 3D
 *tool to calculate centroid and principal axes 
 *of a thresholded stack; originally designed for 16-bit CT scans 
 *of a bone in air so default thresholds are 0 and 4000 HU, but most greyscale images shoudl be handled
 *Outputs stack data aligned to principal axes
 *
 *Copyright 2008 2009 2010 Michael Doube 
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
 *@author Michael Doube
 */

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ImageProcessor;
import ij.plugin.PlugIn;
import ij.measure.Calibration;
import ij.gui.*;

import java.awt.AWTEvent;
import java.awt.Checkbox;
import java.awt.Rectangle;
import java.awt.TextField;
import java.util.Arrays;
import java.util.Vector;

import org.doube.jama.Matrix;
import org.doube.jama.EigenvalueDecomposition;
import org.doube.util.DialogModifier;
import org.doube.util.ImageCheck;
import org.doube.util.ResultInserter;
import org.doube.util.ThresholdGuesser;

public class Moments implements PlugIn, DialogListener {

	private boolean fieldUpdated = false;
	private Calibration cal;

	public void run(String arg) {
		if (!ImageCheck.checkIJVersion())
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
		String pixUnits;
		if (ImageCheck.huCalibrated(imp)) {
			pixUnits = "HU";
			fieldUpdated = true;
		} else
			pixUnits = "grey";
		GenericDialog gd = new GenericDialog("Setup");
		gd.addNumericField("Start Slice:", 1, 0);
		gd.addNumericField("End Slice:", imp.getStackSize(), 0);

		gd.addCheckbox("HU Calibrated", ImageCheck.huCalibrated(imp));
		gd.addNumericField("Bone Min:", min, 1, 6, pixUnits + " ");
		gd.addNumericField("Bone Max:", max, 1, 6, pixUnits + " ");
		gd
				.addMessage("Only pixels >= bone min\n"
						+ "and <= bone max are used.");
		gd.addMessage("Density calibration coefficients");
		gd.addNumericField("Slope", 0, 4, 6, "g.cm^-3 / " + pixUnits + " ");
		gd.addNumericField("Y_Intercept", 1.8, 4, 6, "g.cm^-3");
		gd.addMessage("Only use pixels between clip values:");
		gd.addCheckbox("Align result", true);
		gd.addCheckbox("Show axes", true);
		gd.addDialogListener(this);
		gd.showDialog();
		if (gd.wasCanceled()) {
			return;
		}
		final int startSlice = (int) gd.getNextNumber();
		final int endSlice = (int) gd.getNextNumber();
		boolean isHUCalibrated = gd.getNextBoolean();
		min = gd.getNextNumber();
		max = gd.getNextNumber();
		if (isHUCalibrated) {
			min = cal.getRawValue(min);
			max = cal.getRawValue(max);
		}
		final double m = gd.getNextNumber();
		final double c = gd.getNextNumber();
		final boolean doAlign = gd.getNextBoolean();
		final boolean doAxes = gd.getNextBoolean();

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
			alignToPrincipalAxes(imp, E, centroid, startSlice, endSlice, min,
					max, doAxes).show();
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
		byte[] bytePixels = new byte[w * h];
		short[] shortPixels = new short[w * h];
		float[] floatPixels = new float[w * h];

		Object emptyPixels = new Object();
		if (bitDepth == 8) {
			emptyPixels = bytePixels;
		} else if (bitDepth == 16) {
			emptyPixels = shortPixels;
		} else if (bitDepth == 32) {
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
		IJ.log("Eigenvalues:");
		E.getD().printToIJLog();
		IJ.log("Eigenvectors:");
		E.getV().printToIJLog();

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
	 *            EigenvalueDecomposition of moments of inertia
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
	public ImagePlus alignToPrincipalAxes(ImagePlus imp,
			EigenvalueDecomposition E, double[] centroid, int startSlice,
			int endSlice, double min, double max, boolean doAxes) {
		final ImageStack sourceStack = imp.getImageStack();
		final Rectangle r = imp.getProcessor().getRoi();
		final int rW = r.x + r.width;
		final int rH = r.y + r.height;
		final int rX = r.x;
		final int rY = r.y;
		Calibration cal = imp.getCalibration();
		final double vW = cal.pixelWidth;
		final double vH = cal.pixelHeight;
		final double vD = cal.pixelDepth;
		final int d = sourceStack.getSize();
		ImageStack targetStack = getRotatedStack(E, imp, centroid, startSlice,
				endSlice, min, max);
		final int wT = targetStack.getWidth();
		final int hT = targetStack.getHeight();
		final int dT = targetStack.getSize();
		final double xC = centroid[0];
		final double yC = centroid[1];
		final double zC = centroid[2];
		final double xTc = wT * vW / 2;
		final double yTc = hT * vH / 2;
		final double zTc = dT * vD / 2;
		final double dXc = xC - xTc;
		final double dYc = yC - yTc;
		final double dZc = zC - zTc;

		Matrix eVec = E.getV();
		double[][] eigenVectors = eVec.getArrayCopy();
		// check orientation of eigenvectors and correct them
		if (eigenVectors[2][0] < 0) {
			double[][] eVecA = eVec.getArray();
			// rotate 180 deg
			for (int row = 0; row < 3; row++) {
				for (int col = 0; col < 2; col++) {
					eVecA[row][col] *= -1;
				}
			}
			eigenVectors = eVec.getArrayCopy();
		}
		IJ.log("\nEigenvector matrix");
		eVec.printToIJLog();

		Matrix eVecInv = eVec.inverse();
		IJ.log("\nInverse Eigenvector matrix");
		eVecInv.printToIJLog();
		final double[][] eigenVecInv = eVecInv.getArrayCopy();
		final double eVI00 = eigenVecInv[0][0];
		final double eVI10 = eigenVecInv[1][0];
		final double eVI20 = eigenVecInv[2][0];
		final double eVI01 = eigenVecInv[0][1];
		final double eVI11 = eigenVecInv[1][1];
		final double eVI21 = eigenVecInv[2][1];
		final double eVI02 = eigenVecInv[0][2];
		final double eVI12 = eigenVecInv[1][2];
		final double eVI22 = eigenVecInv[2][2];

		// for each voxel in the target stack,
		// find the corresponding source voxel

		// Cache the sourceStack's processors
		ImageProcessor[] sliceProcessors = new ImageProcessor[d + 1];
		for (int z = 1; z <= d; z++) {
			sliceProcessors[z] = sourceStack.getProcessor(z);
		}
		for (int z = 1; z <= dT; z++) {
			IJ.showStatus("Aligning image stack...");
			IJ.showProgress(z, dT);
			targetStack.setPixels(getEmptyPixels(wT, hT, imp.getBitDepth()), z);
			ImageProcessor targetIP = targetStack.getProcessor(z);
			final double zD = z * vD - zTc;
			final double zDeVI00 = zD * eVI00;
			final double zDeVI01 = zD * eVI01;
			final double zDeVI02 = zD * eVI02;
			for (int y = 0; y < hT; y++) {
				final double yD = y * vH - yTc;
				final double yDeVI10 = yD * eVI10;
				final double yDeVI11 = yD * eVI11;
				final double yDeVI12 = yD * eVI12;
				for (int x = 0; x < wT; x++) {
					final double xD = x * vW - xTc;
					final double xAlign = xD * eVI20 + yDeVI10 + zDeVI00 + xTc;
					final double yAlign = xD * eVI21 + yDeVI11 + zDeVI01 + yTc;
					final double zAlign = xD * eVI22 + yDeVI12 + zDeVI02 + zTc;
					// possibility to do some voxel interpolation instead
					// of just rounding in next 3 lines
					final int xA = (int) Math.floor((xAlign + dXc) / vW);
					final int yA = (int) Math.floor((yAlign + dYc) / vH);
					final int zA = (int) Math.floor((zAlign + dZc) / vD);

					if (xA < rX || xA >= rW || yA < rY || yA >= rH
							|| zA < startSlice || zA > endSlice) {
						continue;
					} else {
						targetIP.set(x, y, sliceProcessors[zA].get(xA, yA));
					}
				}
			}
		}
		if (doAxes) {
			// draw axes on stack
			final int xCent = (int) Math.floor(xTc / vW);
			final int yCent = (int) Math.floor(yTc / vH);
			final int zCent = (int) Math.floor(zTc / vD);
			final int axisColour = Integer.MAX_VALUE;
			for (int z = 1; z <= dT; z++) {
				ImageProcessor axisIP = targetStack.getProcessor(z);
				// z axis
				axisIP.set(xCent, yCent, axisColour);
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
		impTarget.setDisplayRange(imp.getDisplayRangeMin(), imp
				.getDisplayRangeMax());
		return impTarget;
	}

	/**
	 * Find the smallest stack to fit the aligned image
	 * 
	 * @param E
	 *            EigenvalueDecomposition of source image moments of inertia
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
	 * @return ImageStack that will 'just fit' the aligned image
	 */
	private ImageStack getRotatedStack(EigenvalueDecomposition E,
			ImagePlus imp, double[] centroid, int startSlice, int endSlice,
			double min, double max) {
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

		final Matrix V = E.getV();
		final double[][] v = V.getArrayCopy();
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
			final double zCzv00 = zCz * v00;
			final double zCzv01 = zCz * v01;
			final double zCzv02 = zCz * v02;
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
						final double xT = xCx * v20 + yCyv10 + zCzv00;
						final double yT = xCx * v21 + yCyv11 + zCzv01;
						final double zT = xCx * v22 + yCyv12 + zCzv02;

						// keep the biggest value to find the greatest distance
						// in x, y and z
						xTmax = Math.max(xTmax, Math.abs(xT));
						yTmax = Math.max(yTmax, Math.abs(yT));
						zTmax = Math.max(zTmax, Math.abs(zT));
					}
				}
			}
		}
		// TODO this still doesn't quite work properly
		// sometimes axes are in the wrong order, or mapping is not quite right
		// poss longest axis not always lowest moment - max distance vs moment!
		double[] dimensions = { xTmax, yTmax, zTmax };
		Arrays.sort(dimensions);

		xTmax = dimensions[0];
		yTmax = dimensions[1];
		zTmax = dimensions[2];

		int tW = (int) Math.floor(2 * xTmax / vW) + 3;
		int tH = (int) Math.floor(2 * yTmax / vH) + 3;
		int tD = (int) Math.floor(2 * zTmax / vD) + 3;

		ImageStack targetStack = new ImageStack(tW, tH, tD);
		IJ.log("New stack created with dimensions (" + tW + ", " + tH + ", "
				+ tD + ")");

		return targetStack;
	}

	/**
	 * Create a copy of the original image aligned to the tensor defined by a
	 * 3x3 Eigenvector matrix
	 * 
	 * @param imp
	 *            input ImagePlus stack
	 * @param E
	 *            Eigenvalue decomposition containing EigenVectors
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
	public ImagePlus alignImage(ImagePlus imp, EigenvalueDecomposition E,
			boolean doAxes, int startSlice, int endSlice, double min,
			double max, double m, double c) {
		final double[] centroid = getCentroid3D(imp, startSlice, endSlice, min,
				max, m, c);
		ImagePlus alignedImp = alignToPrincipalAxes(imp, E, centroid,
				startSlice, endSlice, min, max, doAxes);
		return alignedImp;
	}

	public boolean dialogItemChanged(GenericDialog gd, AWTEvent e) {
		Vector<?> checkboxes = gd.getCheckboxes();
		Vector<?> nFields = gd.getNumericFields();
		Checkbox box0 = (Checkbox) checkboxes.get(0);
		boolean isHUCalibrated = box0.getState();
		@SuppressWarnings("unused")
		double start = gd.getNextNumber();
		@SuppressWarnings("unused")
		double end = gd.getNextNumber();
		double min = gd.getNextNumber();
		double max = gd.getNextNumber();
		TextField minT = (TextField) nFields.get(2);
		TextField maxT = (TextField) nFields.get(3);
		if (isHUCalibrated && !fieldUpdated) {
			minT.setText("" + cal.getCValue(min));
			maxT.setText("" + cal.getCValue(max));
			fieldUpdated = true;
		}
		if (!isHUCalibrated && fieldUpdated) {
			minT.setText("" + cal.getRawValue(min));
			maxT.setText("" + cal.getRawValue(max));
			fieldUpdated = false;
		}
		if (isHUCalibrated)
			DialogModifier.replaceUnitString(gd, "grey", "HU");
		else
			DialogModifier.replaceUnitString(gd, "HU", "grey");
		return true;
	}

}
