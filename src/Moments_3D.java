/** Moments 3D
 *tool to calculate centroid and principal axes 
 *of a thresholded stack; originally designed for 16-bit CT scans 
 *of a bone in air so default thresholds are 0 and 4000 HU, but most greyscale images shoudl be handled
 *Outputs stack data aligned to principal axes
 *
 *Copyright 2008 2009 Michael Doube 
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
import java.awt.Rectangle;

import org.doube.bonej.ImageCheck;
import org.doube.bonej.ResultInserter;
import org.doube.jama.Matrix;
import org.doube.jama.EigenvalueDecomposition;

public class Moments_3D implements PlugIn {

	public void run(String arg) {
		if (!ImageCheck.checkIJVersion())
			return;
		final ImagePlus imp = IJ.getImage();
		if (null == imp) {
			IJ.noImage();
			return;
		}

		ImageProcessor ip = imp.getProcessor();

		Calibration cal = imp.getCalibration();
		final double[] coeff = cal.getCoefficients();
		// final float[] CTable = cal.getCTable();
		final String valueUnit = cal.getValueUnit();
		double min = 0;
		double max = 4000; // min and maximum bone value in HU
		if (!cal.isSigned16Bit() && !cal.calibrated()) {
			IJ.run("Threshold...");
			new WaitForUserDialog(
					"This image is not density calibrated.\nSet the threshold, then click OK.")
					.show();
			min = ip.getMinThreshold();
			max = ip.getMaxThreshold();
			IJ.log("Image is uncalibrated: using user-determined threshold "
					+ min + " to " + max);
		} else if (coeff[0] == -1000 && coeff[1] == 1.0) {
			// looks like an HU calibrated image
			// convert HU limits to pixel values
			min = cal.getRawValue(min);
			max = cal.getRawValue(max);
			IJ.log("Image looks like it is HU calibrated. Using " + min
					+ " and " + max + " " + valueUnit + " as bone cutoffs");
		} else if (cal.isSigned16Bit() && !cal.calibrated()) {
			new WaitForUserDialog(
					"This image is not density calibrated.\nSet the threshold, then click OK.")
					.show();
			min = ip.getMinThreshold();
			max = ip.getMaxThreshold();
			IJ.log("Image is uncalibrated: using user-determined threshold "
					+ min + " to " + max);
		} else {
			new WaitForUserDialog(
					"This image is not density calibrated.\nSet the threshold, then click OK.")
					.show();
			min = ip.getMinThreshold();
			max = ip.getMaxThreshold();
			IJ.log("Image is uncalibrated: using user-determined threshold "
					+ min + " to " + max);
		}
		GenericDialog gd = new GenericDialog("Setup");
		gd.addNumericField("Start Slice:", 1, 0);
		gd.addNumericField("End Slice:", imp.getStackSize(), 0);
		gd.addMessage("Density calibration constants");
		gd.addNumericField("Slope", 0, 3, 5, "g.cm^-3 / " + valueUnit);
		gd.addNumericField("Y_Intercept", 1.8, 3, 5, "g.cm^-3");
		gd.addMessage("Only use pixels between clip values:");
		gd.addNumericField("Minimum", min, 0, 6, valueUnit);
		gd.addNumericField("Maximum", max, 0, 6, valueUnit);
		gd.addCheckbox("Align result", true);
		gd.addCheckbox("Show axes", true);
		gd.showDialog();
		if (gd.wasCanceled()) {
			return;
		}
		final int startSlice = (int) gd.getNextNumber();
		final int endSlice = (int) gd.getNextNumber();
		final double m = gd.getNextNumber();
		final double c = gd.getNextNumber();
		min = gd.getNextNumber();
		max = gd.getNextNumber();
		final boolean doAlign = gd.getNextBoolean();
		final boolean doAxes = gd.getNextBoolean();

		double[] centroid = findCentroid3D(imp, startSlice, endSlice, min, max,
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
		ri.setResultInRow(imp, "Vol (" + units + "^3)", moments[0]);
		ri.setResultInRow(imp, "Mass (g)", moments[1]);
		ri.setResultInRow(imp, "Icxx (kg.m^2)", moments[2]);
		ri.setResultInRow(imp, "Icyy (kg.m^2)", moments[3]);
		ri.setResultInRow(imp, "Iczz (kg.m^2)", moments[4]);
		ri.setResultInRow(imp, "Icxy (kg.m^2)", moments[5]);
		ri.setResultInRow(imp, "Icxz (kg.m^2)", moments[6]);
		ri.setResultInRow(imp, "Icyz (kg.m^2)", moments[7]);
		ri.setResultInRow(imp, "I1 (kg.m^2)", E.getD().get(2, 2));
		ri.setResultInRow(imp, "I2 (kg.m^2)", E.getD().get(1, 1));
		ri.setResultInRow(imp, "I3 (kg.m^2)", E.getD().get(0, 0));

		if (doAlign)
			alignToPrincipalAxes(imp, E, centroid, startSlice, endSlice, min,
					max, doAxes).show();
		return;
	}

	// cortical bone apparent density (material density * volume fraction) from
	// Mow & Huiskes (2005) p.140
	// using 1.8 g.cm^-3: 1mm^3 = 0.0018 g = 0.0000018 kg = 1.8*10^-6 kg; 1mm^2
	// = 10^-6 m^2
	// conversion coefficient from mm^5 to kg.m^2 = 1.8*10^-12
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
	 * Get a scale factor because density is in g / cm^3 but our units are mm so
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
	private Object getEmptyPixels(int w, int h, int bitDepth) {
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
	public double[] findCentroid3D(ImagePlus imp, int startSlice, int endSlice,
			final double min, final double max, double m, double c) {
		final ImageStack stack = imp.getImageStack();
		final Rectangle r = imp.getProcessor().getRoi();
		final int rW = r.x + r.width;
		final int rH = r.y + r.height;
		final int rX = r.x;
		final int rY = r.y;
		IJ.log("findCentroid3D()");
		IJ.log("rX = "+rX);
		IJ.log("rY = "+rY);
		IJ.log("rW = "+rW);
		IJ.log("rH = "+rH);
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
		IJ.log("calculateMoments()");
		IJ.log("rX = "+rX);
		IJ.log("rY = "+rY);
		IJ.log("rW = "+rW);
		IJ.log("rH = "+rH);
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
					final double testPixel = ip.get(x, y);
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
		// create the intertia tensor matrix
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

		double[] moments = { sumVoxVol, sumVoxMass, Icxx, Icyy, Iczz, Icxy,
				Icxz, Icyz };
		Object[] result = { E, moments };

		return result;
	}

	/**
	 * 
	 * @param imp
	 * @param E
	 * @param centroid
	 * @param startSlice
	 * @param endSlice
	 * @param doAxes
	 * @return
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
		IJ.log("alignToPricipalAxes()");
		IJ.log("rX = "+rX);
		IJ.log("rY = "+rY);
		IJ.log("rW = "+rW);
		IJ.log("rH = "+rH);
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
					final int xA = (int) Math.round((xAlign + dXc) / vW);
					final int yA = (int) Math.round((yAlign + dYc) / vH);
					final int zA = (int) Math.round((zAlign + dZc) / vD);

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
			final int xCent = (int) Math.round(xTc / vW);
			final int yCent = (int) Math.round(yTc / vH);
			final int zCent = (int) Math.round(zTc / vD);
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

		ImagePlus impTarget = new ImagePlus("Aligned" + imp.getTitle(),
				targetStack);
		impTarget.setCalibration(imp.getCalibration());
		impTarget.setDisplayRange(imp.getDisplayRangeMin(), imp
				.getDisplayRangeMax());
		return impTarget;
	}

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
		IJ.log("getRotatedStack()");
		IJ.log("rX = "+rX);
		IJ.log("rY = "+rY);
		IJ.log("rW = "+rW);
		IJ.log("rH = "+rH);
		
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
						// distance from centroid in original coordinate system
						// axes
						// xCx, yCx, zCx
						final double xCx = x * vW - xC;

						// now transform each coordinate
						// transformed coordinate is dot product of original
						// coordinates
						// and eigenvectors
						final double yT= xCx * v20 + yCyv10 + zCzv00;
						final double zT= xCx * v21 + yCyv11 + zCzv01;
						final double xT= xCx * v22 + yCyv12 + zCzv02;

						// keep the biggest value
						xTmax = Math.max(xTmax, Math.abs(xT));
						yTmax = Math.max(yTmax, Math.abs(yT));
						zTmax = Math.max(zTmax, Math.abs(zT));
					}
				}
			}
		}

		int tW = (int) Math.round(2 * xTmax / vW);
		int tH = (int) Math.round(2 * yTmax / vH);
		int tD = (int) Math.round(2 * zTmax / vD);

		IJ.log("Maximum distance in x = " + xTmax);
		IJ.log("Maximum distance in y = " + yTmax);
		IJ.log("Maximum distance in z = " + zTmax);

		ImageStack targetStack = new ImageStack(tW, tH, tD);
		IJ.log("New stack created with dimensions (" + tW + ", " + tH + ", "
				+ tD + ")");

		return targetStack;
	}
}
