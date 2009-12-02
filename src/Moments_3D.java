
/** Moments 3D
 *tool to calculate centroid and principal axes 
 *of a thresholded stack; assumes a 16-bit CT scan 
 *of a bone in air - default thresholds are 0 and 4000 HU
 *Eigen decomposition performed with the Jama package
 *http://math.nist.gov/javanumerics/jama/
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
 *@version 0.3
 */

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ImageProcessor;
import ij.plugin.filter.PlugInFilter;
import ij.measure.Calibration;
import ij.measure.ResultsTable;
import ij.gui.*;
import java.awt.Rectangle;

import org.doube.bonej.ImageCheck;
import org.doube.bonej.ResultInserter;
import org.doube.jama.*;

public class Moments_3D implements PlugInFilter {
	ImagePlus imp;
	protected ImageStack stack;
	public float[] CTable;
	public double[] coeff = { 0, 1 };
	public double vW, vH, vD, m, c;
	public int minT = 0, maxT = 4000; // min and maximum bone value in HU
	public int startSlice = 1, endSlice;
	public boolean doAlign, doAxes;
	public String title, units, valueUnit;

	public int setup(String arg, ImagePlus imp) {
		if (imp == null || imp.getNSlices() < 2) {
			IJ.showMessage("A stack must be open");
			return DONE;
		}
		stack = imp.getStack();
		endSlice = stack.getSize();
		this.imp = imp;
		title = imp.getShortTitle();
		// TODO extend to other image types by using
		// ImageProcessor sliceProcessor = ImageStack.getProcessor(sliceNumber);
		// int value32bit = sliceProcessor.getPixelValue(x,y);
		return DOES_16 + STACK_REQUIRED;
	}

	public void run(ImageProcessor ip) {
		if (!ImageCheck.checkIJVersion())
			return;

		Calibration cal = imp.getCalibration();
		vW = cal.pixelWidth;
		vH = cal.pixelHeight;
		vD = cal.pixelDepth;
		units = cal.getUnits();
		coeff = cal.getCoefficients();
		CTable = cal.getCTable();
		valueUnit = cal.getValueUnit();
		if (!cal.isSigned16Bit() && !cal.calibrated()) {
			IJ.run("Threshold...");
			new WaitForUserDialog(
					"This image is not density calibrated.\nSet the threshold, then click OK.")
					.show();
			minT = (int) ip.getMinThreshold();
			maxT = (int) ip.getMaxThreshold();
			showDialog(valueUnit);
			IJ.log("Image is uncalibrated: using user-determined threshold "
					+ minT + " to " + maxT);
		} else if (coeff[0] == -1000 && coeff[1] == 1.0) {
			// looks like an HU calibrated image
			// convert HU limits to pixel values
			showDialog("Hounsfield units");
			minT = (int) Math.round(cal.getRawValue(minT));
			maxT = (int) Math.round(cal.getRawValue(maxT));
			IJ.log("Image looks like it is HU calibrated. Using " + minT
					+ " and " + maxT + " " + valueUnit + " as bone cutoffs");
		} else if (cal.isSigned16Bit() && !cal.calibrated()) {
			new WaitForUserDialog(
					"This image is not density calibrated.\nSet the threshold, then click OK.")
					.show();
			minT = (int) ip.getMinThreshold();
			maxT = (int) ip.getMaxThreshold();
			showDialog(valueUnit);
			IJ.log("Image is uncalibrated: using user-determined threshold "
					+ minT + " to " + maxT);
		} else {
			new WaitForUserDialog(
					"This image is not density calibrated.\nSet the threshold, then click OK.")
					.show();
			minT = (int) ip.getMinThreshold();
			maxT = (int) ip.getMaxThreshold();
			showDialog(valueUnit);
			IJ.log("Image is uncalibrated: using user-determined threshold "
					+ minT + " to " + maxT);
			// IJ.error("Unrecognised file type");
			// return;
		}
		ResultsTable rt = ResultsTable.getResultsTable();
		rt.incrementCounter();
		rt.addLabel("Label", imp.getTitle());
		double[] centroid = findCentroid3D(stack, startSlice, endSlice);
		if (centroid[0] < 0) {
			IJ
					.error(
							"Empty Stack",
							"No voxels are a"
									+ "vailable for calculation.\nCheck your ROI and threshold.");
			return;
		}
		EigenvalueDecomposition E = calculateMoments(stack, startSlice,
				endSlice, centroid, rt);
		rt.show("Results");
		if (doAlign)
			alignToPrincipalAxes(stack, E, centroid);
		return;
	}

	// cortical bone apparent density (material density * volume fraction) from
	// Mow & Huiskes (2005) p.140
	// using 1.8 g.cm^-3: 1mm^3 = 0.0018 g = 0.0000018 kg = 1.8*10^-6 kg; 1mm^2
	// = 10^-6 m^2
	// conversion coefficient from mm^5 to kg.m^2 = 1.8*10^-12
	// double cc = 1.8*Math.pow(10, -12);
	/*-----------------------------------------------------------------*/
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

	public double voxelDensity(double pixelValue, double m, double c) {
		// density in g / cm^3 but our units are mm
		// so density is 1000* too high
		// linear calibration function, y = mx + c
		double factor = 1;
		if (units.contains("mm")) {
			factor = 1000;
		} else {
			factor = 1;
		}
		double voxelDensity = (m * pixelValue + c) / factor;
		if (voxelDensity < 0)
			voxelDensity = 0;
		return voxelDensity;
	}/* end voxelDensity */

	/*------------------------------------------------------------------*/
	/**
	 * Find the centroid of an ImageStack using startSlice and endSlice as top
	 * and bottom limits and a Rectangle as x and y limits.
	 * 
	 * @param stack
	 * @param startSlice
	 * @param endSlice
	 * @return double[3] centroid
	 * 
	 */

	public double[] findCentroid3D(ImageStack stack, int startSlice,
			int endSlice) {
		Rectangle r = stack.getRoi();
		double voxVol = vW * vH * vD;
		int offset, i;
		int w = stack.getWidth();
		int h = stack.getHeight();
		int sliceSize = w * h;
		double sumx = 0;
		double sumy = 0;
		double sumz = 0;
		double sumMass = 0;
		for (int z = startSlice; z <= endSlice; z++) {
			short[] slicePixels = new short[sliceSize];
			IJ.showStatus("Calculating centroid...");
			IJ.showProgress(z, endSlice);
			slicePixels = (short[]) stack.getPixels(z);
			for (int y = r.y; y < (r.y + r.height); y++) {
				offset = y * w;
				for (int x = r.x; x < (r.x + r.width); x++) {
					i = offset + x;
					int testPixel = slicePixels[i] & 0xffff; // convert signed
																// short to int
					if (testPixel >= minT && testPixel <= maxT) {
						double voxelMass = voxelDensity(testPixel, m, c)
								* voxVol;
						sumMass += voxelMass;
						sumx += (double) x * voxelMass;
						sumy += (double) y * voxelMass;
						sumz += (double) z * voxelMass;
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
		/*
		 * ResultsTable rt = ResultsTable.getResultsTable();
		 * rt.addValue("Xc ("+units+")", centroid[0]);
		 * rt.addValue("Yc ("+units+")", centroid[1]);
		 * rt.addValue("Zc ("+units+")", centroid[2]);
		 */
		ResultInserter ri = ResultInserter.getInstance();
		ri.setResultInRow(imp, "Xc (" + units + ")", centroid[0]);
		ri.setResultInRow(imp, "Yc (" + units + ")", centroid[1]);
		ri.setResultInRow(imp, "Zc (" + units + ")", centroid[2]);
		return centroid;
	}/* end findCentroid3D */

	public short[] makeWorkArray(ImageStack stack, int startSlice, int endSlice) {
		Rectangle r = stack.getRoi();
		int al = stack.getSize() + 1;
		int w = stack.getWidth();
		int h = stack.getHeight();
		int sliceSize = w * h;
		int sl = w * h * al;
		short[] sourceWorkArray = new short[sl];
		// for (int n = 0; n<sourceWorkArray.length; n++) sourceWorkArray[n] =
		// 0;
		for (int s = startSlice; s <= endSlice; s++) {
			short[] slicePixels = (short[]) stack.getPixels(s);
			for (int y = r.y; y < (r.y + r.height); y++) {
				int offset = y * w;
				for (int x = r.x; x < (r.x + r.width); x++) {
					int i = offset + x;
					sourceWorkArray[s * sliceSize + y * w + x] = slicePixels[i];
				}
			}
		}
		return sourceWorkArray;
	}

	public EigenvalueDecomposition calculateMoments(ImageStack stack,
			int startSlice, int endSlice, double[] centroid, ResultsTable rt) {
		// START OF 3D MOMENT CALCULATIONS
		// Our CT scans are not quantitative as they contain sharpening artefact
		// so density (g.cm^-3) is left out of this calculation and estimated at
		// the end
		Rectangle r = stack.getRoi();
		int w = stack.getWidth();
		double cX = centroid[0];
		double cY = centroid[1];
		double cZ = centroid[2];
		double voxVol = vW * vH * vD;
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
			IJ.showProgress(z, endSlice);
			short[] slicePixels = (short[]) stack.getPixels(z);
			for (int y = r.y; y < (r.y + r.height); y++) {
				int offset = y * w;
				for (int x = r.x; x < (r.x + r.width); x++) {
					int i = offset + x;
					int testPixel = slicePixels[i] & 0xffff;
					if (testPixel >= minT && testPixel <= maxT) {
						sumVoxVol += voxVol;
						double voxMass = voxelDensity(testPixel, m, c) * voxVol;
						sumVoxMass += voxMass;
						Icxx += ((y * vH - cY) * (y * vH - cY) + (z * vD - cZ)
								* (z * vD - cZ) + (vH * vH + vD * vD) / 12)
								* voxMass;
						Icyy += ((x * vW - cX) * (x * vW - cX) + (z * vD - cZ)
								* (z * vD - cZ) + (vW * vW + vD * vD) / 12)
								* voxMass;
						Iczz += ((y * vH - cY) * (y * vH - cY) + (x * vW - cX)
								* (x * vW - cX) + (vH * vH + vW * vW) / 12)
								* voxMass;
						Icxy += (x * vW - cX) * (y * vH - cY) * voxMass;
						Icxz += (x * vW - cX) * (z * vD - cZ) * voxMass;
						Icyz += (y * vH - cY) * (z * vD - cZ) * voxMass;
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

		ResultInserter ri = ResultInserter.getInstance();
		ri.setResultInRow(this.imp, "Vol (" + units + "^3)", sumVoxVol);
		ri.setResultInRow(this.imp, "Mass (g)", sumVoxMass);
		ri.setResultInRow(this.imp, "Icxx (kg.m^2)", Icxx);
		ri.setResultInRow(this.imp, "Icyy (kg.m^2)", Icyy);
		ri.setResultInRow(this.imp, "Iczz (kg.m^2)", Iczz);
		ri.setResultInRow(this.imp, "Icxy (kg.m^2)", Icxy);
		ri.setResultInRow(this.imp, "Icxz (kg.m^2)", Icxz);
		ri.setResultInRow(this.imp, "Icyz (kg.m^2)", Icyz);
		ri.setResultInRow(this.imp, "I1 (kg.m^2)", E.getD().get(2, 2));
		ri.setResultInRow(this.imp, "I2 (kg.m^2)", E.getD().get(1, 1));
		ri.setResultInRow(this.imp, "I3 (kg.m^2)", E.getD().get(0, 0));

		return E;
	}

	public void alignToPrincipalAxes(ImageStack stack,
			EigenvalueDecomposition E, double[] centroid) {
		final int al = stack.getSize() + 1;
		final int h = stack.getHeight();
		final int w = stack.getWidth();
		final double xC = centroid[0];
		final double yC = centroid[1];
		final double zC = centroid[2];
		final int sliceSize = w * h;
		final short[] sourceWorkArray = makeWorkArray(stack, startSlice, endSlice);
		short[] targetWorkArray = new short[sourceWorkArray.length];
		// Matrix eVal = E.getD();
		Matrix eVec = E.getV();
		// double[][] eigenValues = eVal.getArrayCopy();
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
		final double[][] eigenVectorsInverse = eVecInv.getArrayCopy();

		// for each voxel in the target stack, find the corresponding source
		// voxel
		for (int z = 1; z < al; z++) {
			IJ.showStatus("Aligning image stack...");
			IJ.showProgress(z, al);
			for (int y = 0; y < h; y++) {
				for (int x = 0; x < w; x++) {
					final double xD = x * vW - xC;
					final double yD = y * vH - yC;
					final double zD = z * vD - zC;
					final double xAlign = xD * eigenVectorsInverse[2][0] + yD
							* eigenVectorsInverse[1][0] + zD
							* eigenVectorsInverse[0][0] + xC;
					final double yAlign = xD * eigenVectorsInverse[2][1] + yD
							* eigenVectorsInverse[1][1] + zD
							* eigenVectorsInverse[0][1] + yC;
					final double zAlign = xD * eigenVectorsInverse[2][2] + yD
							* eigenVectorsInverse[1][2] + zD
							* eigenVectorsInverse[0][2] + zC;
					// possibility to do some voxel interpolation instead
					// of just rounding in next 3 lines
					final int xA = (int) Math.round(xAlign / vW);
					final int yA = (int) Math.round(yAlign / vH);
					final int zA = (int) Math.round(zAlign / vD);
					if (xA >= 0 && xA < w && yA >= 0 && yA < h && zA >= 0
							&& zA < al) {
						final int k = zA * sliceSize + yA * w + xA;
						// possibility of resizing the target stack depending on
						// the
						// dimensions of the aligned image rather than just
						// clipping
						targetWorkArray[z * sliceSize + y * w + x] = sourceWorkArray[k];
					}
				}
			}
		}
		if (doAxes) {
			// draw axes on stack
			int xCent = (int) Math.round(centroid[0] / vW);
			int yCent = (int) Math.round(centroid[1] / vH);
			int zCent = (int) Math.round(centroid[2] / vD);
			int centPos = (int) Math.round(yCent * w + xCent);
			short axisColour = Short.MAX_VALUE;
			for (int z = 0; z < al; z++) {
				int zOffset = z * sliceSize;
				targetWorkArray[zOffset + centPos] = axisColour;
				if (z == zCent) {
					// y axis
					for (int y = 0; y < h; y++)
						targetWorkArray[zOffset + y * w + xCent] = axisColour;
					// x axis
					for (int x = 0; x < w; x++)
						targetWorkArray[zOffset + yCent * w + x] = axisColour;
				}
			}
		}
		ImageStack target = new ImageStack(w, h);
		for (int z = 0; z < al; z++) {
			short[] sliceArray = new short[sliceSize];
			System.arraycopy(targetWorkArray, z * sliceSize, sliceArray, 0,
					sliceSize);
			int n = z + 1;
			target.addSlice("slice " + n, sliceArray);
		}
		ImagePlus impTarget = new ImagePlus(title + "_aligned", target);
		impTarget.setCalibration(imp.getCalibration());
		impTarget.setDisplayRange(imp.getDisplayRangeMin(), imp
				.getDisplayRangeMax());
		impTarget.show();
		return;
	}

	public void showDialog(String valueUnit) {
		GenericDialog gd = new GenericDialog("Setup");
		gd.addNumericField("Start Slice:", startSlice, 0);
		gd.addNumericField("End Slice:", endSlice, 0);
		gd.addMessage("Density calibration constants");
		gd.addNumericField("Slope", 0, 3, 5, "g.cm^-3 / " + valueUnit);
		gd.addNumericField("Y_Intercept", 1.8, 3, 5, "g.cm^-3");
		gd.addMessage("Only use pixels between clip values:");
		gd.addNumericField("Minimum", minT, 0, 6, valueUnit);
		gd.addNumericField("Maximum", maxT, 0, 6, valueUnit);
		gd.addCheckbox("Align result", true);
		gd.addCheckbox("Show axes", true);
		gd.showDialog();
		if (gd.wasCanceled()) {
			return;
		}
		startSlice = (int) gd.getNextNumber();
		endSlice = (int) gd.getNextNumber();
		m = gd.getNextNumber();
		c = gd.getNextNumber();
		minT = (int) gd.getNextNumber();
		maxT = (int) gd.getNextNumber();
		doAlign = gd.getNextBoolean();
		doAxes = gd.getNextBoolean();
	}
}
