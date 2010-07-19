package org.doube.bonej;

/**
 * SliceGeometry plugin for ImageJ
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

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.StackConverter;
import ij.plugin.Duplicator;
import ij.plugin.PlugIn;
import ij.measure.Calibration;
import ij.measure.ResultsTable;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.gui.Wand;
import ij3d.Content;
import ij3d.Image3DUniverse;

import java.awt.AWTEvent;
import java.awt.Checkbox;
import java.awt.Color;
import java.awt.Rectangle;
import java.awt.TextField;
import java.util.ArrayList;
import java.util.List;
import java.util.Vector;

import javax.vecmath.Color3f;
import javax.vecmath.Point3f;

import org.doube.util.DialogModifier;
import org.doube.util.ImageCheck;
import org.doube.util.ThresholdGuesser;

import customnode.CustomPointMesh;

/**
 * <p>
 * Calculate 2D geometrical parameters
 * </p>
 * 
 * @author Michael Doube
 * 
 */

public class SliceGeometry implements PlugIn, DialogListener {
	private Calibration cal;
	private int al, startSlice, endSlice;
	private double vW, vH;// , min, max;
	/** Linear unit of measure */
	private String units;
	/** Do local thickness measurement in 3D */
	private boolean doThickness3D;
	/** Do local thickness measurement in 2D */
	private boolean doThickness2D;
	/** Show slice centroid */
	private boolean doCentroids;
	/** Show principal axes */
	private boolean doAxes;
	/** if true, show annotation in a new window */
	private boolean doCopy;
	/** If true, process the whole stack */
	private boolean doStack;
	/** Number of thresholded pixels in each slice */
	private double[] cslice;
	/** Cross-sectional area */
	private double[] cortArea;
	/** Mean of 3D local thickness in slice */
	private double[] meanCortThick3D;
	/** Maximum 3D local thickness in slice */
	private double[] maxCortThick3D;
	/** Standard deviation of 3D local thickness in slice */
	private double[] stdevCortThick3D;
	/** Mean of 2D local thickness in slice */
	private double[] meanCortThick2D;
	/** Maximum 2D local thickness in slice */
	private double[] maxCortThick2D;
	/** Standard deviation of 2D local thickness in slice */
	private double[] stdevCortThick2D;
	/** normal x distance from parallel axis summed over pixels */
	private double[] Sx;
	/** normal y distance from parallel axis summed over pixels */
	private double[] Sy;
	/** squared normal distances from parallel axis (Iz) */
	private double[] Sxx;
	/** squared normal distances from parallel axis (Iz) */
	private double[] Syy;
	private double[] Sxy;
	private double[] Myy;
	private double[] Mxx;
	private double[] Mxy;
	/** Angle of principal axes */
	private double[] theta;
	/**
	 * 2nd moment of area around minimum principal axis (shorter axis, larger I)
	 */
	private double[] Imax;
	/**
	 * 2nd moment of area around maximum principal axis (longer axis, smaller I)
	 */
	private double[] Imin;
	/** product moment of area, should be 0 if theta calculated perfectly */
	private double[] Ipm;
	/** length of major axis */
	private double[] R1;
	/** length of minor axis */
	private double[] R2;
	/** maximum distance from minimum principal axis (longer) */
	private double[] maxRadMin;
	/** maximum distance from maximum principal axis (shorter) */
	private double[] maxRadMax;
	/** Section modulus around minimum principal axis */
	private double[] Zmax;
	/** Section modulus around maximum principal axis */
	private double[] Zmin;
	/** Maximum diameter */
	private double[] feretMax;
	/** Angle of maximum diameter */
	private double[] feretAngle;
	/** Minimum diameter */
	private double[] feretMin;
	/** List of empty slices. If true, slice contains 0 pixels to analyse */
	private boolean[] emptySlices;
	/** List of slice centroids */
	private double[][] sliceCentroids;
	private double[] integratedDensity;
	private double[] meanDensity;
	private double m;
	private double c;
	private double[][] weightedCentroids;
	private boolean fieldUpdated = false;
	/** List of perimeter lengths */
	private double[] perimeter;
	/** List of maximal distances from centroid */
	private double[] maxRadCentre;
	/** List of polar section moduli */
	private double[] Zpol;
	private boolean do3DAnnotation;

	public void run(String arg) {
		if (!ImageCheck.checkEnvironment())
			return;
		ImagePlus imp = IJ.getImage();
		if (null == imp) {
			IJ.noImage();
			return;
		}

		this.cal = imp.getCalibration();
		this.vW = cal.pixelWidth;
		this.vH = cal.pixelHeight;
		this.units = cal.getUnits();
		this.al = imp.getStackSize() + 1;

		String pixUnits;
		if (ImageCheck.huCalibrated(imp)) {
			pixUnits = "HU";
			fieldUpdated = true;
		} else
			pixUnits = "grey";

		double[] thresholds = ThresholdGuesser.setDefaultThreshold(imp);
		double min = thresholds[0];
		double max = thresholds[1];
		GenericDialog gd = new GenericDialog("Options");

		// guess bone from image title
		int boneID = BoneList.guessBone(imp);
		String[] bones = BoneList.get();
		gd.addChoice("Bone: ", bones, bones[boneID]);

		gd.addCheckbox("2D_Thickness", true);
		gd.addCheckbox("3D_Thickness", false);
		gd.addCheckbox("Draw_Axes", true);
		gd.addCheckbox("Draw_Centroids", true);
		gd.addCheckbox("Annotated_Copy_(2D)", true);
		gd.addCheckbox("3D_Annotation", false);
		gd.addCheckbox("Process_Stack", false);
		// String[] analyses = { "Weighted", "Unweighted", "Both" };
		// gd.addChoice("Calculate: ", analyses, analyses[1]);
		gd.addCheckbox("HU_Calibrated", ImageCheck.huCalibrated(imp));
		gd.addNumericField("Bone_Min:", min, 1, 6, pixUnits + " ");
		gd.addNumericField("Bone_Max:", max, 1, 6, pixUnits + " ");
		gd
				.addMessage("Only pixels >= bone min\n"
						+ "and <= bone max are used.");
		gd.addMessage("Density calibration coefficients");
		gd.addNumericField("Slope", 0, 4, 6, "g.cm^-3 / " + pixUnits + " ");
		gd.addNumericField("Y_Intercept", 1.8, 4, 6, "g.cm^-3");
		gd.addHelp("http://bonej.org/slice");
		gd.addDialogListener(this);
		gd.showDialog();
		String bone = gd.getNextChoice();
		boneID = BoneList.guessBone(bone);
		this.doThickness2D = gd.getNextBoolean();
		this.doThickness3D = gd.getNextBoolean();
		this.doAxes = gd.getNextBoolean();
		this.doCentroids = gd.getNextBoolean();
		this.doCopy = gd.getNextBoolean();
		this.do3DAnnotation = gd.getNextBoolean();
		this.doStack = gd.getNextBoolean();
		if (this.doStack) {
			this.startSlice = 1;
			this.endSlice = imp.getImageStackSize();
		} else {
			this.startSlice = imp.getCurrentSlice();
			this.endSlice = imp.getCurrentSlice();
		}

		boolean isHUCalibrated = gd.getNextBoolean();
		min = gd.getNextNumber();
		max = gd.getNextNumber();
		if (isHUCalibrated) {
			min = cal.getRawValue(min);
			max = cal.getRawValue(max);
		}
		this.m = gd.getNextNumber();
		this.c = gd.getNextNumber();
		if (gd.wasCanceled())
			return;

		if (calculateCentroids(imp, min, max) == 0) {
			IJ.error("No pixels available to calculate.\n"
					+ "Please check the threshold and ROI.");
			return;
		}

		calculateMoments(imp, min, max);
		if (this.doThickness3D)
			calculateThickness3D(imp, min, max);
		if (this.doThickness2D)
			calculateThickness2D(imp, min, max);

		roiMeasurements(imp, min, max);

		// TODO locate centroids of multiple sections in a single plane

		ResultsTable rt = ResultsTable.getResultsTable();
		rt.reset();

		String title = imp.getTitle();
		for (int s = this.startSlice; s <= this.endSlice; s++) {
			rt.incrementCounter();
			rt.addLabel(title);
			rt.addValue("Bone Code", boneID);
			rt.addValue("Slice", s);
			rt.addValue("CSA (" + units + "²)", this.cortArea[s]);
			rt.addValue("X cent. (" + units + ")", this.sliceCentroids[0][s]);
			rt.addValue("Y cent. (" + units + ")", this.sliceCentroids[1][s]);
			rt.addValue("Density", this.meanDensity[s]);
			rt.addValue("wX cent. (" + units + ")",
					this.weightedCentroids[0][s]);
			rt.addValue("wY cent. (" + units + ")",
					this.weightedCentroids[1][s]);
			rt.addValue("Theta (rad)", this.theta[s]);
			rt.addValue("R1 (" + units + ")", this.R1[s]);
			rt.addValue("R2 (" + units + ")", this.R2[s]);
			rt.addValue("Imin (" + units + "^4)", this.Imin[s]);
			rt.addValue("Imax (" + units + "^4)", this.Imax[s]);
			rt.addValue("Ipm (" + units + "^4)", this.Ipm[s]);
			rt.addValue("Zmax (" + units + "³)", this.Zmax[s]);
			rt.addValue("Zmin (" + units + "³)", this.Zmin[s]);
			rt.addValue("Zpol (" + units + "³)", this.Zpol[s]);
			rt.addValue("Feret Min (" + units + ")", this.feretMin[s]);
			rt.addValue("Feret Max (" + units + ")", this.feretMax[s]);
			rt.addValue("Feret Angle (rad)", this.feretAngle[s]);
			rt.addValue("Perimeter (" + units + ")", this.perimeter[s]);
			if (this.doThickness3D) {
				rt.addValue("Max Thick 3D (" + units + ")",
						this.maxCortThick3D[s]);
				rt.addValue("Mean Thick 3D (" + units + ")",
						this.meanCortThick3D[s]);
				rt.addValue("SD Thick 3D (" + units + ")",
						this.stdevCortThick3D[s]);
			}
			if (this.doThickness2D) {
				rt.addValue("Max Thick 2D (" + units + ")",
						this.maxCortThick2D[s]);
				rt.addValue("Mean Thick 2D (" + units + ")",
						this.meanCortThick2D[s]);
				rt.addValue("SD Thick 2D (" + units + ")",
						this.stdevCortThick2D[s]);
			}
		}
		rt.show("Results");

		if (this.doAxes || this.doCentroids) {
			if (!this.doCopy) {
				ImagePlus annImp = annotateImage(imp);
				imp.setStack(null, annImp.getImageStack());
			} else {
				annotateImage(imp).show();
			}
		}
		if (this.do3DAnnotation)
			show3DAxes(imp);
		return;
	}

	/**
	 * Draw centroids and / or principal axes on a copy of the original image
	 * 
	 * @param imp
	 * @return ImagePlus with centroid and / or principal axes drawn
	 */
	private ImagePlus annotateImage(ImagePlus imp) {
		ImageStack stack = imp.getImageStack();
		int w = stack.getWidth();
		int h = stack.getHeight();
		ImageStack annStack = new ImageStack(w, h);
		for (int s = this.startSlice; s <= this.endSlice; s++) {
			ImageProcessor annIP = stack.getProcessor(s).duplicate();
			annIP.setColor(Color.white);
			double cX = this.sliceCentroids[0][s] / this.vW;
			double cY = this.sliceCentroids[1][s] / this.vH;

			if (this.doCentroids && !this.emptySlices[s]) {
				annIP.drawOval((int) Math.floor(cX - 4), (int) Math
						.floor(cY - 4), 8, 8);
			}

			if (this.doAxes && !this.emptySlices[s]) {
				double th = this.theta[s];
				double rMin = this.R1[s];
				double rMax = this.R2[s];
				double thPi = th + Math.PI / 2;

				int x1 = (int) Math.floor(cX - Math.cos(thPi) * 2 * rMin);
				int y1 = (int) Math.floor(cY - Math.sin(thPi) * 2 * rMin);
				int x2 = (int) Math.floor(cX + Math.cos(thPi) * 2 * rMin);
				int y2 = (int) Math.floor(cY + Math.sin(thPi) * 2 * rMin);
				annIP.drawLine(x1, y1, x2, y2);

				x1 = (int) Math.floor(cX - Math.cos(-th) * 2 * rMax);
				y1 = (int) Math.floor(cY + Math.sin(-th) * 2 * rMax);
				x2 = (int) Math.floor(cX + Math.cos(-th) * 2 * rMax);
				y2 = (int) Math.floor(cY - Math.sin(-th) * 2 * rMax);
				annIP.drawLine(x1, y1, x2, y2);
			}
			annStack.addSlice(stack.getSliceLabel(s), annIP);
		}
		ImagePlus ann = new ImagePlus("Annotated_" + imp.getTitle(), annStack);
		ann.setCalibration(imp.getCalibration());
		return ann;
	}

	/**
	 * Display principal axes on a 3D rendered version of the image
	 * 
	 * @param imp
	 *            Original image
	 */
	private void show3DAxes(ImagePlus imp) {
		Calibration cal = imp.getCalibration();
		// copy the data from inside the ROI and convert it to 8-bit
		Duplicator d = new Duplicator();
		ImagePlus roiImp = d.run(imp, 1, imp.getImageStackSize());

		// initialise and show the 3D universe
		Image3DUniverse univ = new Image3DUniverse();
		univ.show();

		double rX = 0;
		double rY = 0;
		if (imp.getRoi() != null) {
			Rectangle roi = imp.getRoi().getBounds();
			rX = roi.getX() * cal.pixelWidth;
			rY = roi.getY() * cal.pixelHeight;
		}

		// list of centroids
		List<Point3f> centroids = new ArrayList<Point3f>();
		// list of axes
		List<Point3f> axes = new ArrayList<Point3f>();
		for (int s = 1; s <= roiImp.getImageStackSize(); s++) {
			if (((Double)this.cortArea[s]).equals(Double.NaN))
				continue;

			final double cX = sliceCentroids[0][s] - rX;
			final double cY = sliceCentroids[1][s] - rY;
			final double cZ = (s - 0.5) * cal.pixelDepth;

			Point3f cent = new Point3f();
			cent.x = (float) cX;
			cent.y = (float) cY;
			cent.z = (float) cZ;
			centroids.add(cent);

			// add the axes to the list
			double th = this.theta[s];
			double rMin = this.R1[s] * cal.pixelWidth;
			double rMax = this.R2[s] * cal.pixelWidth;
			double thPi = th + Math.PI / 2;

			Point3f start1 = new Point3f();
			start1.x = (float) (cX - Math.cos(thPi) * 2 * rMin);
			start1.y = (float) (cY - Math.sin(thPi) * 2 * rMin);
			start1.z = (float) cZ;
			axes.add(start1);

			Point3f end1 = new Point3f();
			end1.x = (float) (cX + Math.cos(thPi) * 2 * rMin);
			end1.y = (float) (cY + Math.sin(thPi) * 2 * rMin);
			end1.z = (float) cZ;
			axes.add(end1);

			Point3f start2 = new Point3f();
			start2.x = (float) (cX - Math.cos(-th) * 2 * rMax);
			start2.y = (float) (cY + Math.sin(-th) * 2 * rMax);
			start2.z = (float) cZ;
			axes.add(start2);

			Point3f end2 = new Point3f();
			end2.x = (float) (cX + Math.cos(-th) * 2 * rMax);
			end2.y = (float) (cY - Math.sin(-th) * 2 * rMax);
			end2.z = (float) cZ;
			axes.add(end2);
		}
		// show the centroids
		CustomPointMesh mesh = new CustomPointMesh(centroids);
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

	/**
	 * Calculate the centroid of each slice
	 * 
	 * @param imp
	 *            Input image
	 * @return double containing sum of pixel count
	 */
	private double calculateCentroids(ImagePlus imp, double min, double max) {
		ImageStack stack = imp.getImageStack();
		Rectangle r = stack.getRoi();
		// 2D centroids
		this.sliceCentroids = new double[2][this.al];
		// pixel counters
		double cstack = 0;
		this.emptySlices = new boolean[this.al];
		this.cslice = new double[this.al];
		this.cortArea = new double[this.al];
		this.integratedDensity = new double[this.al];
		this.meanDensity = new double[this.al];
		this.weightedCentroids = new double[2][this.al];
		final double pixelArea = this.vW * this.vH;
		final int roiXEnd = r.x + r.width;
		final int roiYEnd = r.y + r.height;
		for (int s = this.startSlice; s <= this.endSlice; s++) {
			IJ.showStatus("Calculating centroids...");
			IJ.showProgress(s - this.startSlice, this.endSlice);
			double sumX = 0;
			double sumY = 0;
			int count = 0;
			double sumD = 0;
			double wSumX = 0;
			double wSumY = 0;
			ImageProcessor ip = stack.getProcessor(s);
			for (int y = r.y; y < roiYEnd; y++) {
				for (int x = r.x; x < roiXEnd; x++) {
					final double pixel = (double) ip.get(x, y);
					if (pixel >= min && pixel <= max) {
						count++;
						sumX += x;
						sumY += y;
						final double wP = pixel * this.m + this.c;
						sumD += wP;
						wSumX += x * wP;
						wSumY += y * wP;
					}
				}
			}
			this.cslice[s] = count;
			this.cortArea[s] = count * pixelArea;
			if (count > 0) {
				this.sliceCentroids[0][s] = sumX * this.vW / count;
				this.sliceCentroids[1][s] = sumY * this.vH / count;
				this.integratedDensity[s] = sumD;
				this.meanDensity[s] = sumD / count;
				this.weightedCentroids[0][s] = wSumX * this.vW / sumD;
				this.weightedCentroids[1][s] = wSumY * this.vH / sumD;
				cstack += count;
				this.emptySlices[s] = false;
			} else {
				this.emptySlices[s] = true;
				this.cortArea[s] = Double.NaN;
				this.sliceCentroids[0][s] = Double.NaN;
				this.sliceCentroids[1][s] = Double.NaN;
				this.cslice[s] = Double.NaN;
			}
		}
		return cstack;
	}

	/**
	 * Calculate second moments of area, length and angle of principal axes
	 * 
	 * @param imp
	 */
	private void calculateMoments(ImagePlus imp, double min, double max) {
		final ImageStack stack = imp.getImageStack();
		final Rectangle r = stack.getRoi();
		// START OF Ix AND Iy CALCULATION
		this.Sx = new double[this.al];
		this.Sy = new double[this.al];
		this.Sxx = new double[this.al];
		this.Syy = new double[this.al];
		this.Sxy = new double[this.al];
		this.Myy = new double[this.al];
		this.Mxx = new double[this.al];
		this.Mxy = new double[this.al];
		this.theta = new double[this.al];
		for (int s = this.startSlice; s <= this.endSlice; s++) {
			IJ.showStatus("Calculating Ix and Iy...");
			IJ.showProgress(s, this.endSlice);
			double sxs = 0;
			double sys = 0;
			double sxxs = 0;
			double syys = 0;
			double sxys = 0;
			final int roiXEnd = r.x + r.width;
			final int roiYEnd = r.y + r.height;
			if (!this.emptySlices[s]) {
				ImageProcessor ip = stack.getProcessor(s);
				for (int y = r.y; y < roiYEnd; y++) {
					for (int x = r.x; x < roiXEnd; x++) {
						final double pixel = (double) ip.get(x, y);
						if (pixel >= min && pixel <= max) {
							final double xVw = x * vW;
							final double yVh = y * vH;
							sxs += xVw;
							sys += yVh;
							sxxs += xVw * xVw;
							syys += yVh * yVh;
							sxys += xVw * yVh;
						}
					}
				}
				this.Sx[s] = sxs;
				this.Sy[s] = sys;
				this.Sxx[s] = sxxs;
				this.Syy[s] = syys;
				this.Sxy[s] = sxys;
				this.Myy[s] = this.Sxx[s]
						- (this.Sx[s] * this.Sx[s] / this.cslice[s])
						+ this.cslice[s] * vW * vW / 12;
				// this.cslice[]/12 is for each pixel's own moment
				this.Mxx[s] = this.Syy[s]
						- (this.Sy[s] * this.Sy[s] / this.cslice[s])
						+ this.cslice[s] * vH * vH / 12;
				this.Mxy[s] = this.Sxy[s]
						- (this.Sx[s] * this.Sy[s] / this.cslice[s])
						+ this.cslice[s] * vH * vW / 12;
				if (this.Mxy[s] == 0)
					this.theta[s] = 0;
				else {
					this.theta[s] = Math.atan((this.Mxx[s] - this.Myy[s] + Math
							.sqrt((this.Mxx[s] - this.Myy[s])
									* (this.Mxx[s] - this.Myy[s]) + 4
									* this.Mxy[s] * this.Mxy[s]))
							/ (2 * this.Mxy[s]));
				}
			} else {
				this.theta[s] = Double.NaN;
			}
		}
		// END OF Ix and Iy CALCULATION
		// START OF Imax AND Imin CALCULATION
		this.Imax = new double[this.al];
		this.Imin = new double[this.al];
		this.Ipm = new double[this.al];
		this.R1 = new double[this.al];
		this.R2 = new double[this.al];
		this.maxRadMin = new double[this.al];
		this.maxRadMax = new double[this.al];
		this.maxRadCentre = new double[this.al];
		this.Zmax = new double[this.al];
		this.Zmin = new double[this.al];
		this.Zpol = new double[this.al];
		for (int s = this.startSlice; s <= this.endSlice; s++) {
			IJ.showStatus("Calculating Imin and Imax...");
			IJ.showProgress(s, this.endSlice);
			if (!this.emptySlices[s]) {
				ImageProcessor ip = stack.getProcessor(s);
				double sxs = 0;
				double sys = 0;
				double sxxs = 0;
				double syys = 0;
				double sxys = 0;
				double maxRadMinS = 0;
				double maxRadMaxS = 0;
				double maxRadCentreS = 0;
				final double cosTheta = Math.cos(this.theta[s]);
				final double sinTheta = Math.sin(this.theta[s]);
				final int roiYEnd = r.y + r.height;
				final int roiXEnd = r.x + r.width;
				final double xC = this.sliceCentroids[0][s];
				final double yC = this.sliceCentroids[1][s];
				final double cS = this.cslice[s];
				for (int y = r.y; y < roiYEnd; y++) {
					final double yYc = y * vH - yC;
					for (int x = r.x; x < roiXEnd; x++) {
						final double pixel = (double) ip.get(x, y);
						if (pixel >= min && pixel <= max) {
							final double xXc = x * vW - xC;
							final double xCosTheta = x * vW * cosTheta;
							final double yCosTheta = y * vH * cosTheta;
							final double xSinTheta = x * vW * sinTheta;
							final double ySinTheta = y * vH * sinTheta;
							sxs += xCosTheta + ySinTheta;
							sys += yCosTheta - xSinTheta;
							sxxs += (xCosTheta + ySinTheta)
									* (xCosTheta + ySinTheta);
							syys += (yCosTheta - xSinTheta)
									* (yCosTheta - xSinTheta);
							sxys += (yCosTheta - xSinTheta)
									* (xCosTheta + ySinTheta);
							maxRadMinS = Math.max(maxRadMinS, Math.abs(xXc
									* cosTheta + yYc * sinTheta));
							maxRadMaxS = Math.max(maxRadMaxS, Math.abs(yYc
									* cosTheta - xXc * sinTheta));
							maxRadCentreS = Math.max(maxRadCentreS, Math
									.sqrt(xXc * xXc + yYc * yYc));
						}
					}
				}
				this.Sx[s] = sxs;
				this.Sy[s] = sys;
				this.Sxx[s] = sxxs;
				this.Syy[s] = syys;
				this.Sxy[s] = sxys;
				this.maxRadMin[s] = maxRadMinS;
				this.maxRadMax[s] = maxRadMaxS;
				this.maxRadCentre[s] = maxRadCentreS;
				final double pixelMoments = cS * vW * vH
						* (cosTheta * cosTheta + sinTheta * sinTheta) / 12;
				this.Imax[s] = vW
						* vH
						* (this.Sxx[s] - (this.Sx[s] * this.Sx[s] / cS) + pixelMoments);
				this.Imin[s] = vW
						* vH
						* (this.Syy[s] - (this.Sy[s] * this.Sy[s] / cS) + pixelMoments);
				this.Ipm[s] = this.Sxy[s] - (this.Sy[s] * this.Sx[s] / cS)
						+ pixelMoments;
				this.R1[s] = Math.sqrt(this.Imin[s] / (cS * vW * vH * vW * vH));
				this.R2[s] = Math.sqrt(this.Imax[s] / (cS * vW * vH * vW * vH));
				this.Zmax[s] = this.Imax[s] / this.maxRadMin[s];
				this.Zmin[s] = this.Imin[s] / this.maxRadMax[s];
				this.Zpol[s] = (this.Imax[s] + this.Imin[s])
						/ this.maxRadCentre[s];
			} else {
				this.Imax[s] = Double.NaN;
				this.Imin[s] = Double.NaN;
				this.Ipm[s] = Double.NaN;
				this.R1[s] = Double.NaN;
				this.R2[s] = Double.NaN;
				this.maxRadMin[s] = Double.NaN;
				this.maxRadMax[s] = Double.NaN;
				this.Zmax[s] = Double.NaN;
				this.Zmin[s] = Double.NaN;
				this.Zpol[s] = Double.NaN;
			}
		}
		return;
	}

	/**
	 * Calculate 3D Local Thickness and determine thickness statistics for the
	 * slice
	 * 
	 */
	private void calculateThickness3D(ImagePlus imp, double min, double max) {
		this.maxCortThick3D = new double[this.al];
		this.meanCortThick3D = new double[this.al];
		this.stdevCortThick3D = new double[this.al];
		Rectangle r = imp.getProcessor().getRoi();
		Thickness th = new Thickness();

		// convert to binary
		ImagePlus binaryImp = convertToBinary(imp, min, max);

		ImagePlus thickImp = th.getLocalThickness(binaryImp, false);

		for (int s = this.startSlice; s <= this.endSlice; s++) {
			if (this.emptySlices[s]) {
				this.maxCortThick3D[s] = Double.NaN;
				this.meanCortThick3D[s] = Double.NaN;
				this.stdevCortThick3D[s] = Double.NaN;
				continue;
			}
			FloatProcessor ip = (FloatProcessor) thickImp.getStack()
					.getProcessor(s);
			double sumPix = 0;
			double sliceMax = 0;
			double pixCount = 0;
			final int roiXEnd = r.x + r.width;
			final int roiYEnd = r.y + r.height;
			for (int y = r.y; y < roiYEnd; y++) {
				for (int x = r.x; x < roiXEnd; x++) {
					final float pixel = Float.intBitsToFloat(ip.get(x, y));
					if (pixel > 0) {
						pixCount++;
						sumPix += pixel;
						sliceMax = Math.max(sliceMax, pixel);
					}
				}
			}
			final double sliceMean = sumPix / pixCount;
			this.meanCortThick3D[s] = sliceMean;
			this.maxCortThick3D[s] = sliceMax;

			double sumSquares = 0;
			for (int y = r.y; y < roiYEnd; y++) {
				for (int x = r.x; x < roiXEnd; x++) {
					final float pixel = Float.intBitsToFloat(ip.get(x, y));
					if (pixel > 0) {
						final double d = sliceMean - pixel;
						sumSquares += d * d;
					}
				}
			}
			this.stdevCortThick3D[s] = Math.sqrt(sumSquares / pixCount);
		}
		return;
	}

	/**
	 * Calculate thickness on individual slices using local thickness
	 * 
	 * @param imp
	 */
	private void calculateThickness2D(ImagePlus imp, double min, double max) {
		this.maxCortThick2D = new double[this.al];
		this.meanCortThick2D = new double[this.al];
		this.stdevCortThick2D = new double[this.al];

		int nThreads = Runtime.getRuntime().availableProcessors();
		SliceThread[] sliceThread = new SliceThread[nThreads];
		for (int thread = 0; thread < nThreads; thread++) {
			sliceThread[thread] = new SliceThread(thread, nThreads, imp, min,
					max, this.meanCortThick2D, this.maxCortThick2D,
					this.stdevCortThick2D, this.startSlice, this.endSlice,
					this.emptySlices);
			sliceThread[thread].start();
		}
		try {
			for (int thread = 0; thread < nThreads; thread++) {
				sliceThread[thread].join();
			}
		} catch (InterruptedException ie) {
			IJ.error("A thread was interrupted.");
		}
		return;
	}

	class SliceThread extends Thread {
		final int thread, nThreads, width, height, startSlice, endSlice;

		double min, max;

		double[] meanThick, maxThick, stdevThick;

		boolean[] emptySlices;

		final ImagePlus impT;

		public SliceThread(int thread, int nThreads, ImagePlus imp, double min,
				double max, double[] meanThick, double[] maxThick,
				double[] stdevThick, int startSlice, int endSlice,
				boolean[] emptySlices) {
			this.impT = imp;
			this.min = min;
			this.max = max;
			this.width = this.impT.getWidth();
			this.height = this.impT.getHeight();
			this.thread = thread;
			this.nThreads = nThreads;
			this.meanThick = meanThick;
			this.maxThick = maxThick;
			this.stdevThick = stdevThick;
			this.startSlice = startSlice;
			this.endSlice = endSlice;
			this.emptySlices = emptySlices;
		}

		public void run() {
			for (int s = this.thread + this.startSlice; s <= this.endSlice; s += this.nThreads) {
				if (this.emptySlices[s]) {
					this.meanThick[s] = Double.NaN;
					this.maxThick[s] = Double.NaN;
					this.stdevThick[s] = Double.NaN;
					continue;
				}
				ImageProcessor ip = impT.getImageStack().getProcessor(s);
				ImagePlus sliceImp = new ImagePlus(" " + s, ip);
				Rectangle r = ip.getRoi();
				// binarise
				ImagePlus binaryImp = convertToBinary(sliceImp, min, max);
				Calibration cal = impT.getCalibration();
				binaryImp.setCalibration(cal);
				// calculate thickness
				Thickness th = new Thickness();
				ImagePlus thickImp = th.getLocalThickness(binaryImp, false);
				FloatProcessor thickIp = (FloatProcessor) thickImp
						.getProcessor();
				double sumPix = 0;
				double sliceMax = 0;
				double pixCount = 0;
				final double roiXEnd = r.x + r.width;
				final double roiYEnd = r.y + r.height;
				for (int y = r.y; y < roiYEnd; y++) {
					for (int x = r.x; x < roiXEnd; x++) {
						final float pixel = Float.intBitsToFloat(thickIp.get(x,
								y));
						if (pixel > 0) {
							pixCount++;
							sumPix += pixel;
							sliceMax = Math.max(sliceMax, pixel);
						}
					}
				}
				final double sliceMean = sumPix / pixCount;
				this.meanThick[s] = sliceMean;
				this.maxThick[s] = sliceMax;

				double sumSquares = 0;
				for (int y = r.y; y < roiYEnd; y++) {
					for (int x = r.x; x < roiXEnd; x++) {
						final float pixel = Float.intBitsToFloat(thickIp.get(x,
								y));
						if (pixel > 0) {
							final double d = sliceMean - pixel;
							sumSquares += d * d;
						}
					}
				}
				this.stdevThick[s] = Math.sqrt(sumSquares / pixCount);
			}
			return;
		}
	}

	private ImagePlus convertToBinary(ImagePlus imp, double min, double max) {
		final int w = imp.getWidth();
		final int h = imp.getHeight();
		final int d = imp.getStackSize();
		final ImageStack sourceStack = imp.getImageStack();
		ImageStack binaryStack = new ImageStack(w, h);
		for (int s = 1; s <= d; s++) {
			ImageProcessor sliceIp = sourceStack.getProcessor(s);
			ByteProcessor binaryIp = new ByteProcessor(w, h);
			for (int y = 0; y < h; y++) {
				for (int x = 0; x < w; x++) {
					if (sliceIp.get(x, y) >= min && sliceIp.get(x, y) <= max) {
						binaryIp.set(x, y, 255);
					} else {
						binaryIp.set(x, y, 0);
					}
				}
			}
			binaryStack.addSlice(sourceStack.getSliceLabel(s), binaryIp);
		}
		ImagePlus binaryImp = new ImagePlus("binaryImp", binaryStack);
		binaryImp.setCalibration(imp.getCalibration());
		return binaryImp;
	}

	private void roiMeasurements(ImagePlus imp, double min, double max) {
		Roi initialRoi = imp.getRoi();
		double[] feretValues = new double[3];
		this.feretAngle = new double[this.al];
		this.feretMax = new double[this.al];
		this.feretMin = new double[this.al];
		this.perimeter = new double[this.al];
		int initialSlice = imp.getCurrentSlice();
		// for the required slices...
		for (int s = this.startSlice; s <= this.endSlice; s++) {
			ImageProcessor ip = imp.getImageStack().getProcessor(s);
			Wand w = new Wand(ip);
			w.autoOutline(0, (int) Math.round(this.sliceCentroids[1][s]
					/ this.vH), min, max, Wand.EIGHT_CONNECTED);
			if (this.emptySlices[s] || w.npoints == 0) {
				this.feretMin[s] = Double.NaN;
				this.feretAngle[s] = Double.NaN;
				this.feretMax[s] = Double.NaN;
				this.perimeter[s] = Double.NaN;
			} else {
				int type = Wand.allPoints() ? Roi.FREEROI : Roi.TRACED_ROI;
				PolygonRoi roi = new PolygonRoi(w.xpoints, w.ypoints,
						w.npoints, type);
				feretValues = roi.getFeretValues();
				this.feretMin[s] = feretValues[2] * this.vW;
				this.feretAngle[s] = feretValues[1] * Math.PI / 180;
				this.feretMax[s] = feretValues[0] * this.vW;
				this.perimeter[s] = roi.getLength() * this.vW;
			}
			feretValues = null;
		}
		IJ.setSlice(initialSlice);
		imp.setRoi(initialRoi);
		return;
	}

	public boolean dialogItemChanged(GenericDialog gd, AWTEvent e) {
		Vector<?> checkboxes = gd.getCheckboxes();
		Vector<?> nFields = gd.getNumericFields();
		Checkbox box6 = (Checkbox) checkboxes.get(7);
		boolean isHUCalibrated = box6.getState();
		TextField minT = (TextField) nFields.get(0);
		TextField maxT = (TextField) nFields.get(1);
		double min = 0;
		double max = 0;
		try {
			min = Double.parseDouble(minT.getText());
			max = Double.parseDouble(maxT.getText());
		} catch (Exception ex) {
			IJ.error("You put text in a number field");
		}
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
		DialogModifier.registerMacroValues(gd, gd.getComponents());
		return true;
	}
}
