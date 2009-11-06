/**
 * Slice_Geometry plugin for ImageJ
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

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;
import ij.plugin.PlugIn;
import ij.measure.Calibration;
import ij.measure.ResultsTable;
import ij.gui.*;

import java.awt.Color;
import java.awt.Rectangle;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.doube.bonej.BoneList;
import org.doube.bonej.Thickness;

/**
 * <p>
 * Calculate 2D geometrical parameters
 * </p>
 * 
 * @author Michael Doube
 * 
 */

public class Slice_Geometry implements PlugIn {
	private int boneID, al, startSlice, endSlice;
	private double vW, vH, vD, airHU, min, max;
	private String units, analyse, calString;

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
	private boolean isCalibrated;
	private double[] cslice;
	private double[] cortArea;
	private double[] meanCortThick3D;
	private double[] maxCortThick3D;
	private double[] stdevCortThick3D;
	private double[] meanCortThick2D;
	private double[] maxCortThick2D;
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
	/** Section modulus around maximum principal axis */
	private double[] Zmax;
	/** Section modulus around minimum principal axis */
	private double[] Zmin;
	/** Alternative calculation of Imax */
	private double[] ImaxFast;
	/** Alternative calculation of Imin */
	private double[] IminFast;
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
	Calibration cal;

	public void run(String arg) {
		ImagePlus imp = IJ.getImage();
		if (null == imp) {
			IJ.noImage();
			return;
		}

		this.cal = imp.getCalibration();
		this.vW = this.cal.pixelWidth;
		this.vH = this.cal.pixelHeight;
		this.vD = this.cal.pixelDepth;
		this.units = this.cal.getUnits();
		this.al = imp.getStackSize() + 1;

		IJ.run("Threshold...");
		new WaitForUserDialog("Adjust the threshold, then click OK.").show();
		this.min = (short) imp.getProcessor().getMinThreshold();
		this.max = (short) imp.getProcessor().getMaxThreshold();

		if (!showDialog(imp)) {
			return;
		}

		if (calculateCentroids(imp) == 0) {
			IJ.error("No pixels available to calculate.\n"
					+ "Please check the threshold and ROI.");
			return;
		}

		calculateMoments(imp);
		if (this.doThickness3D)
			calculateThickness3D(imp);
		if (this.doThickness2D)
			calculateThickness2D(imp);

		roiMeasurements(imp);
		// TODO locate centroids of multiple sections in a single plane

		showSliceResults(imp);

		if (this.doAxes || this.doCentroids) {
			if (!this.doCopy) {
				ImagePlus annImp = annotateImage(imp);
				imp.setStack(null, annImp.getImageStack());
			} else {
				annotateImage(imp).show();
			}
		}
	}

	/**
	 * Draw centroids and / or principal axes on a copy of the original image
	 * 
	 * @param imp
	 * @return
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

			if (this.doCentroids) {
				annIP.drawOval((int) Math.floor(cX - 4), (int) Math
						.floor(cY - 4), 8, 8);
			}

			if (this.doAxes) {
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
				annStack.addSlice(stack.getSliceLabel(s), annIP);
			}
		}
		ImagePlus ann = new ImagePlus("Annotated_" + imp.getTitle(), annStack);
		ann.setCalibration(imp.getCalibration());
		return ann;
	}

	protected double calculateCentroids(ImagePlus imp) {
		ImageStack stack = imp.getImageStack();
		// 2D centroids
		this.sliceCentroids = new double[2][this.al];
		// pixel counters
		double cstack = 0;
		Rectangle r = stack.getRoi();
		int w = stack.getWidth();
		this.emptySlices = new boolean[this.al];
		this.cslice = new double[this.al];
		this.cortArea = new double[this.al];
		double pixelArea = this.vW * this.vH;
		for (int s = this.startSlice; s <= this.endSlice; s++) {
			double sumX = 0;
			double sumY = 0;
			this.cslice[s] = 0;
			short[] pixels = (short[]) stack.getPixels(s);
			for (int y = r.y; y < (r.y + r.height); y++) {
				int offset = y * w;
				for (int x = r.x; x < (r.x + r.width); x++) {
					int i = offset + x;
					if (pixels[i] >= this.min && pixels[i] <= this.max) {
						this.cslice[s]++;
						this.cortArea[s] += pixelArea;
						sumX += x * this.vW;
						sumY += y * this.vH;
					}
				}
			}
			if (this.cslice[s] > 0) {
				this.sliceCentroids[0][s] = sumX / this.cslice[s];
				this.sliceCentroids[1][s] = sumY / this.cslice[s];
				cstack += this.cslice[s];
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

	private void calculateMoments(ImagePlus imp) {
		ImageStack stack = imp.getImageStack();
		Rectangle r = stack.getRoi();
		int w = stack.getWidth();
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
		for (int s = 1; s <= stack.getSize(); s++) {
			if (!this.emptySlices[s]) {
				short[] pixels = (short[]) stack.getPixels(s);
				for (int y = r.y; y < (r.y + r.height); y++) {
					int offset = y * w;
					for (int x = r.x; x < (r.x + r.width); x++) {
						int i = offset + x;
						if (pixels[i] >= this.min && pixels[i] <= this.max) {
							this.Sx[s] += x;
							this.Sy[s] += y;
							this.Sxx[s] += x * x;
							this.Syy[s] += y * y;
							this.Sxy[s] += y * x;
						}
					}
				}
				this.Myy[s] = this.Sxx[s]
						- (this.Sx[s] * this.Sx[s] / this.cslice[s])
						+ this.cslice[s] / 12;
				// this.cslice[]/12 is for each pixel's own moment
				this.Mxx[s] = this.Syy[s]
						- (this.Sy[s] * this.Sy[s] / this.cslice[s])
						+ this.cslice[s] / 12;
				this.Mxy[s] = this.Sxy[s]
						- (this.Sx[s] * this.Sy[s] / this.cslice[s])
						+ this.cslice[s] / 12;
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
		this.Zmax = new double[this.al];
		this.Zmin = new double[this.al];
		this.ImaxFast = new double[this.al];
		this.IminFast = new double[this.al];
		for (int s = 1; s <= stack.getSize(); s++) {
			if (!this.emptySlices[s]) {
				short[] pixels = (short[]) stack.getPixels(s);
				this.Sx[s] = 0;
				this.Sy[s] = 0;
				this.Sxx[s] = 0;
				this.Syy[s] = 0;
				this.Sxy[s] = 0;
				for (int y = r.y; y < (r.y + r.height); y++) {
					int offset = y * w;
					for (int x = r.x; x < (r.x + r.width); x++) {
						int i = offset + x;
						if (pixels[i] >= this.min && pixels[i] <= this.max) {
							this.Sx[s] += x * Math.cos(this.theta[s]) + y
									* Math.sin(this.theta[s]);
							this.Sy[s] += y * Math.cos(this.theta[s]) - x
									* Math.sin(this.theta[s]);
							this.Sxx[s] += (x * Math.cos(this.theta[s]) + y
									* Math.sin(this.theta[s]))
									* (x * Math.cos(this.theta[s]) + y
											* Math.sin(this.theta[s]));
							this.Syy[s] += (y * Math.cos(this.theta[s]) - x
									* Math.sin(this.theta[s]))
									* (y * Math.cos(this.theta[s]) - x
											* Math.sin(this.theta[s]));
							this.Sxy[s] += (y * Math.cos(theta[s]) - x
									* Math.sin(theta[s]))
									* (x * Math.cos(theta[s]) + y
											* Math.sin(theta[s]));
							this.maxRadMin[s] = Math.max(this.maxRadMin[s],
									Math.abs((x - this.sliceCentroids[0][s])
											* Math.cos(this.theta[s])
											+ (y - this.sliceCentroids[1][s])
											* Math.sin(theta[s])));
							this.maxRadMax[s] = Math.max(this.maxRadMax[s],
									Math.abs((y - this.sliceCentroids[1][s])
											* Math.cos(this.theta[s])
											- (x - sliceCentroids[0][s])
											* Math.sin(theta[s])));
						}
					}
				}
				this.Imax[s] = this.Sxx[s]
						- (this.Sx[s] * this.Sx[s] / this.cslice[s])
						+ this.cslice[s]
						* (Math.pow(Math.cos(this.theta[s]), 2) + Math.pow(Math
								.sin(this.theta[s]), 2)) / 12;
				this.Imin[s] = this.Syy[s]
						- (this.Sy[s] * this.Sy[s] / this.cslice[s])
						+ this.cslice[s]
						* (Math.pow(Math.cos(this.theta[s]), 2) + Math.pow(Math
								.sin(this.theta[s]), 2)) / 12;
				this.Ipm[s] = this.Sxy[s]
						- (this.Sy[s] * this.Sx[s] / this.cslice[s])
						+ this.cslice[s]
						* (Math.pow(Math.cos(this.theta[s]), 2) + Math.pow(Math
								.sin(this.theta[s]), 2)) / 12;
				this.R1[s] = Math.sqrt(this.Imin[s] / this.cslice[s]);
				this.R2[s] = Math.sqrt(this.Imax[s] / this.cslice[s]);
				this.Zmax[s] = this.Imax[s] / this.maxRadMin[s];
				this.Zmin[s] = this.Imin[s] / this.maxRadMax[s];
				this.ImaxFast[s] = (this.Mxx[s] + this.Myy[s])
						/ 2
						+ Math.sqrt(Math.pow(((this.Mxx[s] - this.Myy[s]) / 2),
								2)
								+ this.Mxy[s] * this.Mxy[s]);
				this.IminFast[s] = (this.Mxx[s] + this.Myy[s])
						/ 2
						- Math.sqrt(Math.pow(((this.Mxx[s] - this.Myy[s]) / 2),
								2)
								+ this.Mxy[s] * this.Mxy[s]);
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
				this.ImaxFast[s] = Double.NaN;
				this.IminFast[s] = Double.NaN;
			}
		}
		return;
	}

	/**
	 * Calculate 3D Local Thickness and determine thickness statistics for the
	 * slice
	 * 
	 */
	private void calculateThickness3D(ImagePlus imp) {
		this.maxCortThick3D = new double[this.al];
		this.meanCortThick3D = new double[this.al];
		this.stdevCortThick3D = new double[this.al];
		Thickness th = new Thickness();
		th.baseImp = imp;
		th.vD = this.vD;
		th.vW = this.vW;
		th.vH = this.vH;
		th.w = imp.getWidth();
		th.h = imp.getHeight();
		th.d = imp.getStackSize();

		// convert to binary
		ImagePlus binaryImp = convertToBinary(imp);

		float[][] workArray = th.GeometrytoDistanceMap(binaryImp);
		th.DistanceMaptoDistanceRidge(workArray);
		th.DistanceRidgetoLocalThickness(workArray);
		ImagePlus thickImp = th
				.LocalThicknesstoCleanedUpLocalThickness(workArray);
		for (int s = this.startSlice; s <= this.endSlice; s++) {
			float[] pixels = (float[]) thickImp.getStack().getPixels(s);
			double sumPix = 0;
			double sliceMax = 0;
			double pixCount = 0;
			for (int p = 0; p < pixels.length; p++) {
				if (pixels[p] > 0)
					pixCount++;
				sumPix += pixels[p];
				sliceMax = Math.max(sliceMax, pixels[p]);
			}
			double sliceMean = sumPix / pixCount;
			this.meanCortThick3D[s] = sliceMean;
			if (pixCount > 0)
				this.maxCortThick3D[s] = sliceMax;
			else
				this.maxCortThick3D[s] = Double.NaN;

			double sumSquares = 0;
			for (int p = 0; p < pixels.length; p++) {
				double pixVal = pixels[p];
				if (pixVal > 0) {
					double d = sliceMean - pixVal;
					sumSquares += d * d;
				}
			}
			this.stdevCortThick3D[s] = Math.sqrt(sumSquares / pixCount);
			// IJ.log("Mean thickness for slice "+(n+1)+" is "+this.meanCortThick[n]+" ("+this.stdevCortThick[n]+")");
		}
		return;
	}

	/**
	 * Calculate thickness on individual slices using local thickness
	 * 
	 * @param imp
	 */
	private void calculateThickness2D(ImagePlus imp) {
		this.maxCortThick2D = new double[this.al];
		this.meanCortThick2D = new double[this.al];
		this.stdevCortThick2D = new double[this.al];

		int nThreads = Runtime.getRuntime().availableProcessors();
		SliceThread[] sliceThread = new SliceThread[nThreads];
		for (int thread = 0; thread < nThreads; thread++) {
			sliceThread[thread] = new SliceThread(thread, nThreads, imp,
					this.meanCortThick2D, this.maxCortThick2D,
					this.stdevCortThick2D, this.startSlice, this.endSlice);
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

		final double[] meanThick, maxThick, stdevThick;

		final ImagePlus impT;

		public SliceThread(int thread, int nThreads, ImagePlus imp,
				double[] meanThick, double[] maxThick, double[] stdevThick,
				int startSlice, int endSlice) {
			this.impT = imp;
			this.width = this.impT.getWidth();
			this.height = this.impT.getHeight();
			this.thread = thread;
			this.nThreads = nThreads;
			this.meanThick = meanThick;
			this.maxThick = maxThick;
			this.stdevThick = stdevThick;
			this.startSlice = startSlice;
			this.endSlice = endSlice;
		}

		public void run() {
			for (int s = this.thread + this.startSlice; s <= this.endSlice; s += this.nThreads) {
				ImageProcessor ip = impT.getImageStack().getProcessor(s);
				ImagePlus sliceImp = new ImagePlus(" " + s, ip);
				// binarise
				ImagePlus binaryImp = convertToBinary(sliceImp);
				Calibration cal = impT.getCalibration();
				binaryImp.setCalibration(cal);
				// calculate thickness
				Thickness th = new Thickness();
				th.baseImp = binaryImp;
				th.vD = cal.pixelDepth;
				th.vW = cal.pixelWidth;
				th.vH = cal.pixelHeight;
				th.w = binaryImp.getWidth();
				th.h = binaryImp.getHeight();
				th.d = binaryImp.getStackSize();
				float[][] workArray = th.GeometrytoDistanceMap(binaryImp);
				th.DistanceMaptoDistanceRidge(workArray);
				th.DistanceRidgetoLocalThickness(workArray);
				ImagePlus thickImp = th
						.LocalThicknesstoCleanedUpLocalThickness(workArray);
				// get thickness stats
				float[] pixels = (float[]) thickImp.getImageStack()
						.getPixels(1);
				double sumPix = 0;
				double sliceMax = 0;
				double pixCount = 0;
				for (int p = 0; p < pixels.length; p++) {
					if (pixels[p] > 0)
						pixCount++;
					sumPix += pixels[p];
					sliceMax = Math.max(sliceMax, pixels[p]);
				}
				double sliceMean = sumPix / pixCount;
				this.meanThick[s] = sliceMean;
				if (pixCount > 0)
					this.maxThick[s] = sliceMax;
				else
					this.maxThick[s] = Double.NaN;

				double sumSquares = 0;
				for (int p = 0; p < pixels.length; p++) {
					double pixVal = pixels[p];
					if (pixVal > 0) {
						double d = sliceMean - pixVal;
						sumSquares += d * d;
					}
				}
				this.stdevThick[s] = Math.sqrt(sumSquares / pixCount);
			}
			return;
		}
	}

	private ImagePlus convertToBinary(ImagePlus imp) {
		int w = imp.getWidth();
		int h = imp.getHeight();
		int d = imp.getStackSize();
		ImageStack sourceStack = imp.getImageStack();
		ImageStack binaryStack = new ImageStack(w, h);
		for (int s = 0; s < d; s++) {
			ImageProcessor sliceIp = sourceStack.getProcessor(s + 1);
			ByteProcessor binaryIp = new ByteProcessor(w, h);
			for (int y = 0; y < h; y++) {
				for (int x = 0; x < w; x++) {
					if (sliceIp.get(x, y) >= this.min
							&& sliceIp.get(x, y) <= this.max) {
						binaryIp.set(x, y, 255);
					} else {
						binaryIp.set(x, y, 0);
					}
				}
			}
			binaryStack.addSlice(sourceStack.getSliceLabel(s + 1), binaryIp);
		}
		ImagePlus binaryImp = new ImagePlus("binaryImp", binaryStack);
		binaryImp.setCalibration(imp.getCalibration());
		return binaryImp;
	}

	// TODO fix this, it's a mess.
	/**
	 * Set up threshold based on HU if image is calibrated or user-based
	 * threshold if it is uncalibrated
	 * 
	 */
	private void setHUCalibration(ImagePlus imp) {
		ImageStack stack = imp.getImageStack();
		this.min = 0; // minimum bone value in HU
		this.max = 4000; // maximum bone value in HU
		double[] coeff = this.cal.getCoefficients();
		if (!this.cal.calibrated()
				|| this.cal == null
				|| (this.cal.getCValue(0) == 0 && this.cal.getCoefficients()[1] == 1)) {
			this.isCalibrated = false;
			this.calString = "Image is uncalibrated\nEnter air and bone pixel values";
			ImageStatistics stats = imp.getStatistics();
			if (stats.min < 50 && stats.min >= 0) {
				this.airHU = 0;
			} else if (stats.min > 31000) {
				this.airHU = 31768;
			} else if (stats.min < -800) {
				this.airHU = -1000;
			} else {
				this.airHU = 0;
			}
		} else {
			this.isCalibrated = true;
			this.calString = "Image is calibrated\nEnter HU below:";
			this.airHU = -1000;
		}
		this.min = this.airHU + 1000;
		this.max = this.airHU + 5000;
		if (this.cal.calibrated()) {
			// convert HU limits to pixel values
			IJ.log("Image is calibrated, using " + this.min + " and "
					+ this.max + " HU as bone cutoffs");
			this.min = (short) Math.round(this.cal.getRawValue(this.min));
			this.max = (short) Math.round(this.cal.getRawValue(this.max));
			IJ.log("Vox Width: " + vW + "; Vox Height: " + vH + " " + units);
			IJ.log("Calibration coefficients:" + coeff[0] + "," + coeff[1]);
			IJ.log("this.min = " + this.min + ", this.max = " + this.max);
		} else {
			IJ.log("Image is not calibrated, using user-determined threshold");
			IJ.run("Threshold...");
			new WaitForUserDialog(
					"This image is not density calibrated.\nSet the threshold, then click OK.")
					.show();
			this.min = (short) stack.getProcessor(imp.getCurrentSlice())
					.getMinThreshold();
			this.max = (short) stack.getProcessor(imp.getCurrentSlice())
					.getMaxThreshold();
		}
		return;
	}

	private boolean showDialog(ImagePlus imp) {
		GenericDialog gd = new GenericDialog("Options");

		gd.addCheckbox("2D Thickness", true);
		gd.addCheckbox("3D Thickness", false);
		gd.addCheckbox("Draw Axes", true);
		gd.addCheckbox("Draw Centroids", true);
		gd.addCheckbox("Annotated Copy", true);
		gd.addCheckbox("Process Stack", false);
		// guess bone from image title
		BoneList bl = new BoneList();
		this.boneID = bl.guessBone(imp);
		String[] bones = bl.getBoneList();
		gd.addChoice("Bone: ", bones, bones[this.boneID]);
		String[] analyses = { "Weighted", "Unweighted", "Both" };
		gd.addChoice("Calculate: ", analyses, analyses[1]);
		gd.addMessage("Set the threshold");
		gd.addNumericField("Air:", this.airHU, 0);
		gd.addNumericField("Bone Min:", this.min, 0);
		gd.addNumericField("Bone Max:", this.max, 0);
		gd.showDialog();
		this.doThickness2D = gd.getNextBoolean();
		this.doThickness3D = gd.getNextBoolean();
		this.doAxes = gd.getNextBoolean();
		this.doCentroids = gd.getNextBoolean();
		this.doCopy = gd.getNextBoolean();
		this.doStack = gd.getNextBoolean();
		if (this.doStack) {
			this.startSlice = 1;
			this.endSlice = imp.getImageStackSize();
		} else {
			this.startSlice = imp.getCurrentSlice();
			this.endSlice = imp.getCurrentSlice();
		}

		String bone = gd.getNextChoice();
		this.boneID = bl.guessBone(bone);
		this.analyse = gd.getNextChoice();
		this.airHU = gd.getNextNumber();
		this.min = gd.getNextNumber();
		this.max = gd.getNextNumber();
		if (gd.wasCanceled()) {
			return false;
		} else {
			return true;
		}
	}

	private void showSliceResults(ImagePlus imp) {
		ResultsTable rt = ResultsTable.getResultsTable();
		rt.reset();

		double unit4 = Math.pow(vW, 4);
		double unit3 = Math.pow(vW, 3);
		String title = imp.getTitle();
		for (int s = this.startSlice; s <= this.endSlice; s++) {
			rt.incrementCounter();
			rt.addLabel(title);
			rt.addValue("Bone Code", this.boneID);
			rt.addValue("Slice", s);
			rt.addValue("CA (" + units + "^2)", this.cortArea[s]);
			rt.addValue("X cent. (" + units + ")", this.sliceCentroids[0][s]);
			rt.addValue("Y cent. (" + units + ")", this.sliceCentroids[1][s]);
			rt.addValue("Theta (rad)", this.theta[s]);
			rt.addValue("R1 (" + units + ")", this.R1[s]);
			rt.addValue("R2 (" + units + ")", this.R2[s]);
			rt.addValue("Imin (" + units + "^4)", this.Imin[s] * unit4);
			rt.addValue("IminFast (" + units + "^4)", this.IminFast[s] * unit4);
			rt.addValue("Imax (" + units + "^4)", this.Imax[s] * unit4);
			rt.addValue("ImaxFast (" + units + "^4)", this.ImaxFast[s] * unit4);
			rt.addValue("Ipm (" + units + "^4)", this.Ipm[s] * unit4);
			rt.addValue("Zmax (" + units + "^3)", this.Zmax[s] * unit3);
			rt.addValue("Zmin (" + units + "^3)", this.Zmin[s] * unit3);
			rt.addValue("Feret Min", this.feretMin[s]);
			rt.addValue("Feret Max", this.feretMax[s]);
			rt.addValue("Feret Angle", this.feretAngle[s]);
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
	}

	private void roiMeasurements(ImagePlus imp) {
		double[] feretValues = new double[3];
		this.feretAngle = new double[this.al];
		this.feretMax = new double[this.al];
		this.feretMin = new double[this.al];
		imp.setActivated();
		Roi r;
		IJ.setThreshold(this.min, this.max);
		// for the required slices...
		for (int s = this.startSlice; s <= this.endSlice; s++) {
			IJ.setSlice(s);
			IJ.doWand(0, (int) Math.round(this.sliceCentroids[1][s] / this.vH));
			r = imp.getRoi();
			if (this.emptySlices[s]) {
				this.feretMin[s] = Double.NaN;
				this.feretAngle[s] = Double.NaN;
				this.feretMax[s] = Double.NaN;
			} else if (r != null) {
				feretValues = r.getFeretValues();
				this.feretMin[s] = feretValues[2];
				this.feretAngle[s] = feretValues[1];
				this.feretMax[s] = feretValues[0];
			} else {
				this.feretMin[s] = Double.NaN;
				this.feretAngle[s] = Double.NaN;
				this.feretMax[s] = Double.NaN;
			}
			r = null;
			feretValues = null;
		}
	}
}