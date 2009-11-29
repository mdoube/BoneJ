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
import ij.plugin.PlugIn;
import ij.measure.Calibration;
import ij.measure.ResultsTable;
import ij.gui.*;

import java.awt.Color;
import java.awt.Rectangle;

import org.doube.bonej.BoneList;
import org.doube.bonej.ImageCheck;
import org.doube.bonej.Thickness;
import org.doube.bonej.ThresholdMinConn;

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
	private double vW, vH, min, max;
	/** Linear unit of measure */
	private String units;
	/** Unused option to do density-weighted calculations */
	// private String analyse;
	/** Message to inform the user what to do with their HU-calibrated image */
	private String calString;
	/** Hounsfield unit value for air is -1000 */
	private static final double airHU = -1000;
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
	/** Is this a Hounsfield Unit calibrated image */
	private boolean isHUCalibrated;
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
		if (!ImageCheck.checkIJVersion())
			return;
		ImagePlus imp = IJ.getImage();
		if (null == imp) {
			IJ.noImage();
			return;
		}

		this.cal = imp.getCalibration();
		this.vW = this.cal.pixelWidth;
		this.vH = this.cal.pixelHeight;
		this.units = this.cal.getUnits();
		this.al = imp.getStackSize() + 1;

		setDefaultThreshold(imp);

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

	protected double calculateCentroids(ImagePlus imp) {
		ImageStack stack = imp.getImageStack();
		Rectangle r = stack.getRoi();
		// 2D centroids
		this.sliceCentroids = new double[2][this.al];
		// pixel counters
		double cstack = 0;
		this.emptySlices = new boolean[this.al];
		this.cslice = new double[this.al];
		this.cortArea = new double[this.al];
		final double pixelArea = this.vW * this.vH;
		final int roiXEnd = r.x + r.width;
		final int roiYEnd = r.y + r.height;
		for (int s = this.startSlice; s <= this.endSlice; s++) {
			IJ.showStatus("Calculating centroids...");
			IJ.showProgress(s - this.startSlice, this.endSlice);
			double sumX = 0;
			double sumY = 0;
			this.cslice[s] = 0;
			this.cortArea[s] = 0;
			ImageProcessor ip = stack.getProcessor(s);
			for (int y = r.y; y < roiYEnd; y++) {
				for (int x = r.x; x < roiXEnd; x++) {
					final double pixel = (double) ip.get(x, y);
					if (pixel >= this.min && pixel <= this.max) {
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

	/**
	 * Calculate second moments of area, length and angle of principal axes
	 * 
	 * @param imp
	 */
	private void calculateMoments(ImagePlus imp) {
		final ImageStack stack = imp.getImageStack();
		final Rectangle r = stack.getRoi();
		final int d = stack.getSize();
		final double pMin = this.min;
		final double pMax = this.max;
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
		for (int s = 1; s <= d; s++) {
			IJ.showStatus("Calculating Ix and Iy...");
			IJ.showProgress(s, d);
			this.Sx[s] = 0;
			this.Sy[s] = 0;
			this.Sxx[s] = 0;
			this.Syy[s] = 0;
			this.Sxy[s] = 0;
			final int roiXEnd = r.x + r.width;
			final int roiYEnd = r.y + r.height;
			if (!this.emptySlices[s]) {
				ImageProcessor ip = stack.getProcessor(s);
				for (int y = r.y; y < roiYEnd; y++) {
					for (int x = r.x; x < roiXEnd; x++) {
						final double pixel = (double) ip.get(x, y);
						if (pixel >= pMin && pixel <= pMax) {
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
		for (int s = 1; s <= d; s++) {
			IJ.showStatus("Calculating Imin and Imax...");
			IJ.showProgress(s, d);
			if (!this.emptySlices[s]) {
				ImageProcessor ip = stack.getProcessor(s);
				this.Sx[s] = 0;
				this.Sy[s] = 0;
				this.Sxx[s] = 0;
				this.Syy[s] = 0;
				this.Sxy[s] = 0;
				final double cosTheta = Math.cos(this.theta[s]);
				final double sinTheta = Math.sin(this.theta[s]);
				final int roiYEnd = r.y + r.height;
				final int roiXEnd = r.x + r.width;
				final double xC = this.sliceCentroids[0][s];
				final double yC = this.sliceCentroids[1][s];
				final double cS = this.cslice[s];
				for (int y = r.y; y < roiYEnd; y++) {
					for (int x = r.x; x < roiXEnd; x++) {
						final double pixel = (double) ip.get(x, y);
						if (pixel >= pMin && pixel <= pMax) {
							final double xCosTheta = x * cosTheta;
							final double yCosTheta = y * cosTheta;
							final double xSinTheta = x * sinTheta;
							final double ySinTheta = y * sinTheta;
							this.Sx[s] += xCosTheta + ySinTheta;
							this.Sy[s] += yCosTheta - xSinTheta;
							this.Sxx[s] += (xCosTheta + ySinTheta)
									* (xCosTheta + ySinTheta);
							this.Syy[s] += (yCosTheta - xSinTheta)
									* (yCosTheta - xSinTheta);
							this.Sxy[s] += (yCosTheta - xSinTheta)
									* (xCosTheta + ySinTheta);
							this.maxRadMin[s] = Math.max(this.maxRadMin[s],
									Math.abs((x - xC) * cosTheta + (y - yC)
											* sinTheta));
							this.maxRadMax[s] = Math.max(this.maxRadMax[s],
									Math.abs((y - yC) * cosTheta - (x - xC)
											* sinTheta));
						}
					}
				}
				final double pixelMoments = cS
						* (cosTheta * cosTheta + sinTheta * sinTheta) / 12;
				this.Imax[s] = this.Sxx[s] - (this.Sx[s] * this.Sx[s] / cS)
						+ pixelMoments;
				this.Imin[s] = this.Syy[s] - (this.Sy[s] * this.Sy[s] / cS)
						+ pixelMoments;
				this.Ipm[s] = this.Sxy[s] - (this.Sy[s] * this.Sx[s] / cS)
						+ pixelMoments;
				this.R1[s] = Math.sqrt(this.Imin[s] / cS);
				this.R2[s] = Math.sqrt(this.Imax[s] / cS);
				this.Zmax[s] = this.Imax[s] / this.maxRadMin[s];
				this.Zmin[s] = this.Imin[s] / this.maxRadMax[s];
				final double a = (this.Mxx[s] + this.Myy[s]) / 2;
				final double b = Math.sqrt((this.Mxx[s] - this.Myy[s])
						* (this.Mxx[s] - this.Myy[s]) / 4 + this.Mxy[s]
						* this.Mxy[s]);
				this.ImaxFast[s] = a + b;
				this.IminFast[s] = a - b;
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
		Rectangle r = imp.getProcessor().getRoi();
		Thickness th = new Thickness();

		// convert to binary
		ImagePlus binaryImp = convertToBinary(imp);

		ImagePlus thickImp = th.getLocalThickness(binaryImp, false);
		for (int s = this.startSlice; s <= this.endSlice; s++) {
			ImageProcessor ip = thickImp.getStack().getProcessor(s);
			double sumPix = 0;
			double sliceMax = 0;
			double pixCount = 0;
			final int roiXEnd = r.x + r.width;
			final int roiYEnd = r.y + r.height;
			for (int y = r.y; y < roiYEnd; y++) {
				for (int x = r.x; x < roiXEnd; x++) {
					final double pixel = (double) ip.get(x, y);
					if (pixel > 0)
						pixCount++;
					sumPix += pixel;
					sliceMax = Math.max(sliceMax, pixel);
				}
			}
			final double sliceMean = sumPix / pixCount;
			this.meanCortThick3D[s] = sliceMean;
			if (pixCount > 0)
				this.maxCortThick3D[s] = sliceMax;
			else
				this.maxCortThick3D[s] = Double.NaN;

			double sumSquares = 0;
			for (int y = r.y; y < roiYEnd; y++) {
				for (int x = r.x; x < roiXEnd; x++) {
					final double d = sliceMean - (double) ip.get(x, y);
					sumSquares += d * d;
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
				Rectangle r = ip.getRoi();
				// binarise
				ImagePlus binaryImp = convertToBinary(sliceImp);
				Calibration cal = impT.getCalibration();
				binaryImp.setCalibration(cal);
				// calculate thickness
				Thickness th = new Thickness();
				ImagePlus thickImp = th.getLocalThickness(binaryImp, false);
				ImageProcessor thickIp = thickImp.getImageStack().getProcessor(
						s);
				double sumPix = 0;
				double sliceMax = 0;
				double pixCount = 0;
				final double roiXEnd = r.x + r.width;
				final double roiYEnd = r.y + r.height;
				for (int y = r.y; y < roiYEnd; y++) {
					for (int x = r.x; x < roiXEnd; x++) {
						final double pixel = (double) thickIp.get(x, y);
						if (pixel > 0)
							pixCount++;
						sumPix += pixel;
						sliceMax = Math.max(sliceMax, pixel);
					}
				}
				final double sliceMean = sumPix / pixCount;
				this.meanThick[s] = sliceMean;
				if (pixCount > 0)
					this.maxThick[s] = sliceMax;
				else
					this.maxThick[s] = Double.NaN;

				double sumSquares = 0;
				for (int y = r.y; y < roiYEnd; y++) {
					for (int x = r.x; x < roiXEnd; x++) {
						final double pixel = (double) thickIp.get(x, y);
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

	private ImagePlus convertToBinary(ImagePlus imp) {
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
					if (sliceIp.get(x, y) >= this.min
							&& sliceIp.get(x, y) <= this.max) {
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

	/**
	 * Set up default thresholds and report them to the user as HU values if the
	 * image has HU calibration or plain values if not. Sets min and max pixel
	 * values as raw pixel values if image is uncalibrated and as HU values if
	 * if it is.
	 * 
	 */
	private void setDefaultThreshold(ImagePlus imp) {
		this.cal = imp.getCalibration();

		double[] coeff = this.cal.getCoefficients();
		if (!this.cal.calibrated() || this.cal == null
				|| (this.cal.getCValue(0) == 0 && coeff[1] == 1)) {
			this.isHUCalibrated = false;
			this.calString = "Image is uncalibrated\nEnter bone threshold values";
			// set some sensible thresholding defaults
			ThresholdMinConn tmc = new ThresholdMinConn();
			int[] histogram = tmc.getStackHistogram(imp);
			final int histoLength = histogram.length;
			int histoMax = histoLength - 1;
			for (int i = histoLength - 1; i >= 0; i--) {
				if (histogram[i] > 0) {
					histoMax = i;
					break;
				}
			}
			this.min = imp.getProcessor().getAutoThreshold(histogram);
			this.max = histoMax;
			if (this.cal.isSigned16Bit()) {
				this.min += Short.MIN_VALUE;
				this.max += Short.MIN_VALUE;
			}

		} else {
			this.isHUCalibrated = true;
			this.calString = "Image has Hounsfield calibration.\nEnter bone HU below:";
			// default bone thresholds are 0- and 4000 HU
			this.min = Slice_Geometry.airHU + 1000;
			this.max = Slice_Geometry.airHU + 5000;
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
		String[] bones = BoneList.getBoneList();
		gd.addChoice("Bone: ", bones, bones[this.boneID]);
		// String[] analyses = { "Weighted", "Unweighted", "Both" };
		// gd.addChoice("Calculate: ", analyses, analyses[1]);
		gd.addMessage(this.calString);
		if (this.isHUCalibrated)
			gd.addMessage("HU for air is " + Slice_Geometry.airHU);
		gd.addNumericField("Bone Min:", this.min, 0);
		gd.addNumericField("Bone Max:", this.max, 0);
		gd
				.addMessage("Only pixels <= bone min\n"
						+ "and >= bone max are used.");
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
		// this.analyse = gd.getNextChoice();
		this.min = gd.getNextNumber();
		this.max = gd.getNextNumber();
		if (this.isHUCalibrated) {
			this.min = this.cal.getRawValue(this.min);
			this.max = this.cal.getRawValue(this.max);
		}
		if (gd.wasCanceled()) {
			return false;
		} else {
			return true;
		}
	}

	private void showSliceResults(ImagePlus imp) {
		ResultsTable rt = ResultsTable.getResultsTable();
		rt.reset();

		final double unit4 = Math.pow(vW, 4);
		final double unit3 = Math.pow(vW, 3);
		String title = imp.getTitle();
		for (int s = this.startSlice; s <= this.endSlice; s++) {
			rt.incrementCounter();
			rt.addLabel(title);
			rt.addValue("Bone Code", this.boneID);
			rt.addValue("Slice", s);
			rt.addValue("CSA (" + units + "^2)", this.cortArea[s]);
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
			rt.addValue("Feret Min (" + units + ")", this.feretMin[s]);
			rt.addValue("Feret Max (" + units + ")", this.feretMax[s]);
			rt.addValue("Feret Angle (rad)", this.feretAngle[s]);
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
		Roi initialRoi = imp.getRoi();
		double[] feretValues = new double[3];
		this.feretAngle = new double[this.al];
		this.feretMax = new double[this.al];
		this.feretMin = new double[this.al];
		int initialSlice = imp.getCurrentSlice();
		// for the required slices...
		for (int s = this.startSlice; s <= this.endSlice; s++) {
			ImageProcessor ip = imp.getImageStack().getProcessor(s);
			Wand w = new Wand(ip);
			w.autoOutline(0, (int) Math.round(this.sliceCentroids[1][s]
					/ this.vH), this.min, this.max, Wand.EIGHT_CONNECTED);
			if (this.emptySlices[s] || w.npoints == 0) {
				this.feretMin[s] = Double.NaN;
				this.feretAngle[s] = Double.NaN;
				this.feretMax[s] = Double.NaN;
			} else {
				int type = Wand.allPoints() ? Roi.FREEROI : Roi.TRACED_ROI;
				Roi roi = new PolygonRoi(w.xpoints, w.ypoints, w.npoints, type);
				feretValues = roi.getFeretValues();
				this.feretMin[s] = feretValues[2] * this.vW;
				this.feretAngle[s] = feretValues[1] * Math.PI / 180;
				this.feretMax[s] = feretValues[0] * this.vW;
			}
			feretValues = null;
		}
		IJ.setSlice(initialSlice);
		imp.setRoi(initialRoi);
		return;
	}
}
