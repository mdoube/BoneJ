package org.doube.bonej;

import java.util.Arrays;

import org.doube.util.DeleteSliceRange;
import org.doube.util.ImageCheck;
import org.doube.util.ThresholdGuesser;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.gui.Plot;
import ij.measure.ResultsTable;
import ij.plugin.PlugIn;

/**
 * <p>
 * Make an informed guess at which slices of an image stack contain the
 * proximal and distal "ends" of a femoral shaft.
 * </p>
 * 
 * @author Nick Powell
 *
 */

public class ShaftGuesser implements PlugIn {
	
	private boolean proximalLow, doGraph;
	
	private int al, startSlice, endSlice, ss, iss;
	/** Smooth over this number of slices. Initially set to 1/50th of stack size */
	private int smoothOver;
	/** Calculate the gradient over this number of slices */
	private int gradientOver;
	/** First and last slice numbers of the femoral shaft */
	private int[] shaftPosition;
	
	/** Median values */
	private double medianMeanCort, mFeretMax, mFeretMin, mPerimeter;
	/** List of eccentricities, based on Feret's min and max */
	private double[] eccentricity;
	/** List of maximum diameters */
	private double[] feretMax;
	/** List of minimum diameters */
	private double[] feretMin;
	/** List of mean cortical thicknesses */
	private double[] meanCortThick2D;
	/** List of perimeter lengths */
	private double[] perimeter;
	/** List of slice numbers */
	private double[] slices;
	/** First and last slices of shaft */
	private double[] shaft = new double[2];
	/** List of gradients at each slice, based on values +/- from that slice */
	private double[] gFeretMax, gFeretMin, gPerimeter, gEccentricity, gMeanCort;
	/** List of smoothed values */
	private double[] sFeretMax, sFeretMin, sPerimeter, sEccentricity, sMeanCort;
	
	public void run(String arg) {
		
		/* ImageJ version checking */
		if (!ImageCheck.checkEnvironment())
			return;
		ImagePlus imp = IJ.getImage();
		if (null == imp) {
			IJ.noImage();
			return;
		}
		
		this.al = imp.getStackSize() + 1;
		this.iss = imp.getImageStackSize();
		
		this.slices = new double[this.al];
		for (int s = this.startSlice; s <= this.iss; s++) {
			slices[s] = (double) s;
		}
		
		double[] thresholds = ThresholdGuesser.setDefaultThreshold(imp);
		double min = thresholds[0];
		double max = thresholds[1];
		
		GenericDialog gd = new GenericDialog("Shaft Guesser Options");
		gd.addCheckbox("Proximal end of bone has lower slice number", true);
		gd.addNumericField("Bone start slice", 1, 0);
		gd.addNumericField("Bone end slice", this.iss, 0);
		gd.addMessage("Choose how to estimate the shaft parameters: ");
		gd.addNumericField("Smooth over +/-", Math.round(this.al / 50), 0, 3, "slices");
		gd.addNumericField("Calculate gradient over +/-", Math.round(this.al / 50), 0, 3, "slices");
		gd.addCheckbox("Graph output", true);
		gd.showDialog();
		
		this.proximalLow = gd.getNextBoolean();
		this.startSlice = (int) gd.getNextNumber();
		this.endSlice = (int) gd.getNextNumber();
		this.smoothOver = (int) gd.getNextNumber();
		this.gradientOver = (int) gd.getNextNumber();
		this.doGraph = gd.getNextBoolean();
		if (gd.wasCanceled())
			return;
		
		/* Initial setup */
		if(!proximalLow) {
			IJ.run("Flip Z");
		}
		
		/* Crop ends */
		if(this.startSlice > 1) {
			DeleteSliceRange dsr1 = new DeleteSliceRange();
			dsr1.deleteSliceRange(imp.getImageStack(), 1, this.startSlice - 1);
			this.iss = imp.getImageStackSize();
			this.endSlice = this.endSlice - (this.startSlice - 1);
			this.startSlice = 1;
		}
		if(this.endSlice < this.iss) {
			DeleteSliceRange dsr2 = new DeleteSliceRange();
			dsr2.deleteSliceRange(imp.getImageStack(), this.endSlice + 1, this.iss);
			this.endSlice = this.iss;
		}
		
		/* Get values from SliceGeometry */
		SliceGeometry sg = new SliceGeometry();
		sg.setParameters(imp);
		sg.calculateCentroids(imp, min, max);
		sg.roiMeasurements(imp, min, max);
		sg.calculateMoments(imp, min, max);
		sg.calculateThickness2D(imp, min, max);
		this.feretMax = sg.getFeretMax();
		this.feretMin = sg.getFeretMin();
		this.perimeter = sg.getPerimeter();
		this.meanCortThick2D = sg.getMeanCorticalThickness2D();
		
		/* Calculate the eccentricities of each slice */
		this.eccentricity = eccentricity(feretMin, feretMax);
		
		/* Smooth out noise */
		this.sEccentricity = smooth(eccentricity, smoothOver);
		this.sFeretMax = smooth(feretMax, smoothOver);
		this.sFeretMin = smooth(feretMin, smoothOver);
		this.sPerimeter = smooth(perimeter, smoothOver);
		this.sMeanCort = smooth(meanCortThick2D, smoothOver);
		
		/* Calculate 'gradients' (changes over a number of slices), preferably after smoothing */
		this.gEccentricity = gradient(sEccentricity, gradientOver);
		this.gFeretMax = gradient(sFeretMax, gradientOver);
		this.gFeretMin = gradient(sFeretMin, gradientOver);
		this.gPerimeter = gradient(sPerimeter, gradientOver);
		this.gMeanCort = gradient(sMeanCort, gradientOver);
		
		/* Calculate median values */
		this.mFeretMax = median(sFeretMax);
		this.mFeretMin = median(sFeretMin);
		this.mPerimeter = median(sPerimeter);
		this.medianMeanCort = median(sMeanCort);
		
		/* Possible outer limits of shaft */
		/* Perimeter measures relative size;
		 * Cortical thickness is an actual bone property;
		 * Eccentricity measures shape.
		 * 
		 * 
		 * 
		 * Poss: don't count slices with NaN...?
		 */
		
		/* Run 1: possible outer limits of shaft */
		this.shaftPosition = shaftLimiter(sPerimeter, mPerimeter, false);
		
		/* Copy image to a new stack (not yet) */
		
		ResultsTable rt = ResultsTable.getResultsTable();
		rt.incrementCounter();
		rt.addValue("Shaft start slice", shaftPosition[0]);
		rt.addValue("Shaft end slice", shaftPosition[1]);
		rt.addValue("Median Feret Max", mFeretMax);
		rt.addValue("Median Feret Min", mFeretMin);
		rt.addValue("Median Perimeter", mPerimeter);
		rt.addValue("Median Mean Cortical Thickness (2D)", medianMeanCort);
		rt.show("Results");
		
		if(doGraph) {
			
			Plot sPerimPlot = new Plot("Smoothed perimeter", "Slice", "Perimeter", this.slices, this.sPerimeter);
			sPerimPlot.show();
			
			Plot sMeanCortPlot = new Plot("Smoothed Mean Cort Thick (2D)", "Slice", "Mean Cort Thickness", this.slices, this.sMeanCort);
			sMeanCortPlot.show();
		}
		
	}
	
	/**
	 * Estimates at which slices the shaft begins and ends, based on numerical limits.
	 * Cycles through slices from the centre, first down, then up, until 
	 * a condition is met (shaftLimit).
	 * Main weakness: assumes the central slice is part of the shaft.
	 * 
	 * @param boneStack
	 * @param shaftLimit
	 * @param isMin Specify whether the shaftLimit value is a minimum (true) or a maximum (false)
	 * @return
	 */
	private int[] shaftLimiter(double[] boneStack, double shaftLimit, boolean isMin) {
		
		int[] shaftEnds = new int[2];
		// Find the central slice. (For stacks with odd number of slices, should round up ImageStackSize/2.)
		int centralSlice = Math.round(boneStack.length / 2);
		
		if(isMin) {
			// Cycle through slices from central slice, first down; then up.
			for(int i = centralSlice; i >= this.startSlice; i--) {
				if(Math.abs(boneStack[i]) < shaftLimit) {
					shaftEnds[0] = i + 1;
					break;
				}
			}
			for(int i = centralSlice; i < this.iss; i++) {
				if(Math.abs(boneStack[i]) < shaftLimit) {
					shaftEnds[1] = i - 1;
					break;
				}
			}
		}
		else {
			// Cycle through slices from central slice, first down; then up.
			for(int i = centralSlice; i >= this.startSlice; i--) {
				if(Math.abs(boneStack[i]) > shaftLimit) {
					shaftEnds[0] = i + 1;
					break;
				}
			}
			for(int i = centralSlice; i < this.iss; i++) {
				if(Math.abs(boneStack[i]) > shaftLimit) {
					shaftEnds[1] = i - 1;
					break;
				}
			}
		}
		
		return shaftEnds;
	}
	
	/**
	 * Eccentricity, 0 < e < 1 (0 is a circle), for a pair of doubles[]
	 * Based on Feret's min and max.
	 * 
	 * See <a href="http://mathworld.wolfram.com/Eccentricity.html">http://mathworld.wolfram.com/Eccentricity.html</a>
	 *
	 * @param feretMin[i] is less than feretMax[i]
	 */
	private double[] eccentricity(double feretMin[], double feretMax[]) {
		
		double[] e = new double[feretMax.length];
		for(int i = 0; i < feretMax.length; i++) {
			e[i] = Math.sqrt(1 - ((feretMin[i] * feretMin[i]) / (feretMax[i] * feretMax[i])));
		}
		return e;
	}
	
	/**
	 * Find the gradient of a double[] by comparing x[i-gradSlices] with x[i+gradSlices].
	 * Close to the ends of a stack, only uses available slices.
	 * 
	 * @param x
	 * @param a number of slices over which to find the gradient: resolution
	 * @return gradient[]
	 */
	private double[] gradient(double[] x, int a) {
		
		double[] gradient = new double[x.length];
		
		for(int s = 0; s < x.length; s++) {
			
			// Don't sum slices < 0
			if(s - a < 0) {
				gradient[s] = (x[s + a] - x[s]) / a;
			}
			
			// Don't sum slices > x.length - 1
			else if(s + a > x.length - 1) {
				gradient[s] = (x[s] - x[s - a]) / a;
			}
			else {
				gradient[s] = (x[s + a] - x[s - a]) / (2 * a);
			}
		}

		return gradient;
	}
	
	/**
	 * Removes noise from a double[] by averaging over the surrounding slices.
	 * Near the ends, only uses available slices.
	 * 
	 * @param x
	 * @param a number of slices +/- over which to average
	 * @return gradient[]
	 */
	private double[] smooth(double[] x, int a) {
		
		double[] average = new double[x.length];
		
		for(int s = 0; s < x.length; s++) {
			
			double sum = 0;
			
			// Don't sum slices < 0
			if(s - a < 0) {
				for(int i = 0; i <= a; i++) {
					sum += x[s+i];
				}
				average[s] = sum / (a/2);
			}
			
			// Don't sum slices > x.length - 1
			else if(s + a > x.length - 1) {
				for(int i = 0; i <= a; i++) {
					sum += x[s-i];
				}
				average[s] = sum / (a/2);
			}
			else {
				for(int i = 0; i <= a; i++) {
					sum += x[s-i];
				}
				for(int i = 0; i <= a; i++) {
					sum += x[s+i];
				}
				average[s] = sum / a;
			}
		}
		
		return average;
	}
	
	/**
	 * Find the median of a double[]
	 */
	public static double median(double[] a) {
		
		double[] b = a.clone();
		int mid = b.length / 2;
		Arrays.sort(b);
		
		if(a.length%2 == 1) {
			return b[mid];
		}
		else {
			return (b[mid - 1] + b[mid]) / 2;
		}
	}
}
