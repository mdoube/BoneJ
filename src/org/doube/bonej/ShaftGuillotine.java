package org.doube.bonej;

import org.doube.bonej.SliceGeometry;
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
 * Uses org.doube.bonej.SliceGeometry to calculate maximum Feret diameter (MFD)/slice.
 * Isolates the shaft numerically, using various calculations on MFD.
 * 
 * Intended for use on femurs.
 * 
 * NB. (Probably) best to align bone first, before carving up.
 * 
 * </p>
 * 
 * @author Nick Powell
 *
 */

public class ShaftGuillotine implements PlugIn {
	
	private boolean proximalLow, doGraph, doFeret, cropEnds, do3DMarkers, doManuallyThreshold;
	private boolean[] emptySlices;
	
	private double shaftSD, boneSD, gradSD, meanMFD, shaftCutOff, min, max;
	private double[] sliceCol, feretMax, feretMin, shaft, boneGradientNoNaNs, feretMaxOnMin, perimeter;
	private double[][] boneGradient, perimeterGradient;
	
	private int al, startSlice, iss, centralSlice, gradSlices;
	private int[] shaftPosition;			// Contains first and last shaft slices
	
	private String shaftOOne = "0.25 x standard deviation of gradient";
	private String shaftOTwo = "Mean(max FD) + 0.5 standard deviations of max FD";
	private String[] shaftViaChoices = {shaftOOne, shaftOTwo};

	public void run(String arg) {
		
		/* ImageJ version checking */
		if (!ImageCheck.checkEnvironment())
			return;
		ImagePlus imp = IJ.getImage();
		if (null == imp) {
			IJ.noImage();
			return;
		}
		
		GenericDialog gd = new GenericDialog("Options");
		gd.addCheckbox("Proximal end of bone has lower slice number", true);
		gd.addCheckbox("Graph output", true);
//		gd.addCheckbox("Crop ends", true);
//		gd.addCheckbox("Max Feret diameter", true);
		gd.addChoice("Shaft selection via", shaftViaChoices, shaftOOne);
		gd.addNumericField("Calculate gradient over # slices: ", 1, 0);
//		gd.addCheckbox("Show 3D Markers at end of shaft", true);
		gd.addCheckbox("Manually threshold", false);
		gd.addNumericField("Manually threshold: min", 30, 0);
		gd.addNumericField("Manually threshold: max", 255, 0);
		gd.showDialog();
		
		this.proximalLow = gd.getNextBoolean();
		this.doGraph = gd.getNextBoolean();
//		this.cropEnds = gd.getNextBoolean();
//		this.doFeret = gd.getNextBoolean();
		String shaftVia = gd.getNextChoice();
		this.gradSlices = (int) gd.getNextNumber();
//		this.do3DMarkers = gd.getNextBoolean();
		this.doManuallyThreshold = gd.getNextBoolean();
		this.min = gd.getNextNumber();
		this.max = gd.getNextNumber();
		
		if(!proximalLow) {
			IJ.run("Flip Z");
		}
		
		if(!doManuallyThreshold) {
			/** Copied from SliceGeometry */
			double[] thresholds = ThresholdGuesser.setDefaultThreshold(imp);
			this.min = thresholds[0];
			this.max = thresholds[1];
		}
		
		this.al = imp.getStackSize() + 1;
		this.startSlice = 1;
		this.iss = imp.getImageStackSize();
		
//		if(doFeret) {
//			
//		}
		
		SliceGeometry sg = new SliceGeometry();
		sg.setParameters(imp);
		sg.calculateCentroids(imp, min, max);
		sg.roiMeasurements(imp, min, max);
		this.feretMax = sg.getFeretMax();
		this.feretMin = sg.getFeretMin();
		this.perimeter = sg.getPerimeter();
		this.emptySlices = sg.getEmptySlices();
		
		this.sliceCol = new double[this.al];
		for (int s = this.startSlice; s <= this.iss; s++) {
			sliceCol[s] = (double) s;
		}
		
		// Ellipsoidness
		this.feretMaxOnMin = new double[feretMax.length];
		for(int i = 0; i < feretMax.length; i++) {
			feretMaxOnMin[i] = feretMax[i] / feretMin[i];
		}
		
		// Gradient
		this.boneGradient = gradient(sliceCol, feretMax, gradSlices);
		this.perimeterGradient = gradient(sliceCol, perimeter, gradSlices);
		
		this.shaftPosition = new int[2];
		
		// Arbitrary: find the middle third of the entire stack.
		this.shaftPosition[0] = Math.round(this.iss / 3);	// NB. Math.round() rounds down
		this.shaftPosition[1] = Math.round(this.iss * 2/3);
		this.shaft = new double[shaftPosition[1] - shaftPosition[0]];
		int j = 0;
		for(int s = this.shaftPosition[0]; s < this.shaftPosition[1]; s++) {
			shaft[j] = this.feretMax[s];
			j++;
		}
		
		// Get rid of NaNs (required for standard deviations)
		this.shaft = removeNaNs(shaft);
		this.feretMax = removeNaNs(feretMax);
		this.boneGradientNoNaNs = removeNaNs(boneGradient[1]);
		
		// Standard deviations
		this.shaftSD = Math.sqrt(variance(shaft));
		this.boneSD = Math.sqrt(variance(feretMax));
		this.gradSD = Math.sqrt(variance(boneGradientNoNaNs));
		
		// Gradient-based: identifies changes in gradient, rather than in raw feretMax
		if(shaftVia == shaftOOne) {
			
			this.shaftCutOff = 0.25 * this.gradSD;
			this.shaftPosition = shaftGuesser(boneGradientNoNaNs, shaftCutOff, gradSlices);
		}
		
		/* Uses Mean(max Feret diameter) and Standard Deviation(max FD).
		 * Also relies on arbitrary 'mid third', using shaftSD.
		 * Measurement is Mean + 1/2 SD */
		if(shaftVia == shaftOTwo) {
			
			this.shaftCutOff = mean(feretMax) + (0.5 * this.shaftSD);
			this.shaftPosition = shaftGuesser(feretMax, shaftCutOff, gradSlices);
		}
		
		ResultsTable rt = ResultsTable.getResultsTable();
		rt.reset();
		rt.incrementCounter();
		rt.addValue("Shaft standard deviation", this.shaftSD);
		rt.addValue("Whole bone standard deviation", this.boneSD);
		rt.addValue("Gradient standard deviation", this.gradSD);
		rt.addValue("Shaft start slice", this.shaftPosition[0]);
		rt.addValue("Shaft end slice", this.shaftPosition[1]);
		rt.addValue("Mean MFD", this.meanMFD);
		rt.addValue("Threshold min", min);
		rt.addValue("Threshold max", max);
		rt.show("Results");
		
		if(doGraph) {
			/** Plot maximum Feret diameter versus slice number */
			Plot feretPlot = new Plot("Feret Plot", "Slice", "Max Feret Diameter", this.sliceCol, this.feretMax);
			feretPlot.show();
			
			Plot gradMFDPlot = new Plot("Gradient (MFD)", "Slice", "Gradient (over " + gradSlices + " slices)", this.boneGradient[0], this.boneGradient[1]);
			gradMFDPlot.show();
			
			Plot ellipsePlot = new Plot("Feret Max / Feret Min", "Slice", "Proportion", this.sliceCol, this.feretMaxOnMin);
			ellipsePlot.show();
			
			Plot perimeterPlot = new Plot("Perimeter", "Slice", "Perimeter", this.sliceCol, this.perimeter);
			perimeterPlot.show();
			
			Plot gradPerimPlot = new Plot("Gradient (Perimeter)", "Slice", "Gradient (over " + gradSlices + " slices)", this.perimeterGradient[0], this.perimeterGradient[1]);
			gradPerimPlot.show();
		}
		
		return;
	}
	
	/**
	 * Attempt to guess at which slices the shaft begins and ends.
	 * Cycles through slices from the centre, until a condition is met (shaftLimit).
	 * 
	 * @param boneStack
	 * @param shaftLimit value boneStack[i] should reach before end of shaft is recorded.
	 * @param slices set to 1 by default: change if boneStack is a fraction of the full bone, eg. 3 if it's 1/3.
	 */
	private int[] shaftGuesser(double boneStack[], double shaftLimit, int slices) {
		
		this.shaftPosition = new int[2];
		
		// Find the central slice. (For stacks with odd number of slices, should round up ImageStackSize/2.)
		this.centralSlice = Math.round(boneStack.length / 2);
		
		// Cycle through slices from central slice, first down; then up.
		for(int k = centralSlice; k >= this.startSlice; k--) {
			if(Math.abs(boneStack[k]) >= shaftLimit) {
				this.shaftPosition[0] = k * slices;
				break;
			}
		}
		for(int k = centralSlice; k < this.iss; k++) {
			if(Math.abs(boneStack[k]) >= shaftLimit) {
				this.shaftPosition[1] = k * slices;
				break;
			}
		}
		
		return shaftPosition;
	}
	
	/**
	 * Crop the ends of the image
	 * 
	 * @param imp
	 * @param min
	 * @param max
	 */
//	private void cropEnds(ImagePlus imp, int min, int max) {
//		
//		GenericDialog gd = new GenericDialog("Crop ends");
//		gd.addMessage("Based on ");
//		gd.showDialog();
//		
//		this.centralSlice = Math.round(imp.getImageStackSize() / 2);
//		
//		for(int i = 1; i <= centralSlice; i++) {
//			
//			if(emptySlices[i]) {
//				
//				
//			}
//		}
//	}
	
	/**
	 * Find the gradient of a double[] by comparing y[i] with y[i+gradSlices]
	 * 
	 * @param x
	 * @param y
	 * @param a number of slices over which to find the gradient: resolution
	 * @return gradient[x][y]
	 */
	private double[][] gradient(double[] x, double[] y, int a) {
		
		int xCount = Math.round(x.length / a);
		double[][] gradient = new double[2][xCount - 1];
		
		for(int s = 0; s < xCount - 1; s++) {
			gradient[0][s] = s + 1;
			gradient[1][s] = (y[(s * a) + a] - y[s * a]) / a;
		}
		
		return gradient;
	}
	
	/**
	 * Remove the NaNs in a double[] and replace with zeros
	 * 
	 * @param a
	 * @return
	 */
	public static double[] removeNaNs(double[] a) {
		
		for(int i = 0; i < a.length; i++) {
			if(Double.isNaN(a[i])) {
				a[i] = 0;
			}
		}
		return a;
	}
	
	/**
	 * Calculate mean of a double[]
	 */
	public static double mean(double[] a) {
		
		double sum = 0;
		for(int i = 0; i < a.length; i++) {
			sum += a[i];
		}
		double mean = sum / a.length;
		
		return mean;
	}
	
	/**
	 * Calculates variance of a double[]
	 * 
	 * @param a
	 * @return var
	 */
	public static double variance(double[] a) {
		
		double mean = mean(a);
		
		double ssdm = 0;
		for(int i = 0; i < a.length; i++) {
			ssdm += (a[i] - mean) * (a[i] - mean);
		}
		
		double variance = ssdm / a.length;
		return variance;
	}
	
	/**
	 * Check if a number is odd, based on bit
	 */
	public boolean isOdd(int num) {
		
		if((num&1) == 1) {
			return true;
		}
		else {
			return false;
		}
	}
}
