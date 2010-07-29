package org.doube.bonej;

import java.util.Arrays;

import org.doube.bonej.SliceGeometry;
import org.doube.geometry.Centroid;
import org.doube.util.ImageCheck;
import org.doube.util.ThresholdGuesser;

import com.sun.xml.internal.bind.v2.runtime.unmarshaller.XsiNilLoader.Array;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.gui.Plot;
import ij.measure.ResultsTable;
import ij.plugin.PlugIn;

/** 
 * <p>
 * Isolates the shaft numerically, using various methods:
 * Uses org.doube.bonej.SliceGeometry to calculate maximum Feret diameter (MFD)/slice.
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
	
	private boolean proximalLow, doGraph, cropEnds, do3DMarkers, doManuallyThreshold;
	private boolean[] emptySlices;
	
	private double shaftSD, boneSD, gradSD, meanMFD, shaftCutOff, min, max, medianCort;
	private double[] slices, feretMax, feretMin, shaft, boneGradientNoNaNs, feretMinOnMax, eccentricity, perimeter, meanCortThick2D, sortedMCT2D, smoothGradE;
	private double[][] boneGradient, perimeterGradient;
	
	private int al, startSlice, iss, centralSlice, gradSlices, smoothSlices;
	private int[] shaftPosition;			// Contains first and last shaft slices
	
	private String opt_0 = "All";
	private String opt_1 = "0.25 x standard deviation of gradient of MFD";
	private String opt_2 = "Mean MFD + 0.5 standard deviations of MFD";
	private String opt_3 = "Shaft has above median mean cortical thickness (2D)";
	private String[] shaftViaChoices = {opt_0, opt_1, opt_2, opt_3};

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
		gd.addChoice("Shaft selection via", shaftViaChoices, opt_0);
		gd.addNumericField("Calculate gradient over # slices: ", 1, 0);
		gd.addNumericField("Smooth over # slices: ", 5, 0);
//		gd.addCheckbox("Show 3D Markers at end of shaft", true);
		gd.addCheckbox("Manually threshold", false);
		gd.addNumericField("Manually threshold: min", 30, 0);
		gd.addNumericField("Manually threshold: max", 255, 0);
		gd.showDialog();
		
		this.proximalLow = gd.getNextBoolean();
		this.doGraph = gd.getNextBoolean();
//		this.cropEnds = gd.getNextBoolean();
		String shaftVia = gd.getNextChoice();
		this.gradSlices = (int) gd.getNextNumber();
		this.smoothSlices = (int) gd.getNextNumber();
//		this.do3DMarkers = gd.getNextBoolean();
		this.doManuallyThreshold = gd.getNextBoolean();
		this.min = gd.getNextNumber();
		this.max = gd.getNextNumber();
		
		if(!proximalLow) {
			IJ.run("Flip Z");
		}
		
		if(!doManuallyThreshold) {
			// Copied from SliceGeometry
			double[] thresholds = ThresholdGuesser.setDefaultThreshold(imp);
			this.min = thresholds[0];
			this.max = thresholds[1];
		}
		
		this.al = imp.getStackSize() + 1;
		this.startSlice = 1;
		this.iss = imp.getImageStackSize();
		
		this.slices = new double[this.al];
		for (int s = this.startSlice; s <= this.iss; s++) {
			slices[s] = (double) s;
		}
		
		// Get values from SliceGeometry
		SliceGeometry sg = new SliceGeometry();
		sg.setParameters(imp);
		sg.calculateCentroids(imp, min, max);
		sg.roiMeasurements(imp, min, max);
		sg.calculateMoments(imp, min, max);
		sg.calculateThickness2D(imp, min, max);
		this.feretMax = sg.getFeretMax();
		this.feretMin = sg.getFeretMin();
		this.perimeter = sg.getPerimeter();
		this.emptySlices = sg.getEmptySlices();
		this.meanCortThick2D = sg.getMeanCorticalThickness2D();
		
		// Median of means of cortical thickness
//		this.medianCort = median(meanCortThick2D);
		
		// Feret Min / Feret Max
		this.feretMinOnMax = new double[feretMax.length];
		for(int i = 0; i < feretMax.length; i++) {
			feretMinOnMax[i] = feretMin[i] / feretMax[i] ;
		}
		
		// Rough eccentricity, based on Feret's min and max
		this.eccentricity = new double[feretMax.length];
		for(int i = 0; i < feretMax.length; i++) {
			eccentricity[i] = eccentricity(feretMin[i], feretMax[i]);
		}
		
		// Gradients
		this.boneGradient = gradient(slices, feretMax, gradSlices);
		this.perimeterGradient = gradient(slices, perimeter, gradSlices);
		this.smoothGradE = smooth(eccentricity, smoothSlices);
		
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
		if(shaftVia == opt_1) {
			
			this.shaftCutOff = 0.25 * this.gradSD;
			this.shaftPosition = shaftGuesser(boneGradientNoNaNs, shaftCutOff, gradSlices);
		}
		
		/* Uses Mean(max Feret diameter) and Standard Deviation(max FD).
		 * Also relies on arbitrary 'mid third', using shaftSD.
		 * Measurement is Mean + 1/2 SD */
		if(shaftVia == opt_2) {
			
			this.shaftCutOff = Centroid.getCentroid(feretMax) + (0.5 * this.shaftSD);
			this.shaftPosition = shaftGuesser(feretMax, shaftCutOff, gradSlices);
		}
		
		if(shaftVia == opt_3) {
			
			this.shaftCutOff = medianCort;
			this.shaftPosition = shaftGuesser(meanCortThick2D, shaftCutOff, 1);
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
		rt.addValue("Median Mean Cortical Thickness 2D", medianCort);
		rt.addValue("SmGrad", smoothGradE[60]);
		rt.show("Results");
		
		if(doGraph) {
			/** Plot maximum Feret diameter versus slice number */
			Plot feretPlot = new Plot("Feret Plot", "Slice", "Max Feret Diameter", this.slices, this.feretMax);
			feretPlot.show();
			
			Plot gradMFDPlot = new Plot("Gradient (MFD)", "Slice", "Gradient (over " + gradSlices + " slices)", this.boneGradient[0], this.boneGradient[1]);
			gradMFDPlot.show();
			
			Plot ellipsePlot = new Plot("Feret Min / Feret Max", "Slice", "Proportion", this.slices, this.feretMinOnMax);
			ellipsePlot.show();
			
			Plot meanCortThick2DPlot = new Plot("Mean Cortical Thickness (2D)", "Slice", "Thickness", this.slices, this.meanCortThick2D);
			meanCortThick2DPlot.show();
			
			Plot perimeterPlot = new Plot("Perimeter", "Slice", "Perimeter", this.slices, this.perimeter);
			perimeterPlot.show();
			
			Plot gradPerimPlot = new Plot("Gradient (Perimeter)", "Slice", "Gradient (over " + gradSlices + " slices)", this.perimeterGradient[0], this.perimeterGradient[1]);
			gradPerimPlot.show();
			
			Plot ePlot = new Plot("Eccentricity", "Slice", "e", this.slices, this.eccentricity);
			ePlot.show();
			
			Plot smoothGradEPlot = new Plot("Smoothed eccentricity", "Slice", "e (over " + smoothSlices + " slices)", this.slices, this.smoothGradE);
			smoothGradEPlot.show();
		}
		
		return;
	}
	
	/**
	 * Rough eccentricity, 0 < e < 1 (0 is a circle), based on Feret's min and max.
	 * See <a href="http://mathworld.wolfram.com/Eccentricity.html">http://mathworld.wolfram.com/Eccentricity.html</a>
	 */
	private double eccentricity(double feretMin, double feretMax) {
		
		double e = Math.sqrt(1 - ((feretMin * feretMin) / (feretMax * feretMax)));
		return e;
	}
	
	/**
	 * Attempt to guess at which slices the shaft begins and ends.
	 * Cycles through slices from the centre, until a condition is met (shaftLimit).
	 * 
	 * @param boneStack
	 * @param shaftLimit value boneStack[i] should reach before end of shaft is recorded.
	 * @param numSlices set to 1 by default: change if boneStack is a fraction of the full bone, eg. 3 if it's 1/3.
	 */
	private int[] shaftGuesser(double boneStack[], double shaftLimit, int numSlices) {
		
		this.shaftPosition = new int[2];
		
		// Find the central slice. (For stacks with odd number of slices, should round up ImageStackSize/2.)
		this.centralSlice = Math.round(boneStack.length / 2);
		
		// Cycle through slices from central slice, first down; then up.
		for(int k = centralSlice; k >= this.startSlice; k--) {
			if(Math.abs(boneStack[k]) >= shaftLimit) {
				this.shaftPosition[0] = k * numSlices;
				break;
			}
		}
		for(int k = centralSlice; k < this.iss; k++) {
			if(Math.abs(boneStack[k]) >= shaftLimit) {
				this.shaftPosition[1] = k * numSlices;
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
	 * Find the gradient of a double[] by averaging over the surrounding slices.
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
	 * Calculates variance of a double[]
	 * 
	 * @param a
	 * @return var
	 */
	public static double variance(double[] a) {
		
		double mean = Centroid.getCentroid(a);
		
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
	
	/**
	 * Find the median of a double[]
	 */
	public static double median(double[] a) {
		
		int mid = a.length / 2;
		Arrays.sort(a);
		
		if(a.length%2 == 1) {
			return a[mid];
		}
		else {
			return (a[mid - 1] + a[mid]) / 2;
		}
	}
}
