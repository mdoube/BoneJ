package org.doube.bonej;

import java.awt.Frame;
import java.util.Arrays;

import org.doube.geometry.Centroid;
import org.doube.util.DeleteSliceRange;
import org.doube.util.ImageCheck;
import org.doube.util.ThresholdGuesser;

import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.gui.Plot;
import ij.measure.Calibration;
import ij.measure.ResultsTable;
import ij.plugin.PlugIn;

/**
 * <p>
 * Make an informed guess at which slices of an image stack contain the
 * proximal and distal "ends" of a femoral shaft.
 * <br><br>
 * Currently requires 8-bit image stacks. May want to align the bone first.
 * </p>
 * 
 * @author Nick Powell
 *
 */

public class ShaftGuesser implements PlugIn {
	
	private Calibration cal;
	private String units, bone, title;
	/** Trabeculae feature size range, 'min-max', in 'unit's */
	private String trabRange;
	
	/** Flip the z-axis if this is false. May be possible to automate detection, 
	 * as distal end of femur generally has a greater volume */
	private boolean proximalLow;
	/** Looks for a longer diaphysis if true */
	private boolean beGenerous;
	private boolean doGraph = false;
	
	private int al, startSlice, endSlice, ss, iss, boneID;
	/** Smooth over this number of slices. Initially set to 1/50th of stack size */
	private int smoothOver;
	/** Calculate the gradient over this number of slices */
	private int gradientOver;
	/** Possible first and last slice numbers of the femoral shaft, from perimeter, eccentricity, mean cortical thickness (2D) */
	private int[][] shaftPositions;
	/** Lower and upper ranges in which the shaft begins and ends */
	private int[][] shaftRange;
	/** First and last slice numbers of the femoral shaft */
	private int[] shaftPosition = new int[2];
	
	/** Measured length, in 'unit's */
	private double boneLength, shaftLength;
	/** Pixel dimensions */
	private double vW, vH, vD;
	/** Trabeculae feature sizes, in 'unit's */
	private double trabMax, trabMin;
	/** Median values */
	private double mEccentricity, mMeanCort, mFeretMax, mFeretMin, mPerimeter, mTrabeculae;
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
	/** List of number of trabeculae features per slice */
	private double[] trabeculae;
	/** List of gradients at each slice, based on values +/- from that slice */
	private double[] gFeretMax, gFeretMin, gPerimeter, gEccentricity, gMeanCort;
	/** List of smoothed values */
	private double[] sFeretMax, sFeretMin, sPerimeter, sEccentricity, sMeanCort, sTrabeculae;
	
	public void run(String arg) {
		
		/* ImageJ version checking */
		if (!ImageCheck.checkEnvironment())
			return;
		ImagePlus imp = IJ.getImage();
		if (null == imp) {
			IJ.noImage();
			return;
		}
		
		/* Seems to be necessary for everything; certainly necessary for AnalyzeParticles */
		if(imp.getBitDepth() != 8) {
			IJ.run("8-bit");
		}
		
		/* Copied from SliceGeometry */
		this.boneID = BoneList.guessBone(imp);
		this.cal = imp.getCalibration();
		this.vW = cal.pixelWidth;
		this.vH = cal.pixelHeight;
		this.vD = cal.pixelDepth;
		this.units = cal.getUnits();
		this.al = imp.getStackSize() + 1;
		this.iss = imp.getImageStackSize();
		
		double[] thresholds = ThresholdGuesser.setDefaultThreshold(imp);
		double min = thresholds[0];
		double max = thresholds[1];
		
		this.trabMin = vW * imp.getWidth() / 200;
		this.trabMax = trabMin * 40;
		
		GenericDialog gd = new GenericDialog("Shaft Guesser Options");
		String[] bones = BoneList.get();
		gd.addChoice("Bone: ", bones, bones[boneID]);
		
		gd.addCheckbox("Proximal end of bone has lower slice number", true);
		gd.addNumericField("Bone start slice", 1, 0);
		gd.addNumericField("Bone end slice", this.iss, 0);
		gd.addMessage("Choose how to estimate the shaft parameters: ");
		gd.addCheckbox("Look for longer shaft", false);
		gd.addNumericField("Smooth over +/-", Math.round(this.al / 50), 0, 3, "slices");
		gd.addNumericField("Calculate gradient over +/-", Math.round(this.al / 50), 0, 3, "slices");
		gd.addMessage("Choose trabeculae parameters: ");
		gd.addNumericField("Min particle size", trabMin, 0, 3, units);
		gd.addNumericField("Max particle size", trabMax, 0, 3, units);
		gd.addCheckbox("Graph output", doGraph);
		gd.showDialog();
		
		this.bone = gd.getNextChoice();
		this.boneID = BoneList.guessBone(bone);
		
		this.proximalLow = gd.getNextBoolean();
		this.startSlice = (int) gd.getNextNumber();
		this.endSlice = (int) gd.getNextNumber();
		this.beGenerous = gd.getNextBoolean();
		this.smoothOver = (int) gd.getNextNumber();
		this.gradientOver = (int) gd.getNextNumber();
		this.trabMin = gd.getNextNumber();
		this.trabMax = gd.getNextNumber();
		this.trabRange = trabMin + "-" + trabMax;
		this.doGraph = gd.getNextBoolean();
		if (gd.wasCanceled())
			return;
		
		/* Initial setup */
		if(!proximalLow) {
			IJ.run("Flip Z");
		}
		
		/* Crop ends (copy image to a new stack first (not yet)) */
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
			this.iss = imp.getImageStackSize();
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
		/*Ideal addition would be: cortical vs. trabecular fraction */
		
		/* Trabecular calculations */
//		this.trabeculae = quanitfyTrabeculae(imp, trabMin, trabMax);
		
		/* Calculate the eccentricities of each slice */
		this.eccentricity = eccentricity(feretMin, feretMax);
		
		/* Smooth out noise */
		this.sEccentricity = smooth(eccentricity, smoothOver);
		this.sFeretMax = smooth(feretMax, smoothOver);
		this.sFeretMin = smooth(feretMin, smoothOver);
		this.sPerimeter = smooth(perimeter, smoothOver);
		this.sMeanCort = smooth(meanCortThick2D, smoothOver);
//		this.sTrabeculae = smooth(trabeculae, smoothOver);
		
		/* Calculate 'gradients' (changes over a number of slices), preferably after smoothing */
		this.gEccentricity = gradient(sEccentricity, gradientOver);
		this.gFeretMax = gradient(sFeretMax, gradientOver);
		this.gFeretMin = gradient(sFeretMin, gradientOver);
		this.gPerimeter = gradient(sPerimeter, gradientOver);
		this.gMeanCort = gradient(sMeanCort, gradientOver);
		
		/* Calculate median values */
		this.mEccentricity = median(sEccentricity);
		this.mFeretMax = median(sFeretMax);
		this.mFeretMin = median(sFeretMin);
		this.mPerimeter = median(sPerimeter);
		this.mMeanCort = median(sMeanCort);
		
		/* Perimeter measures relative size;
		 * Cortical thickness is an actual bone property;
		 * Eccentricity measures shape.
		 * 
		 * Poss: if we increase the median cutoff by 10%, does this increase the selection length by more than 5% of the bone?
		 * Poss: don't count slices with NaN...?
		 * Poss: gradient for reliability: low is no; high is yes.
		 * 
		 * The gradients actually tell me something about how significant the features at that particular point are.
		 * 
		 * Peak detection for lesser and greater trochanter?
		 */
		
		/* Reset start and end slices based on outer limits calculated above */
//		this.startSlice = shaftPosition[0];
//		this.endSlice = shaftPosition[1];
		
		/* Run 1: possible outer limits of shaft */
		this.shaftPositions = new int[4][2];
		this.shaftPositions[0] = shaftLimiter(sPerimeter, mPerimeter, false);
		this.shaftPositions[1] = shaftLimiter(sEccentricity, mEccentricity, false);
		this.shaftPositions[2] = shaftLimiter(sMeanCort, mMeanCort, true);
//		this.shaftPositions[3] = shaftLimiter(sTrabeculae, mTrabeculae, false);
		
		/* Run 2: refine shaft limits: by shape */
//		if(this.eccentricity[shaftPosition[0]] > mEccentricity) {
//			
//		}
		
		/* Check for broken values: proximal greater than distal */
		for(int i = 0; i < shaftPositions.length; i++) {
			this.shaftPositions[i] = checkShaftSanity(shaftPositions[i]);
		}
		
		/* Attempt to get over notches in the shaft.
		 * If we're being generous, increase or decrease the median by 10% 
		 * - if this increases the shaft length by >= 10%, keep the new values */
//		if(beGenerous) {
//			this.mPerimeter = this.mPerimeter * 1.1;
//			this.mEccentricity = this.mEccentricity * 1.1;
//			this.mMeanCort = this.mMeanCort * 1.1;
//			
//			if(shaftLimiter(sPerimeter, mPerimeter, false)[1] - shaftLimiter(sPerimeter, mPerimeter, false)[0] >= 1.1 * (this.shaftPositions[0][1] - this.shaftPositions[0][0])) {
//				this.shaftPositions[0] = shaftLimiter(sPerimeter, mPerimeter, false);
//			}
//			if(shaftLimiter(sEccentricity, mEccentricity, false)[1] - shaftLimiter(sEccentricity, mEccentricity, false)[0] >= 1.1 * (this.shaftPositions[1][1] - this.shaftPositions[1][0])) {
//				this.shaftPositions[1] = shaftLimiter(sEccentricity, mEccentricity, false);
//			}
//			if(shaftLimiter(sMeanCort, mMeanCort, false)[1] - shaftLimiter(sMeanCort, mMeanCort, false)[0] >= 1.1 * (this.shaftPositions[2][1] - this.shaftPositions[2][0])) {
//				this.shaftPositions[2] = shaftLimiter(sMeanCort, mMeanCort, false);
//			}
//		}
		
		/* Switch columns and rows for inRange, below */
		this.shaftRange = new int[shaftPositions[0].length][shaftPositions.length];
		for(int i = 0; i < shaftRange.length; i++) {
			for(int j = 0; j < shaftRange[0].length; j++) {
				shaftRange[i][j] = shaftPositions[j][i];
			}
		}
		
		/* inRange tests for 'close by' numbers, and returns the pair it finds.
		 * beGenerous: gives min value (lower slice num) of proximal end and 
		 * max value of distal end (i.e. longer shaft); otherwise vice versa.
		 * Could also average the pair. */
		if(beGenerous) {
			this.shaftPosition[0] = inRange(shaftRange[0], smoothOver)[0];
			this.shaftPosition[1] = inRange(shaftRange[1], smoothOver)[1];
		}
		else {
			this.shaftPosition[0] = inRange(shaftRange[0], smoothOver)[1];
			this.shaftPosition[1] = inRange(shaftRange[1], smoothOver)[0];
		}
		
		/* If inRange fails to find a pair, average the values from perimeter, etc. */
		if(shaftPosition[0] == 0) {
			shaftPosition[0] = Centroid.getCentroid(shaftRange[0]);
		}
		if(shaftPosition[1] == 0) {
			shaftPosition[1] = Centroid.getCentroid(shaftRange[1]);
		}
		
		this.shaftLength = (shaftPosition[1] - shaftPosition[0]) * vD;
		this.boneLength = this.iss * vD;
		
		ResultsTable rt = ResultsTable.getResultsTable();
		rt.incrementCounter();
		rt.addLabel("Title", title);
		rt.addValue("Bone Code", boneID);
		
		rt.addValue("Median Eccentricity", mEccentricity);
		rt.addValue("Median Feret Max", mFeretMax);
		rt.addValue("Median Feret Min", mFeretMin);
		rt.addValue("Median Perimeter", mPerimeter);
		rt.addValue("Median Mean Cortical Thickness (2D)", mMeanCort);
		
		rt.addValue("Shaft start slice (Perim)", shaftPositions[0][0]);
		rt.addValue("Shaft end slice (Perim)", shaftPositions[0][1]);
		rt.addValue("Shaft start gradient (Perim)", gPerimeter[shaftPositions[0][0]]);
		rt.addValue("Shaft end gradient (Perim)", gPerimeter[shaftPositions[0][1]]);
		
		rt.addValue("Shaft start slice (Ecc)", shaftPositions[1][0]);
		rt.addValue("Shaft end slice (Ecc)", shaftPositions[1][1]);
		rt.addValue("Shaft start gradient (Ecc)", gEccentricity[shaftPositions[1][0]]);
		rt.addValue("Shaft end gradient (Ecc)", gEccentricity[shaftPositions[1][1]]);
		
		rt.addValue("Shaft start slice (Cort)", shaftPositions[2][0]);
		rt.addValue("Shaft end slice (Cort)", shaftPositions[2][1]);
		rt.addValue("Shaft start gradient (Cort)", gMeanCort[shaftPositions[2][0]]);
		rt.addValue("Shaft end gradient (Cort)", gMeanCort[shaftPositions[2][1]]);
		
		rt.addValue("Proximal shaft end", shaftPosition[0]);
		rt.addValue("Distal shaft end", shaftPosition[1]);
		rt.addValue("Diaphysis length (" + units + ")", shaftLength);
		rt.addValue("Bone length (" + units + ")", boneLength);
		rt.addValue("Diaphysis proportion", shaftLength / boneLength);
		
		rt.show("Results");
		
		if(doGraph) {
			
			/* For plotting */
			this.slices = new double[imp.getStackSize() + 1];
			for (int s = this.startSlice; s <= this.iss; s++) {
				slices[s] = (double) s;
			}
			
			Plot sEccPlot = new Plot("Smoothed eccentricity", "Slice", "Eccentricity", this.slices, this.sEccentricity);
			sEccPlot.show();
			
			Plot sPerimPlot = new Plot("Smoothed perimeter", "Slice", "Perimeter", this.slices, this.sPerimeter);
			sPerimPlot.show();
			
			Plot sMeanCortPlot = new Plot("Smoothed Mean Cort Thick (2D)", "Slice", "Mean Cort Thickness", this.slices, this.sMeanCort);
			sMeanCortPlot.show();
			
			Plot gEccPlot = new Plot("Gradient eccentricity", "Slice", "Eccentricity", this.slices, this.gEccentricity);
			gEccPlot.show();
			
			Plot gPerimPlot = new Plot("Gradient perimeter", "Slice", "Perimeter", this.slices, this.gPerimeter);
			gPerimPlot.show();
			
			Plot gMeanCortPlot = new Plot("Gradient Mean Cort Thick (2D)", "Slice", "Mean Cort Thickness", this.slices, this.gMeanCort);
			gMeanCortPlot.show();
		}
	}
	
	public int[] getShaftPosition() {
		return this.shaftPosition;
	}
	
	/**
	 * Estimates at which slices the shaft begins and ends, based on numerical limits.
	 * Cycles through slices from the centre, first down, then up, until 
	 * a condition is met (shaftLimit).
	 * Main weaknesses: assumes the central slice is part of the shaft. Also assumes
	 * a distinctive "u-shaped" profile of the shaft (or inverse when isMin is true).
	 * 
	 * @param boneStack
	 * @param shaftLimit
	 * @param isMin Specify whether the shaftLimit value is a minimum (true) or a maximum (false)
	 * @return
	 */
	private int[] shaftLimiter(double[] boneStack, double shaftLimit, boolean isMin) {
		
		int[] shaftEnds = new int[2];
		/* Find the central slice. (For stacks with odd number of slices, 
		 * should round up ImageStackSize/2.)
		 * An improvement would be to start at the slice containing the centroid
		 * of the bone, but slower to run.
		 */
		int centralSlice = Math.round((this.endSlice - this.startSlice) / 2);
		
		// Cycle through slices from central slice, first down; then up.
		for(int i = centralSlice; i >= this.startSlice; i--) {
			if(isMin && Math.abs(boneStack[i]) < shaftLimit) {
				shaftEnds[0] = i + 1;
				break;
			}
			else if(!isMin && Math.abs(boneStack[i]) > shaftLimit) {
				shaftEnds[0] = i + 1;
				break;
			}
		}
		for(int i = centralSlice; i < this.endSlice; i++) {
			if(isMin && Math.abs(boneStack[i]) < shaftLimit) {
				shaftEnds[1] = i - 1;
				break;
			}
			else if(!isMin && Math.abs(boneStack[i]) > shaftLimit) {
				shaftEnds[1] = i - 1;
				break;
			}
		}
		
		return shaftEnds;
	}
	
	/**
	 * Checks whether a shaft's start and end points are the wrong way around.
	 * 
	 * @param a
	 * @return b, empty array, if so. (Else return input array a).
	 */
	private int[] checkShaftSanity(int[] a) {
		
		if(a[0] > a[1]) {
			int[] b = new int[a.length];
			return b;
		}
		else return a;
	}
	
	/* 
	 * Want: maximum ROI ie outer region of bone
	 * - could say: here's the centroid of the shaft
	 * - here's the average shaft diameter or max over whole shaft,
	 * - ROI is circular region with this diameter from centroid
	 * Ideally though, select outer edge of bone, or even mid-cortical(!?) as ROI border.
	 * 
	 * Need to guess at size of trabeculae relative to size of image or bone iself
	 * 
	 * Actually want: count of  holes formed by trabeculae ("shafts"? "connecting bits")
	 * on each slice.
	 * Plot by slice. See what it looks like.
	 * 
	 * Possible? to then draw an outer region around all the points it finds, so:
	 * - within this region: trabeculae;
	 * - nb medullary canal: use distance between points? so want doughnut around medullary canal
	 * - then compare total cortical area with total trabeculae area... 
	 *
	 * */
	/**
	 * Quantify the trabecular particles in each slice.
	 * Requires a binary image.
	 * 
	 * Flow is:
	 * Copy to new image
	 * Set 8-bit
	 * Binarise (threshold)
	 * Watershed
	 * Set min& max particle sizes
	 * Analyse
	 * 
	 * @param imp
	 * @param trabMin
	 * @param trabMax
	 * @return list of numbers of trabecular particles, within the given range, per slice.
	 */
	private double[] quanitfyTrabeculae(ImagePlus imp, double trabMin, double trabMax) {
		
		
		
		double[] trabList = new double[imp.getStackSize()];
		
		/* Gives results but can't use them */
//		IJ.run("Analyze Particles...", "size=" + trabMin + "-" + trabMax + "circularity=0.00-1.00 show=Masks display clear include summarize stack");
		
//		Frame frame = WindowManager.getFrame("Summary");
//		IJ.saveAs("Text", "\\home\\ij_results\\Summary of femur.xls");
		return trabList;
	}
	
	/**
	 * Estimate of eccentricity, 0 < e < 1 (0 is a circle), for a pair of doubles[]
	 * Based on Feret's min and max, so this is not *true* eccentricity,
	 * which requires the axes to be perpendicular (true anyway in a true ellipse).
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
	public static double[] smooth(double[] x, int a) {
		
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
	 * Tests for 'close' number pairs. Discover whether two ints are within 
	 * a given range. Will theoretically work for any int[] length.
	 * Currently finds the pair furthest apart (constrained to the range).
	 * 
	 * @param a
	 * @param diff, the inclusive range
	 * @return sortedList, the furthest-apart pair within the range, lowest first; zeros if none.
	 */
	private int[] inRange(int[] a, int diff) {
		
		int[] sortedList = new int[2];
		int[] z = a.clone();
		
		// Sort first as will break out of loop on finding the pair widest apart.
		Arrays.sort(z);
		
		for(int i = 0; i < z.length; i++) {
			for(int j = 0; j < (z.length - (i + 1)); j++) {
				if(z[z.length - (j+1)] - z[i] <= diff
						&& z[z.length - (j+1)] - z[i] >= sortedList[1] - sortedList[0]) {
					
					sortedList[0] = z[i];
					sortedList[1] = z[z.length - (j+1)];
				}
			}
		}
		return sortedList;
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
	
	/**
	 * Calculate the variance of a double[]
	 * 
	 * @param a
	 * @return variance
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
}
