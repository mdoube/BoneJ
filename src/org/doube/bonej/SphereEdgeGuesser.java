package org.doube.bonej;

import java.awt.AWTEvent;
import java.awt.Checkbox;
import java.awt.Color;
import java.awt.TextField;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.Vector;

import org.doube.geometry.Centroid;
import org.doube.geometry.FitEllipsoid;
import org.doube.geometry.FitSphere;
import org.doube.geometry.Trig;
import org.doube.geometry.Vectors;
import org.doube.util.DialogModifier;
import org.doube.util.ImageCheck;
import org.doube.util.ThresholdGuesser;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.gui.ImageCanvas;
import ij.gui.ImageWindow;
import ij.gui.Plot;
import ij.gui.WaitForUserDialog;
import ij.measure.Calibration;
import ij.measure.ResultsTable;
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;

/**
 * <p>Aims to drastically reduce the complexity of, and time taken to, 
 * fit a sphere to a sphere-like object such as a femoral head.  Requires 
 * one manual input point inside the sphere.</p>
 * 
 * <p>Potential issues: false positives such as holes, well-defined trabeculae, etc.
 * There is some attempt to account for these, as well as image artefacts.</p>
 * 
 * <p>Can switch resolution to 'unit's (generally mm), although fitSphere likes 
 * the current format and is picky.</p>
 * 
 * <p>Arbitrary elements: Along individual vectors, boundaries are found using 
 * medians of pixel values along that vector - the limits for these are in 
 * the boundaryLimiter method.</p>
 * 
 * <p>If the resultant sphere has an unlikely Z coordinate (slice), fitSphere 
 * may not have received a wide enough range of z coordinates.</p>
 * 
 * @author Nick Powell
 *
 */
public class SphereEdgeGuesser implements PlugIn, DialogListener, MouseListener {
	
	private ImageProcessor[] sliceProcessors = null;
	private ImageCanvas canvas;
	private Calibration cal;
	
	/** Show the final points used to generate the sphere on a new image */
	private boolean showPoints = false;
	/** Show a sample plot of getPixelValue at each point along a vector */
	private boolean doPlot = false;
	/** Use the centre found from run 1 as the initial point for run 2 (etc.) */
	private boolean doRecursion = true;
	/** Use method of looking away from the centroid of the bone */
	private boolean ignoreCentroidDirection = false;
	/** For HU units (copied from Neck Shaft Angle) */
	private boolean fieldUpdated = false;
	/** Exclude vectors pointing towards the centroid of the entire bone */
	private boolean useBoneCentroid = true;
	/** Also exclude vectors pointing towards centroid of slice containing initialPoint */
	private boolean useSliceCentroid = false;
	/** Attempt to fit an ellipsoid rather than a sphere. Uses FitEllipsoid.yuryPetrov(points). */
	private boolean fitEllipsoid = false;
	
	/** Linear unit of measure */
	private String units;
	/** Pixel dimensions in 'unit's */
	private double vH, vW, vD;
	/** Unit vector sizes, -1 < value < 1 */
	private double uX, uY, uZ;
	/** Unit vectors of the vector joining initialPoint - Centroid */
	private double[] uCW, uCS;
	/** Standard deviation values of the centre point and radius */
	private double[] sd;
	/** Exclude vectors within this multiple of pi radians from the vector joining the initial point to the centroid of the entire bone. (0.5 is a hemisphere.) */
	private double excludeWholeC = 0.5;
	/** Exclude vectors within this multiple of pi radians from the vector joining the initial point to the centroid of the slice containing the initial point. (0.5 is a hemisphere.) */
	private double excludeSliceC = 0.75;
	
	/** Initial (iX), current (nX) coordinates (in pixels or 'unit's) */
	private int iX, iY, iZ, nX, nY, nZ;
	/** Image dimensions in 'unit's */
	private int d, h, w;
	/** Number of times to re-run the code */
	private int runNum = 100;
	/** Number of vectors to create */
	private int numVectors = 500;
	/** Weight uZ by this factor */
	private double bias_uZ = 0;
	/** Multiple of standard deviation from the mean, within which the 3D points 
	 * may fluctuate their values +/-. If the code finds too few points, increase 
	 * this value. Default is 1; sometimes set to 2. */
	private double sd3DMult = 1;
	/** Multiple of standard deviation from the mean to subtract from the mean 
	 * in order to limit the radius of the points in 2D. Default is 0.3: increase 
	 * to restrict rules further. */
	private double sd2DMult = 0.3;
	
	/** Pixel values along a vector. Expands for vector length. */
	private ArrayList<Float> vectorPixelValues;
	/** Array (known length) of ArrayLists (differing lengths) 
	 * containing pixel values along each vector */
	private ArrayList<Float>[] pixelValues;
	
	/** Centroid of the bone. W: whole stack; S: slice containing initialPoint. */
	private double[] centroidW, centroidS;
	/** When sphere fitting, holds (x, y, z) centre and radius, in same units as those given to fitSphere. 
	 * When ellipsoid fitting, holds (x, y, z) centre and 3 radii. */
	private double[] dimensions;
	/** Fill with mean dimensions from all runs */
	private double[] meanDimensions;
	/** User's chosen start point */
	private double[] initialPoint = new double[3];
	private double[] initialPointUser = new double[3];
	/** List of median pixel values of each vector */
	private double[] medianValues;
	/** List of number of steps to bone boundary (median crossing). */
	private int[] boundarySteps;
	/** Distance in 'unit's to the bone boundary given by boundarySteps */
	private double[] boundaryDistances;
	/** List of magnitudes of each unit vector, in 'unit's  */
	private double[] magnitudes;
	/** List of distances to coordinates from a centroid */
	private double[] coOrdDistances;
	/** List of biases to apply to uZ */
	private double[] biases;
	/** Sample (get pixel values) ahead along vectors these percentages of the vector. 
	 * To jump ahead, the pixel values at the two lowest percentages must be beneath 
	 * the limit - or just the largest value. */
	private double[] sampleAheadPc = {5,10,20}; // So far, these haven't needed adjustment, but they could be, for large variations in trabeculae sizes.
	
	/** Mean of distances to each coordinate */
	private double meanCoOrdDistance;
	/** Standard deviation of all the distances */
	private double sdCoOrdDistance;
	/** Adjusted limit (from median) */
	private double adjustedLimit;
	
	/** Random unit vectors (x, y, z) */
	private double[][] unitVectors;
	/** Ragged array of pixel values for each vector */
	private double[][] pxVals;
	/** Ragged array of x-axis values for plotting pixel value profiles along vectors */
	private double[][] xValues;
	
	/** Single set of boundary coordinates x, y, z */
	private double[] coOrds;
	/** ArrayList (of length limited by meanDistance) of boundary coordinates x, y, z */
	private ArrayList<double[]> coOrdinates, newCoOrdinates;
	/** Array of boundary coordinates x, y, z. T: transposed; Ref: refined. */
	private double[][] coOrdArray;
	
	public void run(String arg) {
		
		/* Copied from SphereFitter */
		if (!ImageCheck.checkEnvironment())
			return;
		final ImagePlus imp = IJ.getImage();
		if (null == imp) {
			IJ.noImage();
			return;
		}
//		ImageCheck ic = new ImageCheck();
//		if (ic.isMultiSlice(imp)) {
//			IJ.error("Stack required for sphere fitting");
//			return;
//		}
		
		ImageWindow win = imp.getWindow();
		this.canvas = win.getCanvas();
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
		
		/* Optionally show a dialogue */
		GenericDialog gd = new GenericDialog("Options");
		gd.addNumericField("Create", numVectors, 0, 4, "vectors");
		gd.addNumericField("Run (to refine)", runNum, 0, 3, "times");
		gd.addMessage("Sphere fitting: ");
		gd.addNumericField("Limit 2D: mean - ", sd2DMult, 2, 5, "standard deviations");
		gd.addNumericField("Limit 3D: mean +/-", sd3DMult, 2 , 5, "standard deviations");
//		gd.addNumericField("Sample ahead along vector by", sampleAheadPc[0], 0, 4, "% ");
//		gd.addNumericField("and", sampleAheadPc[1], 0, 4, "% ");
//		gd.addNumericField("or", sampleAheadPc[2], 0, 4, "% ");
		gd.addCheckbox("Vectors from initial point avoid...", ignoreCentroidDirection);
		gd.addNumericField("cone +/-", excludeWholeC, 2 , 5, "pi radians from...");
		gd.addCheckbox("...whole bone centroid", useBoneCentroid);
		gd.addNumericField("cone +/-", excludeSliceC, 2 , 5, "pi radians from...");
		gd.addCheckbox("...initial slice centroid", useSliceCentroid);
		gd.addCheckbox("HU Calibrated", ImageCheck.huCalibrated(imp));
		gd.addNumericField("Bone Min:", min, 1, 6, pixUnits + " ");
		gd.addNumericField("Bone Max:", max, 1, 6, pixUnits + " ");
		gd.addCheckbox("Use recursion", doRecursion);
		gd.addCheckbox("Show points", showPoints);
		gd.addCheckbox("Plot a vector's profile", doPlot);
		gd.addCheckbox("Fit an ellipsoid instead of a sphere", fitEllipsoid);
//		gd.addDialogListener();
		gd.showDialog();
		if(gd.wasCanceled()) { return; }
		this.numVectors = (int) gd.getNextNumber();
		this.runNum = (int) gd.getNextNumber();
		this.sd2DMult = gd.getNextNumber();
		this.sd3DMult = gd.getNextNumber();
//		this.sampleAheadPc[0] = gd.getNextNumber();
//		this.sampleAheadPc[1] = gd.getNextNumber();
//		this.sampleAheadPc[2] = gd.getNextNumber();
		this.ignoreCentroidDirection = gd.getNextBoolean();
		this.useBoneCentroid = gd.getNextBoolean();
		this.excludeWholeC = gd.getNextNumber();
		this.useSliceCentroid = gd.getNextBoolean();
		this.excludeSliceC = gd.getNextNumber();
		boolean isHUCalibrated = gd.getNextBoolean();
		if (isHUCalibrated) {
			min = cal.getRawValue(min);
			max = cal.getRawValue(max);
		}
		min = gd.getNextNumber();
		max = gd.getNextNumber();
		this.doRecursion = gd.getNextBoolean();
		this.showPoints = gd.getNextBoolean();
		this.doPlot = gd.getNextBoolean();
		this.fitEllipsoid = gd.getNextBoolean();
		
		if(fitEllipsoid) {
			if(!ignoreCentroidDirection) {
				IJ.log("Switching to centroid-based method.");
				ignoreCentroidDirection = true;
			}
			if(doRecursion) {
				IJ.log("Switching off recursion.");
				doRecursion = false;
			}
		}
		
		/* User clicks to set initial point */
		getInitialPointUser(canvas);

		/* Confirm initial point and measurements */
//		GenericDialog gd1 = new GenericDialog("Info");
//		gd1.addMessage("Initial point (x,y,z): (" + iX + ", " + iY + ", " + iZ + ")");
//		gd1.addMessage("Pixel sizes (vW,vH,vD): (" + vW + ", " + vH + ", " + vD + ")");
//		gd1.addMessage("imp.getCurrentSlice: " + imp.getCurrentSlice() + "");
//		gd1.showDialog();
//		IJ.log("Pixel sizes (vW,vH,vD): " + vW + ", " + vH + ", " + vD + "");
//		IJ.log("Initial point: " + iX + ", " + iY + ", " + iZ + "");
		
		if(ignoreCentroidDirection) {
			findCentroids(imp, min, max);
		}
		
		findSphere(imp, runNum, numVectors, fitEllipsoid, ignoreCentroidDirection);

		/* Show statistical results */
		fillResultsTable("stats", meanDimensions, sd, fitEllipsoid).show("Results");
		
		annotateCentre(imp, meanDimensions).show();
		
		if(doPlot) {
			/* For plotting pixel values, which are helpful if the wrong shapes 
			 * are being included (i.e. boundaryLimiter needs adjusting). */
			this.xValues = new double[pixelValues.length][];
			for(int j = 0; j < xValues.length; j++) {
				xValues[j] = new double[pixelValues[j].size()];
				for (int k = 0; k < pixelValues[j].size(); k++) {
					xValues[j][k] = k;
				}
			}
			
			Plot aPlot = new Plot("Pixel Values along vector " + 0 + "", "Distance (pixels)", "Pixel value", xValues[0], pxVals[0]);
			aPlot.show();
		}
		
		
		return;
	}
	
	
	/**
	 * @return sphereDim[4] Sphere centre x, y, z coordinates and radius.
	 */
	public double[] getMeanDimensions() {
		return this.meanDimensions;
	}
	public double[] getStandardDeviations() {
		return this.sd;
	}
	
	/**
	 * Run this method to pop up a dialog and set up a MouseListener so the user
	 * can set an initial point from which to project vectors.
	 * 
	 * @param x
	 * @param y
	 * @param z
	 */
	public void getInitialPointUser(ImageCanvas canvas) {
		
		// required for mousePressed to work
		this.canvas = canvas;
		
		// remove stale MouseListeners
		MouseListener[] l = canvas.getMouseListeners();
		for (int n = 0; n < l.length; n++) {
			canvas.removeMouseListener(l[n]);
		}
		// add a new MouseListener
		canvas.addMouseListener(this);
		new WaitForUserDialog("Click (roughly) in the centre of\n" 
									+ "the femoral head, then hit \'OK\'").show();
		canvas.removeMouseListener(this);
	}
	
	public void mousePressed(MouseEvent e) {
		ImagePlus imp = IJ.getImage();
		int x = canvas.offScreenX(e.getX());
		int y = canvas.offScreenY(e.getY());
		int z = imp.getCurrentSlice();
		setInitialPoint(x, y, z);
		setInitialPointUser(x, y, z);
	}

	public void mouseReleased(MouseEvent e) { }
	public void mouseExited(MouseEvent e) { }
	public void mouseClicked(MouseEvent e) { }
	public void mouseEntered(MouseEvent e) { }
	public void mouseMoved(MouseEvent e) { }
	
	public void setInitialPoint(double x, double y, double z) {
		this.initialPoint[0] = x;
		this.initialPoint[1] = y;
		this.initialPoint[2] = z;
		
		this.iX = (int) Math.floor(initialPoint[0]);
		this.iY = (int) Math.floor(initialPoint[1]);
		this.iZ = (int) Math.floor(initialPoint[2]);
	}
	
	public void setInitialPointUser(double x, double y, double z) {
		this.initialPointUser[0] = x;
		this.initialPointUser[1] = y;
		this.initialPointUser[2] = z;
	}
	
	/**
	 * Method to run most of the code.
	 * 
	 * @param imp
	 * @param runNum
	 * 				<br /> int, the number of times to re-project vectors and fit a
	 * 						sphere, from which to take the mean dimensions.
	 * @param numVectors
	 * 				<br /> int, the number of vectors to project during each run.
	 * 						If ignoreCentroidDirection is false, this number of 
	 * 						vectors will actually be projected in both 2D and 3D.
	 * @param fitEllipsoid
	 * 				<br /> boolean, attempt to fit an ellipsoid instead of a sphere
	 * 						(sphere is default, if this is false).
	 * @param ignoreCentroidDirection
	 * 				<br /> boolean, use method of ignoring methods in the direction
	 * 						of particular centroids.
	 */
	public void findSphere(ImagePlus imp, int runNum, int numVectors, boolean fitEllipsoid, boolean ignoreCentroidDirection) {
		
		/* Set up array of ImageProcessors for reference later */
		ImageStack stack = imp.getStack();
		this.sliceProcessors = new ImageProcessor[imp.getStackSize() + 1];
		for(int s = 1; s < sliceProcessors.length; s++) {
			sliceProcessors[s] = stack.getProcessor(s);
		}
		
		w = imp.getWidth();
		h = imp.getHeight();
		d = imp.getStackSize();
		
		cal = imp.getCalibration();
		vW = cal.pixelWidth;
		vH = cal.pixelHeight;
		vD = cal.pixelDepth;
		units = cal.getUnits();
		
		this.biases = new double[3];
		biases[0] = 0;			// Run 1: vectors in current slice only (2D)
		biases[1] = 1;			// Run 2: 3D vectors
		biases[2] = 1;			// Or: ignore vectors towards the centroid
		
		if(fitEllipsoid) {
			dimensions = new double[6];
		}
		else {
			dimensions = new double[4];
		}
		double[][] allDimensions = new double[runNum][dimensions.length];
		
		/* Run a number of times, based on either the same initial point, or 
		 * an updating one (doRecursion). */
		for(int r = 0; r < runNum; r++) {
			
//			IJ.log("Initial point: " + iX + ", " + iY + ", " + iZ + "");
			
			/* Results from 2D, relying upon the user's guess at the centroid, 
			 * are used for 3D - by limiting coordinate distances. */
			for(int i = 0; i < biases.length; i++) {
				
				if(ignoreCentroidDirection) {
					i = 2;
				}
				if(!ignoreCentroidDirection && i == 2) {
					break;
				}
				
				bias_uZ = biases[i];
				
				this.pixelValues = projectVectors(numVectors, bias_uZ);
				
				findBoundaries(pixelValues);
				
				/* Add boundary coordinates to ArrayList (send method a new ArrayList) */
				newCoOrdinates = addNewCoOrdinatesToArrayList(new ArrayList<double[]>(), bias_uZ);
				
				switch (i) {
					case 0 : {	// 2D
						coOrdinateDistances(newCoOrdinates);
						
						/* Calculate mean coordinate distance (weighted by outliers). */
						meanCoOrdDistance = Centroid.getCentroid(this.coOrdDistances);
						sdCoOrdDistance = Math.sqrt(ShaftGuesser.variance(this.coOrdDistances));
						
						refineCoOrdinates(newCoOrdinates, i);
						coOrdinateDistances(newCoOrdinates);
						
						/* Recalculate mean coordinate distance (less weighted by outliers). */
						meanCoOrdDistance = Centroid.getCentroid(this.coOrdDistances);
						sdCoOrdDistance = Math.sqrt(ShaftGuesser.variance(this.coOrdDistances));
						
						/* Fill with new values from 2D; append 3D values. */
						coOrdinates = new ArrayList<double[]>();
						
						break;
					}
					case 1 : {	// 3D
						coOrdinateDistances(newCoOrdinates);
						refineCoOrdinates(newCoOrdinates, i);
						break;
					}
					
					case 2 : {	// Ignore centroid direction (for ellipsoids)
						
						refineCoOrdinates(newCoOrdinates, i);
						
						/* Fill with new values. */
						coOrdinates = new ArrayList<double[]>();
						break;
					}
				}
				
				/* Append the refined newCoOrdinates to coOrdinates */
				coOrdinates.addAll(newCoOrdinates);
			}
			
			/* Copy to array */
			this.coOrdArray = (double[][]) coOrdinates.toArray(new double[coOrdinates.size()][3]);
			
			if(showPoints) {
				/* Show refined points */
				annotateImage(imp, coOrdArray).show();
			}
			
			/* Give fitSphere something it can understand. */
			for(int n = 0; n < coOrdArray.length; n++) {
				coOrdArray[n][0] = coOrdArray[n][0] * vW;
				coOrdArray[n][1] = coOrdArray[n][1] * vH;
				coOrdArray[n][2] = coOrdArray[n][2] * vD;
//				IJ.log("Feeding fitSphere: " + coOrdArray[n][0] + ", " + coOrdArray[n][1] + ", " + coOrdArray[n][2] + "");
			}
			
			/* Fit sphere to points. This feeds fitSphere with units, like 
			 * SphereFitter. See ResultsTable for the current output format. */
			if(!fitEllipsoid) {
				try {
					dimensions = FitSphere.fitSphere(coOrdArray);
				} catch (IllegalArgumentException ia) {
					IJ.showMessage(ia.getMessage());
					return;
				} catch (RuntimeException re) {
					IJ.showMessage("Can't fit sphere to points.\n"
							+ "Rules may be too strict. Add more points and try again.");
					return;
				}
			}
			else {
				Object[] ellipsoid = new Object[5];
				try {
					ellipsoid = FitEllipsoid.yuryPetrov(coOrdArray); 
				} catch (IllegalArgumentException ia) {
					IJ.showMessage(ia.getMessage());
					return;
				} catch (RuntimeException re) {
					IJ.showMessage("Can't fit ellipsoid to points.\n"
							+ "Rules may be too strict. Add more points and try again.");
					return;
				}
				
				double[] ellipseCentroid = (double[]) ellipsoid[0];
				double[] ellipseRadii = (double[]) ellipsoid[1];
				dimensions[0] = ellipseCentroid[0];
				dimensions[1] = ellipseCentroid[1];
				dimensions[2] = ellipseCentroid[2];
				dimensions[3] = ellipseRadii[0];
				dimensions[4] = ellipseRadii[1];
				dimensions[5] = ellipseRadii[2];
			}
			
			
			/* Show results for each run */
			fillResultsTable("run_" + r + "", dimensions, new double[dimensions.length], fitEllipsoid).show("Results");
			
			allDimensions[r] = dimensions;
			
			/* Recursion */
			if(doRecursion) {
				/* Set initial point based on fitSphere's output */
				setInitialPoint(allDimensions[r][0] / vW, allDimensions[r][1] / vH, Math.floor(allDimensions[r][2] / vD));

				/* Check new sphere covers user's chosen initial point: 
				 * if not, reset to user's chosen initial point. Z is the most 
				 * likely axis to require adjustment. */
				if(r > 0) {
					if(Math.abs(allDimensions[r][0] - initialPointUser[0] * vW) > allDimensions[r][3]
					    || Math.abs(allDimensions[r][1] - initialPointUser[1] * vH) > allDimensions[r][3]
					    || Math.abs(allDimensions[r][2] - initialPointUser[2] * vD) > allDimensions[r][3]
					    || Math.abs(allDimensions[r][1] - initialPointUser[1] * vH) > allDimensions[r][3]	
						|| allDimensions[r][2] < 0 || allDimensions[r][2] > imp.getStackSize() * vD
						|| allDimensions[r][0] < 0 || allDimensions[r][0] > w
						|| allDimensions[r][1] < 0 || allDimensions[r][1] > h) {
						
						IJ.log("Initial point reset (run " + r + ").");
						
						setInitialPoint(initialPointUser[0], initialPointUser[1], Math.floor(initialPointUser[2]));
					}
				}
			}
		
		}	// end multiple runs
		
		/* Math for stats */
		double[][] allDimensionsT = new double[dimensions.length][runNum];
		sd = new double[dimensions.length];
		this.meanDimensions = new double[dimensions.length];
		allDimensionsT = transposeArray(allDimensions);
		for(int i = 0; i < dimensions.length; i++) {
			meanDimensions[i] = Centroid.getCentroid(allDimensionsT[i]);
			sd[i] = Math.sqrt(ShaftGuesser.variance(allDimensionsT[i]));
		}
	}
	
	
	public void findCentroids(ImagePlus imp, double min, double max) {
		
		Moments m = new Moments();
		/* Find centroid (of whole stack) */
		this.centroidW = m.getCentroid3D(imp, 1, imp.getImageStackSize(), min, max, 0, 1);
		/* Find centroid (of slice containing initialPoint) */
		this.centroidS = m.getCentroid3D(imp, iZ, iZ, min, max, 0, 1);
//		int midSlice = (int) imp.getImageStackSize() / 4;
//		if(iZ > midSlice) {
//			centroidS = m.getCentroid3D(imp, imp.getImageStackSize() - midSlice, iZ, min, max, 0, 1);
//		}
//		else {
//			centroidS = m.getCentroid3D(imp, iZ, midSlice, min, max, 0, 1);
//		}
		
		if (centroidW[0] < 0 || centroidS[0] < 0) {
			IJ.error("Empty Stack", "No voxels available for calculation."
					+ "\nCheck your ROI and threshold.");
			return;
		}
		
		/* Get unit vectors of initialPoint-Centroid direction */
		double iCDist = Trig.distance3D(iX * vW, iY * vH, iZ * vD, centroidW[0], centroidW[1], centroidW[2]);
		this.uCW = new double[3];
		uCW[0] = (iX * vW - centroidW[0]) / iCDist;
		uCW[1] = (iY * vH - centroidW[1]) / iCDist;
		uCW[2] = (iZ * vD - centroidW[2]) / iCDist;
	
		iCDist = Trig.distance3D(iX * vW, iY * vH, iZ * vD, centroidS[0], centroidS[1], centroidS[2]);
		this.uCS = new double[3];
		uCS[0] = (iX * vW - centroidS[0]) / iCDist;
		uCS[1] = (iY * vH - centroidS[1]) / iCDist;
		uCS[2] = (iZ * vD - centroidS[2]) / iCDist;

		/* To discard a 3D hemisphere of vectors around this vector, 
		 * discard if acos(uC (dot) uV) > 90 degrees (pi/2 radians). */
	}
	
	
	/**
	 * <p>Create random 3D unit vectors, with the Z direction weighted 
	 * by multiplying by bias_uZ (so for single slice, set bias_uZ = 0).</p>
	 * 
	 * <p>Fills both pixelValues and magnitudes.</p>
	 * 
	 * @param numVectors number of vectors to create
	 * @param bias_uZ weighting given to uZ
	 */
	private ArrayList<Float>[] projectVectors(int numVectors, double bias_uZ) {
		
		/* Create array of random 3D unit vectors */
		this.unitVectors = Vectors.random3D(numVectors);
		
		/* Array of ArrayLists */
		this.pixelValues = new ArrayList[unitVectors.length];
//		for(int p = 0; p < pixelValues.length; p++) {
//			pixelValues[p] = new ArrayList<Float>();
//		}
		this.magnitudes = new double[unitVectors.length];
		
		/* Cycle through all vectors */
		for(int i = 0; i < unitVectors.length; i++) {
			
			uX = unitVectors[i][0];
			uY = unitVectors[i][1];
			uZ = unitVectors[i][2] * bias_uZ;
			
			pixelValues[i] = getPixelsAlongVector(uX, uY, uZ);
			
			this.magnitudes[i] = Trig.distance3D(uX * vW, uY * vH, uZ * vD);
//			this.magnitudes[i] = Math.sqrt(Math.pow(uX * vW, 2) + Math.pow(uY * vH, 2) + Math.pow(uZ * vD, 2));
		}
		
		return pixelValues;
	}
	
	
	/**
	 * <p>Find pixel values at all pixels along a vector, from the initial point 
	 * determined by the user, to the border of the image.</p>
	 * 
	 * <p>Uses getPixelValue(int x, int y).</p>
	 * 
	 * @param uX unit vector in X direction
	 * @param uY unit vector in Y direction
	 * @param uZ unit vector in Z direction
	 * @return ArrayList<Float> vectorPixelValues
	 */
	private ArrayList<Float> getPixelsAlongVector(double uX, double uY, double uZ) {
		
		/* Use ArrayList as we don't know how long each line (j) will be.
		 * Unfortunately, ArrayList can only be filled with Objects. */
		this.vectorPixelValues = new ArrayList<Float>();
		
		int j = 0;
		while(j > -1) {
			
			this.nX = (int) Math.floor(iX + (j * uX));		// for 'unit' resolution, divide by vW, etc.
			this.nY = (int) Math.floor(iY + (j * uY));
			this.nZ = (int) Math.round(iZ + (j * uZ));
			
			/* Break out if outside the image */
			if(nX < 0 || nX > w || nY < 0 || nY > h || nZ < 1 || nZ > d) {
				j = -2;
				break;
			}
			
			Float nowValue = new Float(sliceProcessors[nZ].getPixelValue(nX, nY));
			
			vectorPixelValues.add(nowValue);
			j++;
		}
		
		return vectorPixelValues;
	}
	
	
	/**
	 * <p>Find the bone boundaries along each vector, based on the median 
	 * pixel values encountered by each individual vector.</p>
	 * 
	 * <p>If final 'boundary points' are found to be incorrectly lying on, 
	 * e.g. large scale image artefacts, adjustedLimit here probably needs 
	 * further adjusting.</p>
	 */
	private void findBoundaries(ArrayList<Float>[] pixelValues) {
		
		/* Some analysis possible on each vector */
		this.boundarySteps = new int[pixelValues.length];
		this.medianValues = new double[pixelValues.length];
		this.boundaryDistances = new double[pixelValues.length];
		
		/* Ragged array of (double) pixel values for every vector (i) */
		this.pxVals = new double[pixelValues.length][];
		for(int i = 0; i < pxVals.length; i++) {
			pxVals[i] = new double[pixelValues[i].size()];
			
			for(int j = 0; j < pixelValues[i].size(); j++) {
				pxVals[i][j] = (double) (Float) pixelValues[i].toArray()[j];
			}
			
			/* Guessing limit and adjustment based on the typical plots seen */
			medianValues[i] = ShaftGuesser.median(pxVals[i]);
			if(medianValues[i] < 0) {
				adjustedLimit = medianValues[i] * 0.5;
			}
			else {
				adjustedLimit = medianValues[i] * 1.25;
			}
			
			boundarySteps[i] = boundaryLimiter(pxVals[i], adjustedLimit, (pxVals[i].length - 1), false);
			boundaryDistances[i] = boundarySteps[i] * magnitudes[i];
		}
	}
	
	
	/**
	 * <p>Along a vector containing pixel values beginning within the suspected 
	 * sphere (femoral head), finds where the median pixel value (of this vector) 
	 * is crossed. Actually starts from the *end* of the vector, i.e. outside 
	 * the sphere. This is taken to be the most significant bone boundary which the 
	 * vector crosses, i.e, the outer surface of the bone.</p>
	 * 
	 * <p>In case narrow surfaces (e.g. scanning artefacts) are found accidentally, 
	 * the code looks 2%, 5% and 10% further along the vector: if this was a narrow 
	 * peak, it is discounted.</p>
	 * 
	 * <p>Modified ShaftGuesser.shaftLimiter</p>
	 * 
	 * @param pxlValues
	 * @param limit
	 * @param startPoint begin here: usually (pxlValues.length - 1)
	 * @param isMin Specify whether the limit is a minimum (true) or a maximum (false).
	 * @return boundaryPosition (int): the number of *steps* (pixels) along the vector 
	 * from initialPoint.
	 */
	private int boundaryLimiter(double[] pxlValues, double limit, int startPoint, boolean isMin) {
		
		int boundaryPosition = 0;
		
//		Arrays.sort(sampleAheadPc);
		int pc1 = (int) Math.floor((pxlValues.length) * sampleAheadPc[0] / 100);
		int pc2 = (int) Math.floor((pxlValues.length) * sampleAheadPc[1] / 100);
		int pc3 = (int) Math.floor((pxlValues.length) * sampleAheadPc[2] / 100);
		
		/* Start from the end of the vector */
		for(int i = startPoint; i >= 0; i--) {
			if(isMin && pxlValues[i] < limit) {
				boundaryPosition = i;
				/* Look ahead */
				if(i - pc3 > 0 && ((pxlValues[i - pc1] > limit && pxlValues[i - pc2] > limit) || (pxlValues[i - pc3] > limit))) {
					boundaryPosition = boundaryLimiter(pxlValues, limit, (i - pc2), isMin);
				}
				break;
			}
			else if(!isMin && pxlValues[i] > limit) {
				boundaryPosition = i;
				/* Look ahead */
				if(i - pc3 > 0 && ((pxlValues[i - pc1] < limit && pxlValues[i - pc2] < limit) || (pxlValues[i - pc3] < limit))) {
					boundaryPosition = boundaryLimiter(pxlValues, limit, (i - pc2), isMin);
				}
				else { break; }
			}
		}
		return boundaryPosition;
	}
	
	
	/**
	 * Get complete ArrayList of boundary coordinates. Fills a new ArrayList; 
	 * append it to an existing one outside the method, if necessary, using 
	 * arrayListName.addAll(addCoOrdinatesToArrayList(new ArrayList<double[]>()))
	 * 
	 * @param newCoOrdinates ArrayList
	 * @return ArrayList<double[]> newCoOrdinates with new coordinates appended.
	 */
	private ArrayList<double[]> addNewCoOrdinatesToArrayList(ArrayList<double[]> newCoOrdinates, double bias_uZ) {
		
		/* Fill, or add to, ArrayList of boundary coordinates */
		for(int i = 0; i < this.boundarySteps.length; i++) {
			
			this.coOrds = new double[3];
			coOrds[0] = (int) Math.floor(iX + (boundarySteps[i] * unitVectors[i][0]));
			coOrds[1] = (int) Math.floor(iY + (boundarySteps[i] * unitVectors[i][1]));
			coOrds[2] = (int) Math.floor(iZ + (boundarySteps[i] * unitVectors[i][2]) * bias_uZ);
			
			newCoOrdinates.add(coOrds);
		}
		return newCoOrdinates;
	}
	
	/**
	 * <p>Calculate coordinate distances in 2D (will be slightly lower than boundaryDistances 
	 * as coOrdinates are Math.floored). <i>Must</i> be re-done every time coOrdinates 
	 * is refined, as this sets up an array.</p>
	 * 
	 * @param coOrdinates
	 * @return
	 */
	private double[] coOrdinateDistances(ArrayList<double[]> coOrdinates) {
		
		this.coOrdDistances = new double[coOrdinates.size()];
		
		Iterator<double[]> itA = coOrdinates.iterator();
		int i = 0;
		while (itA.hasNext()) {
			double[] a = itA.next();
//			coOrdDistances[i] = Trig.distance3D(a[0] * vW, a[1] * vH, a[2] * vD, initialPoint[0], initialPoint[1], initialPoint[2]);
			coOrdDistances[i] = Math.sqrt(Math.pow((a[0] - initialPoint[0]) * vW, 2) 
					+ Math.pow((a[1] - initialPoint[1]) * vH, 2) 
					+ Math.pow((a[2] - initialPoint[2]) * vD, 2));
			i++;
		}
		
		return coOrdDistances;
	}
	
	/**
	 * <p>Coordinates refined based on means and standard deviations: 
	 * outliers removed.</p>
	 * 
	 * <p>
	 * <b>case 0: </b>2D: large outliers removed
	 * </br>
	 * <b>case 1: </b>3D: means from refined 2D results are used to refine 3D results
	 * </p>
	 * 
	 * @param caseNum sets which math to use on the data
	 * @param coOrdinates ArrayList containing all coordinates, x,y,z, to be refined
	 */
	private void refineCoOrdinates(ArrayList<double[]> coOrdinates, int caseNum) {
		
		/* Do different things based on which run this is */
		switch (caseNum) {
		
		/* Refine based on distance to the points */
		
		/* 2D */
		case 0 : {
			
			/* Remove large outliers */
			Iterator<double[]> itB = coOrdinates.iterator();
			int j = 0;
			while (itB.hasNext()) {
				double[] b = itB.next();
				if(coOrdDistances[j] > meanCoOrdDistance + (sdCoOrdDistance * -sd2DMult)) {
					itB.remove();
				}
				j++;
			}
			
			break;
		}
		
		/* 3D relies upon 2D */
		case 1 : {
			
			/* Refine based on 2D distances */
			Iterator<double[]> itC = coOrdinates.iterator();
			int i = 0;
			while (itC.hasNext()) {
				double[] c = itC.next();
				if(coOrdDistances[i] > meanCoOrdDistance + (sdCoOrdDistance * sd3DMult)
						 || coOrdDistances[i] < meanCoOrdDistance + (sdCoOrdDistance * -sd3DMult)) {
					itC.remove();
				}
				i++;
			}
			
			break;
		}
		
		/* Or, refine based upon vector directions being not within a hemisphere (90 
		 * degrees in all directions) of the vector connecting the initialPoint with 
		 * the centroid (of this half of the bone, so it's not too exaggerated). */
		
		case 2 : {
			
			/* Haven't removed any generated coOrdinates at this point (still have 
			 * numVectors), so must still have numVectors unit vectors. */
			Iterator<double[]> itD = coOrdinates.iterator();
			int i = 0;
			while (itD.hasNext()) {
				double[] d = itD.next();
				if(Math.acos(uCW[0] * unitVectors[i][0]
				            + uCW[1] * unitVectors[i][1] 
				            + uCW[2] * unitVectors[i][2]) > Math.PI * excludeWholeC
				            || this.useSliceCentroid 
				            && Math.acos(uCS[0] * unitVectors[i][0]
				            + uCS[1] * unitVectors[i][1] 
				            + uCS[2] * unitVectors[i][2]) > Math.PI * excludeSliceC) {
					itD.remove();
					
//					double toACos = uCW[0] * unitVectors[i][0] + uCW[1] * unitVectors[i][1] + uCW[2] * unitVectors[i][2];
//					IJ.log("I am removing a coordinate!");
//					IJ.log("MathdotPI / 2: " + (Math.PI / 2) + "");
//					IJ.log("uC; uV: " + uCX + ", " + uCY + ", " + uCZ + "; " + unitVectors[i][0] + ", " + unitVectors[i][1] + ", " + unitVectors[i][2] + "");
//					IJ.log("Pre-Math.acos: " + toACos + "");
//					IJ.log("Math.acos: " + Math.acos(uCX * unitVectors[i][0] + uCY * unitVectors[i][1] + uCZ * unitVectors[i][2]) + "");
				}
				
				i++;
			}
			
			break;
		}
		
		
		default: break;
		}
	}
	
	
	/**
	 * Provided a 2D inputArray[a][b] with constant number of columns per row,
	 * i.e. not a ragged array, transposes this array so that row and column 
	 * values are switched.
	 * 
	 * @param inputArray
	 * @return transposedArray[b][a]
	 */
	public double[][] transposeArray(double[][] inputArray) {
		
		double[][] transposedArray = new double[inputArray[0].length][inputArray.length];
		
		/* Cycle through inputArray's rows, then columns */
		for(int i = 0; i < inputArray.length; i++) {
			for(int j = 0; j < inputArray[0].length; j++) {
				transposedArray[j][i] = inputArray[i][j];
			}
		}
		
		return transposedArray;
	}
	
	
	/**
	 * Draw points used to calculate the sphere on a copy of the original image.
	 * Modified from SliceGeometry.
	 * 
	 * @param imp
	 * @param coOrdinateArray double[n][3] containing n (x, y, z) coordinates.
	 * @return ImagePlus with points drawn
	 */
	private ImagePlus annotateImage(ImagePlus imp, double[][] coOrdinateArray) {
		ImageStack stack = imp.getImageStack();
		int w = stack.getWidth();
		int h = stack.getHeight();
		ImageStack annStack = new ImageStack(w, h);
		for (int s = 1; s < imp.getStackSize() + 1; s++) {
			ImageProcessor annIP = stack.getProcessor(s).duplicate();
			annIP.setColor(Color.white);
			
			for(int i = 0; i < coOrdinateArray.length; i++) {
				if(s == (int) coOrdinateArray[i][2]) {
					annIP.drawDot((int) Math.floor(coOrdinateArray[i][0]), 
							(int) Math.floor(coOrdinateArray[i][1]));
					annIP.drawOval((int) Math.floor(coOrdinateArray[i][0]), 
							(int) Math.floor(coOrdinateArray[i][1]), 4, 4);
				}
			}
			
			annStack.addSlice(stack.getSliceLabel(s), annIP);
		}
		ImagePlus ann = new ImagePlus("Annotated_" + imp.getTitle(), annStack);
		ann.setCalibration(imp.getCalibration());
		return ann;
	}
	
	
	/**
	 * Fill a ResultsTable with values.
	 * 
	 * @param label
	 * @param dimensions
	 * @param deviations
	 * @param fitEllipsoid
	 * @return
	 */
	private ResultsTable fillResultsTable(String label, double[] dimensions, double[] deviations, boolean fitEllipsoid) {
		
		ResultsTable rt = ResultsTable.getResultsTable();
		rt.incrementCounter();
		rt.addLabel("Type", label);
		rt.addValue("Points_used", coOrdinates.size());
		rt.addValue("X_centroid_(" + this.units + ")", dimensions[0]);
		rt.addValue("Y_centroid_(" + this.units + ")", dimensions[1]);
		rt.addValue("Z_centroid_(approx._slice)", dimensions[2] / vD);
		rt.addValue("Radius_1_(" + this.units + ")", dimensions[3]);
		if(fitEllipsoid) {
			rt.addValue("Radius_2_(" + this.units + ")", dimensions[4]);
			rt.addValue("Radius_3_(" + this.units + ")", dimensions[5]);
		}
		rt.addValue("sdX_(" + this.units + ")", deviations[0]);
		rt.addValue("sdY_(" + this.units + ")", deviations[1]);
		rt.addValue("sdZ_(" + this.units + ")", deviations[2]);
		rt.addValue("sdR_(" + this.units + ")", deviations[3]);
		if(fitEllipsoid) {
			rt.addValue("sdZ_(" + this.units + ")", deviations[2]);
			rt.addValue("sdR_(" + this.units + ")", deviations[3]);
		}
		return rt;
	}
	
	
	/**
	 * Draw lines showing centre and radius of the sphere on a copy of the 
	 * original image.
	 * 
	 * @param imp
	 * @param coOrdinateArray
	 * @return
	 */
	public ImagePlus annotateCentre(ImagePlus imp, double[] dimensions) {
		ImageStack stack = imp.getImageStack();
		int w = stack.getWidth();
		int h = stack.getHeight();
		ImageStack annStack = new ImageStack(w, h);
		for (int s = 1; s < imp.getStackSize() + 1; s++) {
			ImageProcessor annIP = stack.getProcessor(s).duplicate();
			annIP.setColor(Color.white);
			
			if(s == (int) Math.round(dimensions[2] / vD)) {
				
				int x1 = (int) Math.floor((dimensions[0] / vW) - dimensions[3] / vW);
				int y1 = (int) Math.floor(dimensions[1] / vH);
				int x2 = (int) Math.floor((dimensions[0] / vW) + dimensions[3] / vW);
				int y2 = (int) Math.floor(dimensions[1] / vH);
				annIP.drawLine(x1, y1, x2, y2);

				x1 = (int) Math.floor(dimensions[0] / vW);
				y1 = (int) Math.floor((dimensions[1] / vH) + dimensions[3] / vH);
				x2 = (int) Math.floor(dimensions[0] / vW);
				y2 = (int) Math.floor((dimensions[1] / vH) - dimensions[3] / vH);
				annIP.drawLine(x1, y1, x2, y2);
			}
			
			annStack.addSlice(stack.getSliceLabel(s), annIP);
		}
		ImagePlus ann = new ImagePlus("Annotated_" + imp.getTitle(), annStack);
		ann.setCalibration(imp.getCalibration());
		return ann;
	}
	
	
	/* Copied from Neck Shaft Angle: only .get(#)s changed */
	public boolean dialogItemChanged(GenericDialog gd, AWTEvent e) {
		Vector<?> checkboxes = gd.getCheckboxes();
		Vector<?> nFields = gd.getNumericFields();
		Checkbox box0 = (Checkbox) checkboxes.get(3);
		boolean isHUCalibrated = box0.getState();
		double min = 0;
		double max = 0;
		TextField minT = (TextField) nFields.get(6);
		TextField maxT = (TextField) nFields.get(7);
		try{
			min = Double.parseDouble(minT.getText());
			max = Double.parseDouble(maxT.getText());
		} catch (Exception ex){
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
