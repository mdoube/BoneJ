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
public class SphereEdgeGuesser implements PlugIn, MouseListener {
	
	private ImageProcessor[] sliceProcessors = null;
	private ImageCanvas canvas;
	private Calibration cal;
	
	/** Only measure a single slice of a multi-slice image */
	private boolean singleSlice = false;
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
	
	/** Linear unit of measure */
	private String units;
	/** Pixel dimensions in 'unit's */
	private double vH, vW, vD;
	/** Unit vector sizes, -1 < value < 1 */
	private double uX, uY, uZ;
	/** Unit vectors of the vector joining initialPoint - Centroid */
	private double uCX, uCY, uCZ;
	/** Standard deviation values of the centre point and radius */
	private double sdX, sdY, sdZ, sdR;
	
	/** Midslice for Moments to get centroid of the near half the bone */
	private int midSlice;
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
	private ArrayList[] pixelValues;
	
	/** Centroid of the bone. W: whole stack; S: slice. */
	private double[] centroidW, centroidS;
	/** Holds (x, y, z) centre and radius, in same units as those given to fitSphere */
	private double[] sphereDim = new double[4];
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
	/** So far, these haven't needed adjustment, but they could be, for large 
	 * variations in trabeculae sizes. Sample (get pixel values) ahead along 
	 * vectors these percentages of the vector. 
	 * To jump ahead, the pixel values at the two lowest percentages must be beneath 
	 * the limit - or just the largest value. */
	private double[] sampleAheadPc = {2,5,10};
	
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
		
		/* Set up array of ImageProcessors for reference later */
		ImageStack stack = imp.getStack();
		this.sliceProcessors = new ImageProcessor[imp.getStackSize() + 1];
		for(int s = 1; s < sliceProcessors.length; s++) {
			sliceProcessors[s] = stack.getProcessor(s);
		}
		
		w = imp.getWidth();
		h = imp.getHeight();
		d = imp.getStackSize();
		
		ImageWindow win = imp.getWindow();
		this.canvas = win.getCanvas();
		cal = imp.getCalibration();
		
		vW = cal.pixelWidth;
		vH = cal.pixelHeight;
		vD = cal.pixelDepth;
		units = cal.getUnits();
		
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
		gd.addNumericField("Limit 2D: mean - ", sd2DMult, 2, 5, "standard deviations");
		gd.addNumericField("Limit 3D: mean +/-", sd3DMult, 2 , 5, "standard deviations");
//		gd.addNumericField("Sample ahead along vector by", sampleAheadPc[0], 0, 4, "% ");
//		gd.addNumericField("and", sampleAheadPc[1], 0, 4, "% ");
//		gd.addNumericField("or", sampleAheadPc[2], 0, 4, "% ");
		gd.addCheckbox("Use centroid method", ignoreCentroidDirection);
		gd.addCheckbox("HU Calibrated", ImageCheck.huCalibrated(imp));
		gd.addNumericField("Bone Min:", min, 1, 6, pixUnits + " ");
		gd.addNumericField("Bone Max:", max, 1, 6, pixUnits + " ");
		gd.addCheckbox("Use recursion", doRecursion);
		gd.addCheckbox("Single slice (don't fit sphere)", singleSlice);
		gd.addCheckbox("Show points", showPoints);
		gd.addCheckbox("Plot a vector's profile", doPlot);
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
		boolean isHUCalibrated = gd.getNextBoolean();
		if (isHUCalibrated) {
			min = cal.getRawValue(min);
			max = cal.getRawValue(max);
		}
		min = gd.getNextNumber();
		max = gd.getNextNumber();
		this.doRecursion = gd.getNextBoolean();
		this.singleSlice = gd.getNextBoolean();
		this.showPoints = gd.getNextBoolean();
		this.doPlot = gd.getNextBoolean();
		
		if(singleSlice) {
			bias_uZ = 0;
		}
		
		/* User clicks to set initial point */
		
		// remove stale MouseListeners
		MouseListener[] l = this.canvas.getMouseListeners();
		for (int n = 0; n < l.length; n++) {
			this.canvas.removeMouseListener(l[n]);
		}
		// add a new MouseListener
		this.canvas.addMouseListener(this);
		new WaitForUserDialog("Click (roughly) in the centre of\n" 
									+ "the femoral head, then hit \'OK\'").show();
		this.canvas.removeMouseListener(this);
		
		/* Confirm initial point and measurements */
//		GenericDialog gd1 = new GenericDialog("Info");
//		gd1.addMessage("Initial point (x,y,z): (" + iX + ", " + iY + ", " + iZ + ")");
//		gd1.addMessage("Pixel sizes (vW,vH,vD): (" + vW + ", " + vH + ", " + vD + ")");
//		gd1.addMessage("imp.getCurrentSlice: " + imp.getCurrentSlice() + "");
//		gd1.showDialog();
//		IJ.log("Pixel sizes (vW,vH,vD): " + vW + ", " + vH + ", " + vD + "");
		IJ.log("Initial point: " + iX + ", " + iY + ", " + iZ + "");
		
		if(ignoreCentroidDirection) {
			/* Find centroid (of whole stack) */
			Moments m = new Moments();
			this.midSlice = imp.getImageStackSize();
			if(iZ < midSlice) {
				centroidW = m.getCentroid3D(imp, 1, midSlice, min, max, 0, 1);
			}
			else {
				centroidW = m.getCentroid3D(imp, midSlice, imp.getImageStackSize(), min, max, 0, 1);
			}
			if (centroidW[0] < 0) {
				IJ.error("Empty Stack", "No voxels available for calculation."
						+ "\nCheck your ROI and threshold.");
				return;
			}
			
			/* Get unit vectors of initialPoint-Centroid direction */
			double iCDist = Trig.distance3D(iX * vW, iY * vH, iZ * vD, centroidW[0], centroidW[1], centroidW[2]);
			uCX = (iX * vW - centroidW[0]) / iCDist;
			uCY = (iY * vH - centroidW[1]) / iCDist;
			uCZ = (iZ * vD - centroidW[2]) / iCDist;
			
			/* Essentially, to discard a 3D hemisphere of vectors around this vector, 
			 * just discard if acos(uC (dot) uV) > 90 degrees (pi/2 radians). */
			
			IJ.log("Centroid (whole stack): " + centroidW[0] + ", " + centroidW[1] + ", " + centroidW[2] + "");
			IJ.log("Centroid (initial slice): " + centroidW[0] + ", " + centroidW[1] + ", " + centroidW[2] + "");
			IJ.log("unit Vectors: " + uCX + ", " + uCY + ", " + uCZ + "");
			IJ.log("iCDist: " + iCDist + "");
		}
		
		
		
		this.biases = new double[3];
		biases[0] = 0;			// Run 1: vectors in current slice only (2D)
		biases[1] = 1;			// Run 2: 3D vectors
		biases[2] = 1;			// Or: ignore vectors towards the centroid
		
		double[][] sphereDims = new double[runNum][4];
		
		
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
				
				projectVectors(numVectors, bias_uZ);
				
				findBoundaries();
				
				/* Add boundary coordinates to ArrayList (send method a new ArrayList) */
				newCoOrdinates = addNewCoOrdinatesToArrayList(new ArrayList<double[]>(), bias_uZ);
				
				switch (i) {
					case 0 : {
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
					case 1 : {
						coOrdinateDistances(newCoOrdinates);
						refineCoOrdinates(newCoOrdinates, i);
						break;
					}
					
					case 2 : {
						
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
				/* Show all detected boundaries */
//				annotateImage(imp, coOrdArray).show();
				/* Show points used by FitSphere */
				annotateImage(imp, coOrdArray).show();
			}
			
			/* Give fitSphere something it can understand. */
			for(int n = 0; n < coOrdArray.length; n++) {
				coOrdArray[n][0] = coOrdArray[n][0] * vW;
				coOrdArray[n][1] = coOrdArray[n][1] * vH;
				coOrdArray[n][2] = coOrdArray[n][2] * vD;
//				IJ.log("Feeding fitSphere: " + coOrdArray[n][0] + ", " + coOrdArray[n][1] + ", " + coOrdArray[n][2] + "");
			}
			
			/* Fit sphere to points.
			 * Nb. This feeds fitSphere with units, like SphereFitter. 
			 * See ResultsTable for the current output format. */
			if(!singleSlice) {
				try {
					sphereDim = FitSphere.fitSphere(coOrdArray);
				} catch (IllegalArgumentException ia) {
					IJ.showMessage(ia.getMessage());
					return;
				} catch (RuntimeException re) {
					IJ.showMessage("Can't fit sphere to points.\n"
								+ "Rules may be too strict. Try again.");
					return;
				}
			}
			
			/* Show results for each run */
			ResultsTable rt = ResultsTable.getResultsTable();
			rt.incrementCounter();
			rt.addLabel("Type", "run_" + r + "");
			rt.addValue("Points_used", coOrdinates.size());
			rt.addValue("X_centroid_(" + units + ")", sphereDim[0]);
			rt.addValue("Y_centroid_(" + units + ")", sphereDim[1]);
			rt.addValue("Z_centroid_(approx._slice)", sphereDim[2] / vD);
			rt.addValue("Radius_(" + units + ")", sphereDim[3]);
			rt.addValue("sdX_(" + units + ")", sdX);
			rt.addValue("sdY_(" + units + ")", sdY);
			rt.addValue("sdZ_(" + units + ")", sdZ);
			rt.addValue("sdR_(" + units + ")", sdR);
			rt.show("Results");
			
			sphereDims[r] = sphereDim;
			
			/* Recursion */
			if(doRecursion) {
				/* Set initial point based on fitSphere's output */
				setInitialPoint(sphereDims[r][0] / vW, sphereDims[r][1] / vH, Math.floor(sphereDims[r][2] / vD));

				/* Check new sphere covers user's chosen initial point: 
				 * if not, reset to user's chosen initial point. Z is the most 
				 * likely axis to require adjustment. */
				if(r > 0) {
					if(sphereDims[r][0] - sphereDims[r][3] > initialPointUser[0] * vW
						|| sphereDims[r][0] + sphereDims[r][3] < initialPointUser[0] * vW
						|| sphereDims[r][1] - sphereDims[r][3] > initialPointUser[1] * vH
						|| sphereDims[r][1] + sphereDims[r][3] < initialPointUser[1] * vH
						|| sphereDims[r][2] - sphereDims[r][3] > initialPointUser[2] * vD
						|| sphereDims[r][2] + sphereDims[r][3] < initialPointUser[2] * vD) {
						
						IJ.log("Initial point reset (run " + r + ").");
						
						setInitialPoint(initialPointUser[0], initialPointUser[1], Math.floor(initialPointUser[2]));
					}
				}
			}
		
		}
		
		double[][] sphereDimsT = new double[4][runNum];
		sphereDimsT = transposeArray(sphereDims);
		sphereDim[0] = Centroid.getCentroid(sphereDimsT[0]);
		sphereDim[1] = Centroid.getCentroid(sphereDimsT[1]);
		sphereDim[2] = Centroid.getCentroid(sphereDimsT[2]);
		sphereDim[3] = Centroid.getCentroid(sphereDimsT[3]);
		
		sdX = Math.sqrt(ShaftGuesser.variance(sphereDimsT[0]));
		sdY = Math.sqrt(ShaftGuesser.variance(sphereDimsT[1]));
		sdZ = Math.sqrt(ShaftGuesser.variance(sphereDimsT[2]));
		sdR = Math.sqrt(ShaftGuesser.variance(sphereDimsT[3]));
		
		/* Show statistical results */
		ResultsTable rt2 = ResultsTable.getResultsTable();
		rt2.incrementCounter();
		rt2.addLabel("Type", "stats");
		rt2.addValue("Points_used", coOrdinates.size());
		rt2.addValue("X_centroid_(" + units + ")", sphereDim[0]);
		rt2.addValue("Y_centroid_(" + units + ")", sphereDim[1]);
		rt2.addValue("Z_centroid_(approx._slice)", sphereDim[2] / vD);
		rt2.addValue("Radius_(" + units + ")", sphereDim[3]);
		rt2.addValue("sdX_(" + units + ")", sdX);
		rt2.addValue("sdY_(" + units + ")", sdY);
		rt2.addValue("sdZ_(" + units + ")", sdZ);
		rt2.addValue("sdR_(" + units + ")", sdR);
		rt2.show("Results");
		
		annotateCentre(imp).show();
		
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
	public double[] getSphereDim() {
		return this.sphereDim;
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
	 * <p>Create random 3D unit vectors, with the Z direction weighted 
	 * by multiplying by bias_uZ (so for single slice, set bias_uZ = 0).</p>
	 * 
	 * <p>Fills both pixelValues and magnitudes.</p>
	 * 
	 * @param numVectors number of vectors to create
	 * @param bias_uZ weighting given to uZ
	 */
	private void projectVectors(int numVectors, double bias_uZ) {
		
		/* Create array of random 3D unit vectors */
		this.unitVectors = Vectors.random3D(numVectors);
		
		/* Array of ArrayLists */
		this.pixelValues = new ArrayList[unitVectors.length];
		
		this.magnitudes = new double[unitVectors.length];
		
		/* Cycle through all vectors */
		for(int i = 0; i < unitVectors.length; i++) {
			
			uX = unitVectors[i][0];
			uY = unitVectors[i][1];
			uZ = unitVectors[i][2] * bias_uZ;
			
			this.pixelValues[i] = getPixelsAlongVector(uX, uY, uZ);
			
			this.magnitudes[i] = Trig.distance3D(uX * vW, uY * vH, uZ * vD);
//			this.magnitudes[i] = Math.sqrt(Math.pow(uX * vW, 2) + Math.pow(uY * vH, 2) + Math.pow(uZ * vD, 2));
		}
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
			
			nX = (int) Math.floor(iX + (j * uX));		// for 'unit' resolution, divide by vW, etc.
			nY = (int) Math.floor(iY + (j * uY));
			nZ = (int) Math.round(iZ + (j * uZ));
			
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
	private void findBoundaries() {
		
		/* Some analysis possible on each vector */
		this.boundarySteps = new int[pixelValues.length];
		this.medianValues = new double[pixelValues.length];
		this.boundaryDistances = new double[pixelValues.length];
		
		/* Ragged array of (double) pixel values for every vector */
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
		
		/* Old method */
		int boundaryPosition = 0;
		int twoPc = (int) Math.floor(pxlValues.length * 0.02);
		int fivePc = (int) Math.floor(pxlValues.length * 0.05);
		int tenPc = (int) Math.floor(pxlValues.length * 0.1);
		
		/* New method can use user input */
		Arrays.sort(sampleAheadPc);
		int pcOne = (int) Math.floor((pxlValues.length) * sampleAheadPc[0] / 100);
		int pcTwo = (int) Math.floor((pxlValues.length) * sampleAheadPc[1] / 100);
		int pcThree = (int) Math.floor((pxlValues.length) * sampleAheadPc[2] / 100);
		
		/* Start from the end of the vector */
		for(int i = startPoint; i >= 0; i--) {
			if(isMin && pxlValues[i] < limit) {
				boundaryPosition = i;
				/* Look ahead */
				if(i - pcThree > 0 && (pxlValues[i - pcThree] > limit || pxlValues[i - pcTwo] > limit) && pxlValues[i - pcOne] > limit) {
					boundaryPosition = boundaryLimiter(pxlValues, limit, (i - pcOne), isMin);
				}
				else { break; }
			}
			else if(!isMin && pxlValues[i] > limit) {
				boundaryPosition = i;
				/* Look ahead */
				if(i - pcThree > 0 && (pxlValues[i - pcThree] < limit || pxlValues[i - pcTwo] < limit) &&  pxlValues[i - pcOne] < limit) {
					boundaryPosition = boundaryLimiter(pxlValues, limit, (i - pcOne), isMin);
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
			
			if(singleSlice) {
				coOrds[2] = (int) Math.floor(iZ + (boundarySteps[i] * bias_uZ));
			}
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
				if(Math.acos(uCX * unitVectors[i][0]
				            + uCY * unitVectors[i][1] 
				            + uCZ * unitVectors[i][2]) > Math.PI / 2) {
					itD.remove();
					
					double toACos = uCX * unitVectors[i][0] + uCY * unitVectors[i][1] + uCZ * unitVectors[i][2];
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
	 * Draw lines showing centre and radius of the sphere on a copy of the 
	 * original image.
	 * 
	 * @param imp
	 * @param coOrdinateArray
	 * @return
	 */
	private ImagePlus annotateCentre(ImagePlus imp) {
		ImageStack stack = imp.getImageStack();
		int w = stack.getWidth();
		int h = stack.getHeight();
		ImageStack annStack = new ImageStack(w, h);
		for (int s = 1; s < imp.getStackSize() + 1; s++) {
			ImageProcessor annIP = stack.getProcessor(s).duplicate();
			annIP.setColor(Color.white);
			
			if(s == (int) Math.round(sphereDim[2] / vD)) {
				
				int x1 = (int) Math.floor((sphereDim[0] / vW) - sphereDim[3] / vW);
				int y1 = (int) Math.floor(sphereDim[1] / vH);
				int x2 = (int) Math.floor((sphereDim[0] / vW) + sphereDim[3] / vW);
				int y2 = (int) Math.floor(sphereDim[1] / vH);
				annIP.drawLine(x1, y1, x2, y2);

				x1 = (int) Math.floor(sphereDim[0] / vW);
				y1 = (int) Math.floor((sphereDim[1] / vH) + sphereDim[3] / vH);
				x2 = (int) Math.floor(sphereDim[0] / vW);
				y2 = (int) Math.floor((sphereDim[1] / vH) - sphereDim[3] / vH);
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
		Checkbox box0 = (Checkbox) checkboxes.get(1);
		boolean isHUCalibrated = box0.getState();
		double min = 0;
		double max = 0;
		TextField minT = (TextField) nFields.get(4);
		TextField maxT = (TextField) nFields.get(5);
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
