package org.doube.bonej;

import java.awt.Color;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Random;

import org.doube.geometry.Centroid;
import org.doube.geometry.FitCircle;
import org.doube.geometry.FitSphere;
import org.doube.geometry.Trig;
import org.doube.geometry.Vectors;
import org.doube.util.ImageCheck;
import org.doube.util.RoiMan;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.gui.ImageCanvas;
import ij.gui.ImageWindow;
import ij.gui.Plot;
import ij.gui.ProfilePlot;
import ij.gui.Roi;
import ij.gui.WaitForUserDialog;
import ij.measure.Calibration;
import ij.measure.ResultsTable;
import ij.plugin.PlugIn;
import ij.plugin.filter.Profiler;
import ij.plugin.frame.RoiManager;
import ij.process.BinaryProcessor;
import ij.process.ByteProcessor;
import ij.process.ColorProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.ShortProcessor;

/**
 * <p>Aims to drastically reduce the complexity of, and time taken to, 
 * fit a sphere to a sphere-like object such as a femoral head.  Requires 
 * one input point inside the sphere.</p>
 * 
 * <p>Potential issues: false positives such as holes, well-defined trabeculae, etc.
 * There is some attempt to account for these, as well as image artefacts.</p>
 * 
 * <p>Can switch resolution to 'unit's (generally mm) with commented out lines.</p>
 * 
 * @author Nick Powell
 *
 */
public class SphereEdgeGuesser implements PlugIn, MouseListener {
	
	private ImageProcessor[] sliceProcessors = null;
	private ImageCanvas canvas;
	
	/** For testing: only measure a single slice of a multi-slice image */
	private boolean singleSlice = false;
	/** Show the final points used to generate the sphere on a new image */
	private boolean showPoints = true;
	/** Show a sample plot */
	private boolean doPlot = false;
	
	/** Linear unit of measure */
	private String units;
	/** Pixel dimensions in 'unit's */
	private double vH, vW, vD;
	/** Unit vector sizes, -1 < value < 1 */
	private double uX, uY, uZ;
	/** Standard deviation values */
	private double sdX, sdY, sdZ;
	
	/** Initial (iX), current (nX) coordinates (in pixels or 'unit's) */
	private int iX, iY, iZ, nX, nY, nZ;
	/** Image dimensions in 'unit's */
	private int d, h, w;
	
	/** Pixel values along a vector. Expands for vector length. */
	private ArrayList<Float> vectorPixelValues;
	/** Array (known length) of ArrayLists (differing lengths) 
	 * containing pixel values along each vector */
	private ArrayList[] pixelValues;
	
	/** Holds (x, y, z) centre and radius, in same units as those given to fitSphere */
	private double[] sphereDim = new double[4];
	/** User's chosen start point */
	private double[] initialPoint = new double[3];
	/** List of median pixel values of each vector */
	private double[] medianValues;
	/** List of number of steps to bone boundary (median crossing). */
	private int[] boundarySteps;
	/** Distance in 'unit's to the bone boundary given by boundarySteps */
	private double[] boundaryDistances;
	/** List of magnitudes of each unit vector, in 'unit's  */
	private double[] magnitudes;
	/** 3D centroid */
	private double[] centroid;
	/** List of distances to coordinates from a centroid */
	private double[] coOrdDistances;
	
	/** Mean distance */
	private double meanDistance, meanCoOrdDistance;
	/** Standard deviation of all the distances */
	private double sdDistance, sdCoOrdDistance;
	/** Adjusted limit (from median) */
	private double adjustedLimit;
	
	/** Number of vectors to create */
	private int numVectors = 1000;
	/** Weight uZ by this factor */
	private double bias_uZ = 0;
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
	private double[][] coOrdArray, coOrdArrayT, coOrdArrayRef, coOrdArrayRef2;
	
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
//		if (!ic.isMultiSlice(imp)) {
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
		Calibration cal = imp.getCalibration();
		
		vW = cal.pixelWidth;
		vH = cal.pixelHeight;
		vD = cal.pixelDepth;
		units = cal.getUnits();
		
		/* Optionally show a dialogue */
//		GenericDialog gd = new GenericDialog("Options");
//		gd.addNumericField("Create", numVectors, 0, 4, "vectors");
//		gd.addNumericField("Use vectors with <=", stepLimit, 0, 2, "* mean steps");
//		gd.addCheckbox("Single slice (don't fit sphere)", this.singleSlice);
//		gd.addCheckbox("Show points", this.showPoints);
//		gd.addCheckbox("Plot a vector's profile", this.doPlot);
//		gd.showDialog();
//		if(gd.wasCanceled()) { return; }
//		this.numVectors = (int) gd.getNextNumber();
//		this.stepLimit = (int) gd.getNextNumber();
//		this.singleSlice = gd.getNextChoice();
//		this.showPoints = gd.getNextChoice();
//		this.doPlot = gd.getNextChoice();
		
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

		
		double biases[] = new double[2];
		biases[0] = 0;			// Run 1: vectors in current slice only
		biases[1] = 1;			// Run 2: 3D vectors
		
		for(int i = 0; i < biases.length; i++) {
			
			bias_uZ = biases[i];
			
			projectVectors(numVectors, bias_uZ);
			
			findBoundaries();
			
			/* Add boundary coordinates to ArrayList (send method a new ArrayList) */
//			newCoOrdinates = new ArrayList<double[]>();
			newCoOrdinates = addNewCoOrdinatesToArrayList(new ArrayList<double[]>());
			
			switch (i) {
			
				case 0 : {
					
					refineCoOrdinates(newCoOrdinates, 0);
					
					/* Calculate mean coordinate distance (weighted by outliers). */
					meanCoOrdDistance = Centroid.getCentroid(this.coOrdDistances);
					sdCoOrdDistance = Math.sqrt(ShaftGuesser.variance(this.coOrdDistances));
					
					refineCoOrdinates(newCoOrdinates, 1);
					refineCoOrdinates(newCoOrdinates, 0);
					
					/* Recalculate mean coordinate distance (less weighted by outliers). */
					meanCoOrdDistance = Centroid.getCentroid(this.coOrdDistances);
					sdCoOrdDistance = Math.sqrt(ShaftGuesser.variance(this.coOrdDistances));
					
					break;
				}
				case 1 : {
					refineCoOrdinates(newCoOrdinates, 0);
					refineCoOrdinates(newCoOrdinates, 2);
					break;
				}
			}
			
			/* Append the refined newCoOrdinates to coOrdinates */
			if(coOrdinates == null) {
				coOrdinates = new ArrayList<double[]>();
			}
			coOrdinates.addAll(newCoOrdinates);
			
			/* Reset initial point based on the above */
//			setInitialPoint((int) centroid[0], (int) centroid[1], (int) centroid[2]);
		}
		
		Iterator<double[]> it = coOrdinates.iterator();
		int i = 0;
		while (it.hasNext()) {
			double[] a = it.next();
			i++;
		}
		
		this.coOrdArrayRef = (double[][]) coOrdinates.toArray(new double[coOrdinates.size()][3]);
		
		if(showPoints) {
			/* Show all detected boundaries */
//			annotateImage(imp, coOrdArray).show();
			/* Show points used by FitSphere */
			annotateImage(imp, coOrdArrayRef).show();
		}
		
		for(int n = 0; n < coOrdArrayRef.length; n++) {
			IJ.log("Ref: " + coOrdArrayRef[n][0] + ", " + coOrdArrayRef[n][1] + ", " + coOrdArrayRef[n][2] + "");
			coOrdArrayRef[n][0] = coOrdArrayRef[n][0] * vW;
			coOrdArrayRef[n][1] = coOrdArrayRef[n][1] * vH;
			coOrdArrayRef[n][2] = coOrdArrayRef[n][2] * vD;
		}
		
		/* Fit sphere to points.
		 * Nb. This feeds fitSphere with pixels, whereas SphereFitter feeds it 'unit's.
		 * The output is thus in pixels. */
		if(!singleSlice) {
			try {
				this.sphereDim = FitSphere.fitSphere(coOrdArrayRef);
			} catch (IllegalArgumentException ia) {
				IJ.showMessage(ia.getMessage());
				return;
			} catch (RuntimeException re) {
				IJ.showMessage("Can't fit sphere to points.\n"
							+ "Rules may be too strict. Try again.");
				return;
			}
		}
		
		/* Show results */
		ResultsTable rt = ResultsTable.getResultsTable();
		rt.incrementCounter();
		rt.addValue("Points used", coOrdinates.size());
		rt.addValue("X centroid (" + units + ")", sphereDim[0]);
		rt.addValue("Y centroid (" + units + ")", sphereDim[1]);
		rt.addValue("Z centroid (approx. slice)", sphereDim[2] / vD);
//		rt.addValue("Z centroid (" + units + ")", sphereDim[2]);
		rt.addValue("Radius (" + units + ")", sphereDim[3]);
		rt.addValue("sdX (pixels)", sdX);
		rt.addValue("sdY (pixels)", sdY);
		rt.addValue("sdZ (pixels)", sdZ);
//		rt.addValue("mX (mm)", meanCoOrds[0] * vW);
//		rt.addValue("mY (mm)", meanCoOrds[1] * vH);
//		rt.addValue("mZ (slice)", meanCoOrds[2]);
		rt.show("Results");
		
		if(doPlot) {
			/* For plotting pixel values */
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
		
//		this.run(arg);

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
//		Calibration cal = imp.getCalibration();
		int x = canvas.offScreenX(e.getX());
		int y = canvas.offScreenY(e.getY());
		int z = imp.getCurrentSlice();
//		final double[] initialPoint = { x * vW, y * vH, z * vD };	// for 'unit' resolution
		setInitialPoint(x, y, z);
	}

	public void mouseReleased(MouseEvent e) { }
	public void mouseExited(MouseEvent e) { }
	public void mouseClicked(MouseEvent e) { }
	public void mouseEntered(MouseEvent e) { }
	public void mouseMoved(MouseEvent e) { }
	
	public void setInitialPoint(int x, int y, int z) {
		this.initialPoint[0] = x;
		this.initialPoint[1] = y;
		this.initialPoint[2] = z;
		
		this.iX = (int) Math.floor(initialPoint[0]);
		this.iY = (int) Math.floor(initialPoint[1]);
		this.iZ = (int) Math.floor(initialPoint[2]);
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
	 * the head. This is taken to be the most significant bone boundary which the 
	 * vector crosses, i.e, the outer surface of the bone.</p>
	 * 
	 * <p>In case narrow surfaces (e.g. scanning artefacts) are found accidentally, 
	 * the code looks 5% further along the vector: if this was a narrow peak, 
	 * it is discounted.</p>
	 * 
	 * <p>Modified ShaftGuesser.shaftLimiter</p>
	 * 
	 * @param pxlValues
	 * @param limit
	 * @param startPoint begin here: usually (pxlValues.length - 1)
	 * @param isMin Specify whether the limit is a minimum (true) or a maximum (false).
	 * @return
	 */
	private int boundaryLimiter(double[] pxlValues, double limit, int startPoint, boolean isMin) {
		
		int boundaryPosition = 0;
		int twoPc = (int) Math.floor(pxlValues.length * 0.02);
		int fivePc = (int) Math.floor(pxlValues.length * 0.05);
		int tenPc = (int) Math.floor(pxlValues.length * 0.1);
		
		/* Start from the end of the vector */
		for(int i = startPoint; i >= 0; i--) {
			if(isMin && pxlValues[i] < limit) {
				boundaryPosition = i;
				/* Look ahead */
				if(i - tenPc > 0 && (pxlValues[i - tenPc] > limit || pxlValues[i - fivePc] > limit) && pxlValues[i - twoPc] > limit) {
					boundaryPosition = boundaryLimiter(pxlValues, limit, (i - twoPc), isMin);
				}
				else { break; }
			}
			else if(!isMin && pxlValues[i] > limit) {
				boundaryPosition = i;
				/* Look ahead */
				if(i - tenPc > 0 && (pxlValues[i - tenPc] < limit || pxlValues[i - fivePc] < limit) &&  pxlValues[i - twoPc] < limit) {
					boundaryPosition = boundaryLimiter(pxlValues, limit, (i - twoPc), isMin);
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
	 * @param coOrdinates
	 * @return ArrayList<double[]> coOrdinates with new coordinates appended.
	 */
	private ArrayList<double[]> addNewCoOrdinatesToArrayList(ArrayList<double[]> newCoOrdinates) {
		
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
	 * <p>Includes math for refinement.</p>
	 * 
	 * @param caseNum Run number, sets which math to use on the data
	 */
	private void refineCoOrdinates(ArrayList<double[]> coOrdinates, int caseNum) {
		
		/* Do different things based on which run this is */
		switch (caseNum) {
		
		/* Coordinate distances */
		case 0 : {
			
			/* Calculate coordinate distances (will be slightly lower than 
			 * boundaryDistances as coOrdinates are Math.floored). */
			this.coOrdDistances = new double[coOrdinates.size()];
			
			Iterator<double[]> itA = coOrdinates.iterator();
			int i = 0;
			while (itA.hasNext()) {
				double[] a = itA.next();
//				coOrdDistances[i] = Trig.distance3D(a, initialPoint) * vW;
				coOrdDistances[i] = Math.sqrt(Math.pow((a[0] - initialPoint[0]) * vW, 2) 
						+ Math.pow((a[1] - initialPoint[1]) * vH, 2) 
						+ Math.pow((a[2] - initialPoint[2]) * vD, 2));
				i++;
			}

			break;
		}
		
		/* 2D */
		case 1 : {
			
			/* Remove large outliers */
			Iterator<double[]> itB = coOrdinates.iterator();
			int j = 0;
			while (itB.hasNext()) {
				double[] b = itB.next();
				if(coOrdDistances[j] > meanCoOrdDistance + (sdCoOrdDistance * 0.5)) {
					itB.remove();
				}
				j++;
			}
			
			break;
		}
		
		/* 3D relies upon 2D */
		case 2 : {
			
			/* Refine based on case 1 */
			Iterator<double[]> itA = coOrdinates.iterator();
			int i = 0;
			while (itA.hasNext()) {
				double[] a = itA.next();
				if(coOrdDistances[i] > meanCoOrdDistance + (sdCoOrdDistance * 1)
						 || coOrdDistances[i] < meanCoOrdDistance + (sdCoOrdDistance * -1)) {
					itA.remove();
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
}
