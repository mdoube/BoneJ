package org.doube.bonej;

import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Random;

import org.doube.geometry.Centroid;
import org.doube.geometry.FitSphere;
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
	
	/** Linear unit of measure */
	private String units;
	/** Pixel dimensions in 'unit's */
	private double vH, vW, vD;
	/** Unit vector sizes, -1 < value < 1 */
	private double uX, uY, uZ;
	/** Initial (iX), current (nX) and final (fX) coordinates (in pixels or 'unit's) */
	private int iX, iY, iZ, nX, nY, nZ, fX, fY, fZ;
	/** Image dimensions in 'unit's */
	private int d, h, w;
	/** Slice number of current z coordinate, nZ */
	private int currentSlice;
	
	/** Pixel values along a vector. Expands for vector length. */
	private ArrayList<Float> vectorPixelValues;
	/** Array (known length) of ArrayLists (differing lengths) 
	 * containing pixel values along each vector */
	private ArrayList[] pixelValues;
	
	/** Holds (x, y, z) centre and radius, in same units as those given to fitSphere */
	private double[] sphereDim = new double[4];
	
	/** List of distances along each vector from initial point, before median is crossed */
	private double[] distances;
	/** User's chosen start point */
	private double[] initialPoint = new double[3];
	/** List of median pixel values of each vector */
	private double[] medianValues;
	/** List of number of steps to bone boundary (median crossing). */
	private int[] boundarySteps;
	
	/** Mean distance */
	private double meanDistance;
	/** Mean steps */
	private double meanSteps;
	/** Adjusted limit (from median) */
	private double adjustedLimit;
	
	/** Number of vectors to create */
	private int numVectors = 1000;
	/** Only use vectors which take <= this multiple of the mean number of steps 
	 * taken by all vectors to reach the bone surface. */
	private double stepLimit = 0.3;
	/** Random unit vectors (x, y, z) */
	private double[][] unitVectors;
	/** Ragged array of pixel values for each vector */
	private double[][] pxVals;
	/** Ragged array of x-axis values for plotting pixel value profiles along vectors */
	private double[][] xValues;
	
	/** Coordinates x, y, z  of the boundary  */
	private double[] coOrds;
	/** ArrayList (of length limited by meanDistance) of boundary coordinates x, y, z */
	private ArrayList<double[]> coOrdinates;
	/** Array of x, y, z coordinates */
	private double[][] coOrdArray;
	
	public void run(String arg) {
		
		/**
		 * Copy image
		 * Binarise???
		 * 
		 */
		
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
//		gd.showDialog();
//		if(gd.wasCanceled()) {
//			return;
//		}
//		this.numVectors = (int) gd.getNextNumber();
//		this.stepLimit = (int) gd.getNextNumber();
		
		// remove stale MouseListeners
		MouseListener[] l = this.canvas.getMouseListeners();
		for (int n = 0; n < l.length; n++) {
			this.canvas.removeMouseListener(l[n]);
		}
		// add a new MouseListener
		this.canvas.addMouseListener(this);
		new WaitForUserDialog("Click inside the femoral head,\n" 
									+ "then hit \'OK\'").show();
		this.canvas.removeMouseListener(this);
		
		iX = (int) Math.floor(initialPoint[0]);
		iY = (int) Math.floor(initialPoint[1]);
		iZ = (int) Math.floor(initialPoint[2]);
		
		/* Confirm initial point */
//		GenericDialog gd1 = new GenericDialog("Info");
//		gd1.addMessage("Initial point (x,y,z): (" + iX + ", " + iY + ", " + iZ + ")");
//		gd1.addMessage("Pixel sizes (vW,vH,vD): (" + vW + ", " + vH + ", " + vD + ")");
//		gd1.addMessage("imp.getCurrentSlice: " + imp.getCurrentSlice() + "");
//		gd1.showDialog();
		
		ImageStack stack = imp.getStack();
		
		/* Set up array of ImageProcessors for reference later */
		this.sliceProcessors = new ImageProcessor[imp.getStackSize() + 1];
		for(int s = 1; s < sliceProcessors.length; s++) {
			sliceProcessors[s] = stack.getProcessor(s);
		}
		
		w = imp.getWidth();
		h = imp.getHeight();
		d = imp.getStackSize();
//		w = (int) Math.floor(imp.getWidth() * vW);			// for 'unit' resolution
//		h = (int) Math.floor(imp.getHeight() * vH);
//		d = (int) Math.floor(imp.getStackSize() * vD);
		
//		GenericDialog gd2 = new GenericDialog("Info");
//		gd2.addMessage("w,h,d: " + w + "," + h + "," + d + "");
//		gd2.showDialog();
		
		/* Create array of random 3D unit vectors */
		this.unitVectors = Vectors.random3D(numVectors);
		
		/* Array of ArrayLists */
		this.pixelValues = new ArrayList[unitVectors.length];
		
		/* Cycle through all vectors */
		for(int i = 0; i < unitVectors.length; i++) {
			
			uX = unitVectors[i][0];
			uY = unitVectors[i][1];
			uZ = unitVectors[i][2];
//			uZ = 0;						// for testing with single slice
			
			/* Use ArrayList as we don't know how long each line (j) will be.
			 * Unfortunately, ArrayList can only be filled with Objects. */
			vectorPixelValues = new ArrayList<Float>();
			
			int j = 0;
			while(j > -1) {
				
				nX = (int) Math.floor(iX + (j * uX));
				nY = (int) Math.floor(iY + (j * uY));
				nZ = (int) Math.round(iZ + (j * uZ));
				
				/* Break out if outside the image */		// also || (nZ / vD) < 1 if using 'unit's resolution
				if(nX < 0 || nX > w || nY < 0 || nY > h || nZ < 1 || nZ > d) {
					j = -2;
					break;
				}
				
				int nowX = nX;
				int nowY = nY;
				this.currentSlice = nZ;
//				int nowX = (int) Math.floor(nX / vW);		// for 'unit' resolution
//				int nowY = (int) Math.floor(nY / vH);
//				this.currentSlice = (int) Math.floor(nZ / vD);
				
//				IJ.log("(nX,nY,nZ): (" + nX + "," + nY + "," + nZ + "), slice: " + currentSlice + "; nowX, nowY: " + nowX + "," + nowY + "");
				
				Float nowValue = new Float(sliceProcessors[currentSlice].getPixelValue(nowX, nowY));
				
//				IJ.log("(nX,nY,nZ): (" + nX + "," + nY + "," + nZ + "); slice: " + currentSlice + "; int: " + nowValue + "");
				
				vectorPixelValues.add(nowValue);
				j++;
			}

			pixelValues[i] = vectorPixelValues;
		}

		/* For plotting */
		this.xValues = new double[pixelValues.length][];
		for(int i = 0; i < xValues.length; i++) {
			xValues[i] = new double[pixelValues[i].size()];
			for (int j = 0; j < pixelValues[i].size(); j++) {
				xValues[i][j] = j;
			}
		}
		
		/* Some analysis possible on each vector */
		this.medianValues = new double[pixelValues.length];
		this.boundarySteps = new int[pixelValues.length];
		this.distances = new double[pixelValues.length];
		
		/* Ragged array of (double) pixel values for every vector */
		this.pxVals = new double[pixelValues.length][];
		
		for(int i = 0; i < pxVals.length; i++) {
			pxVals[i] = new double[pixelValues[i].size()];
			for(int j = 0; j < pixelValues[i].size(); j++) {
				pxVals[i][j] = (double) (Float) pixelValues[i].toArray()[j];
			}
			medianValues[i] = ShaftGuesser.median(pxVals[i]);
			
			/* Adjust our guessing limit based on the typical plots seen */
			if(medianValues[i] < 0) {
				adjustedLimit = medianValues[i] * 0.5;
			}
			else {
				adjustedLimit = medianValues[i] * 1.25;
			}
			boundarySteps[i] = boundaryLimiter(pxVals[i], adjustedLimit, (pxVals[i].length - 1), false);
			distances[i] = boundarySteps[i] * Math.sqrt((unitVectors[i][0] * unitVectors[i][0]) + (unitVectors[i][1] * unitVectors[i][1]) + (unitVectors[i][2] * unitVectors[i][2]));
		}
		
		Plot aPlot = new Plot("Pixel Values along vector " + 1 + "", "Distance (" + units + ")", "Pixel value", xValues[0], pxVals[0]);
		aPlot.show();
		
		/* */
		coOrdinates = new ArrayList<double[]>();
		this.meanDistance = Centroid.getCentroid(distances);
		this.meanSteps = Centroid.getCentroid(boundarySteps);
		
		IJ.log("meanSteps: " + meanSteps + "");
		
		/* Refinement */
		
		/* Get ArrayList of boundary coordinates, 
		 * limited by meanSteps to boundary edge. */ 
		for(int i = 0; i < boundarySteps.length; i++) {
//			if(true) {
			if(boundarySteps[i] <= (meanSteps * stepLimit)
					&& boundarySteps[i] >= (meanSteps * stepLimit / 2)) {
				
				this.coOrds = new double[3];
				
				fX = (int) Math.floor(iX + (boundarySteps[i] * unitVectors[i][0]));
				fY = (int) Math.floor(iY + (boundarySteps[i] * unitVectors[i][1]));
				fZ = (int) Math.floor(iZ + (boundarySteps[i] * unitVectors[i][2]));
//				fZ = (int) Math.floor(iZ + (boundarySteps[i] * 0));		// for testing with single slice
				
				coOrds[0] = fX;
				coOrds[1] = fY;
				coOrds[2] = fZ;
//				coOrds[0] = (int) Math.floor(fX / vW);	// for 'unit' resolution
//				coOrds[1] = (int) Math.floor(fY / vH);
//				coOrds[2] = (int) Math.floor(fZ / vD);
				coOrdinates.add(coOrds);
				
//				IJ.log("(fX,fY,fZ): (" + fX + ", " + fY + ", " + fZ + ")");
//				IJ.log("coOrdinates.size(): " + coOrdinates.size() + "; values: " + coOrdinates.get(i)[0] + "");
				IJ.log("boundarySteps: " + boundarySteps[i] + "; coOrds: " + fX + ", " + fY + ", "  + fZ + "; mm X, Y: " + fX * vW + ", " + fY * vH + ", " + fZ + ", ");
			}
		}
		
		/* Transfer to array (and transpose) for math (need rows) */
		this.coOrdArray = transposeArray((double[][]) coOrdinates.toArray(new double[coOrdinates.size()][3]));
		
		/* More refinement: standard deviations */
		double sdX = Math.sqrt(ShaftGuesser.variance(coOrdArray[0]));
		double sdY = Math.sqrt(ShaftGuesser.variance(coOrdArray[1]));
		double sdZ = Math.sqrt(ShaftGuesser.variance(coOrdArray[2]));
		
		double[] meanCoOrds = new double[3];
		meanCoOrds[0] = Centroid.getCentroid(coOrdArray[0]);
		meanCoOrds[1] = Centroid.getCentroid(coOrdArray[1]);
		meanCoOrds[2] = Centroid.getCentroid(coOrdArray[2]);
		
		/* Iterate through coOrdinates ArrayList and refine */
//		Iterator<double[]> it = coOrdinates.iterator();
//		while (it.hasNext()) {
//			double[] a = it.next();
//			if(Math.abs(a[0] - meanCoOrds[0]) > sdX * 0.5
//					|| Math.abs(a[1] - meanCoOrds[1]) > sdY * 0.5
//					|| Math.abs(a[2] - meanCoOrds[2]) > sdZ * 0.5) {
//				it.remove();
//			}
//		}
		
		
		
		
		
		
		/* Fit sphere to points.
		 * Nb. This feeds fitSphere with pixels, whereas SphereFitter feeds it 'unit's.
		 * The output is thus in pixels. */
		this.sphereDim = FitSphere.fitSphere((double[][]) coOrdinates.toArray(new double[coOrdinates.size()][3]));
		
//		GenericDialog gdSphereDims = new GenericDialog("Sphere dimensions");
//		gdSphereDims.addMessage("Centre (x, y, z): " + sphereDim[0] + ", " + sphereDim[1] + ", " + sphereDim[2] + "; radius: " + sphereDim[3] + "");
//		gdSphereDims.showDialog();
		
		/* Show results */
		ResultsTable rt = ResultsTable.getResultsTable();
		rt.incrementCounter();
		rt.addValue("coOrdinates.size", coOrdinates.size());
		rt.addValue("X centroid (" + units + ")", sphereDim[0] * vW);
		rt.addValue("Y centroid (" + units + ")", sphereDim[1] * vH);
		rt.addValue("Z centroid (" + units + ")", sphereDim[2] * vD);
		rt.addValue("Z centroid (approx. slice)", sphereDim[2]);
		rt.addValue("Radius (pixels)", sphereDim[3]);
		rt.addValue("sdX (pixels)", sdX);
		rt.addValue("sdY (pixels)", sdY);
		rt.addValue("sdZ (pixels)", sdZ);
		rt.addValue("mX (mm)", meanCoOrds[0] * vW);
		rt.addValue("mY (mm)", meanCoOrds[1] * vH);
		rt.addValue("mZ (mm)", meanCoOrds[2] * vD);
		rt.show("Results");

		return;
	}
	
	/**
	 * 
	 * @return
	 */
	public double[] getSphereDim() {
		return this.sphereDim;
	}
	public void setInitialPoint(double[] a) {
		this.initialPoint = a;
	}
	
	public void mousePressed(MouseEvent e) {
		ImagePlus imp = IJ.getImage();
		Calibration cal = imp.getCalibration();
		int x = canvas.offScreenX(e.getX());
		int y = canvas.offScreenY(e.getY());
		int z = imp.getCurrentSlice();
		final double[] initialPoint = { x, y, z };
//		final double[] initialPoint = { x * vW, y * vH, z * vD };	// for 'unit' resolution
		setInitialPoint(initialPoint);
	}

	public void mouseReleased(MouseEvent e) { }
	public void mouseExited(MouseEvent e) { }
	public void mouseClicked(MouseEvent e) { }
	public void mouseEntered(MouseEvent e) { }
	public void mouseMoved(MouseEvent e) { }

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
	public int boundaryLimiter(double[] pxlValues, double limit, int startPoint, boolean isMin) {
		
		int boundaryPosition = 0;
		int fivePc = (int) Math.floor(pxlValues.length * 0.05);
		int tenPc = (int) Math.floor(pxlValues.length * 0.1);
		
		/* Start from the end of the vector */
		for(int i = startPoint; i >= 0; i--) {
			if(isMin && pxlValues[i] < limit) {
				boundaryPosition = i;
				/* Look 5% ahead */
//				if(i - fivePc > 0 && pxlValues[i - fivePc] > limit) {
//					boundaryPosition = boundaryLimiter(pxlValues, limit, (i - fivePc), isMin);
//				}
				if(i - tenPc > 0 && pxlValues[i - tenPc] > limit&& pxlValues[i - fivePc] > limit) {
					boundaryPosition = boundaryLimiter(pxlValues, limit, (i - tenPc), isMin);
				}
				else { break; }
			}
			else if(!isMin && pxlValues[i] > limit) {
				boundaryPosition = i;
				/* Look 5% ahead */
//				if(i - fivePc > 0 && pxlValues[i - fivePc] < limit) {
//					boundaryPosition = boundaryLimiter(pxlValues, limit, (i - fivePc), isMin);
//				}
				if(i - tenPc > 0 && pxlValues[i - tenPc] < limit && pxlValues[i - fivePc] < limit) {
					boundaryPosition = boundaryLimiter(pxlValues, limit, (i - tenPc), isMin);
				}
				else { break; }
			}
		}
		
		return boundaryPosition;
	}
	
	/**
	 * Provided a 2D inputArray[][] with constant number of columns per row,
	 * i.e. not a ragged array, transposes this array so that row and column 
	 * values are switched.
	 * 
	 * @param inputArray
	 * @return transposedArray[][]
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
}
