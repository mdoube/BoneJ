package org.doube.bonej;

import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.util.ArrayList;
import java.util.Random;

import org.doube.geometry.Centroid;
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
 * Aims to drastically reduce the complexity of and time take to fit a sphere.
 * Requires one input point inside the sphere.
 * 
 * Potential issues: false positives such as holes, well-defined trabeculae, etc.
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
	/** Initial (iX), current (nX) and final (fX) coordinates in 'unit's */
	private int iX, iY, iZ, nX, nY, nZ, fX, fY, fZ;
	/** Image dimensions in 'unit's */
	private int d, h, w;
	/** Unit vector sizes, -1 < value < 1 */
	private double uX, uY, uZ;
	/** Slice number of current z coordinate, nZ */
	private int currentSlice;
	
	/** Pixel values along a vector. Expands for vector length. */
	private ArrayList<Float> vectorPixelValues;
	/** Array (known length) of ArrayLists (differing lengths) 
	 * containing pixel values along each vector */
	private ArrayList[] pixelValues;
	
	private double[] sphereDim = new double[4];
	
	/** List of distances along each vector from initial point, before median is crossed */
	private double[] distances;
	/** User's chosen start point */
	private double[] initialPoint = new double[3];
	/** List of median pixel values of each vector */
	private double[] medianValues;
	/** List of steps to bone boundary (median crossing) */
	private int[] boundarySteps;
	
	/** Mean distance */
	private double meanDistance;
	/** Adjusted limit (from median) */
	private double adjustedLimit;
	
	/** Number of vectors to create */
	private int numVectors;
	/** Random unit vectors (x, y, z) */
	private double[][] unitVectors;
	/** Ragged array of pixel values for each vector */
	private double[][] pxVals;
	/** Ragged array of x-axis values for plotting pixel value profiles along vectors */
	private double[][] xValues;
	
	/** Coordinates x, y, z  of the boundary  */
	private int[] coOrds;
	/** ArrayList (of length limited by meanDistance) of boundary coordinates */
	private ArrayList<int[]> coOrdinates;
	
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
		
		GenericDialog gd = new GenericDialog("Options");
		gd.addNumericField("Create", 1, 0, 4, "vectors");
		gd.showDialog();
		if(gd.wasCanceled()) {
			return;
		}
		this.numVectors = (int) gd.getNextNumber();
		
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
		GenericDialog gd1 = new GenericDialog("Info");
		gd1.addMessage("Initial point (x,y,z): (" + iX + "," + iY + "," + iZ + ")");
		gd1.addMessage("Pixel sizes (vW,vH,vD): (" + vW + "," + vH + "," + vD + ")");
		gd1.addMessage("imp.getCurrentSlice: " + imp.getCurrentSlice() + "");
		gd1.showDialog();
		
		ImageStack stack = imp.getStack();
		
		/* Set up array of ImageProcessors for reference later */
		this.sliceProcessors = new ImageProcessor[imp.getStackSize() + 1];
		for(int s = 1; s < sliceProcessors.length; s++) {
			sliceProcessors[s] = stack.getProcessor(s);
		}
		
		w = (int) Math.floor(imp.getWidth() * vW);
		h = (int) Math.floor(imp.getHeight() * vH);
		d = (int) Math.floor(imp.getStackSize() * vD);
		
		GenericDialog gd2 = new GenericDialog("Info");
		gd2.addMessage("w,h,d: " + w + "," + h + "," + d + "");
		gd2.showDialog();
		
		/* Create array of random 3D unit vectors */
		this.unitVectors = Vectors.random3D(numVectors);
		
		/* Array of ArrayLists */
		this.pixelValues = new ArrayList[unitVectors.length];
		
		/* Cycle through all vectors */
		for(int i = 0; i < unitVectors.length; i++) {
			
			uX = unitVectors[i][0];
			uY = unitVectors[i][1];
			uZ = unitVectors[i][2];
			
			/* Use ArrayList as we don't know how long each line (j) will be.
			 * Unfortunately, ArrayList can only be filled with Objects. */
			vectorPixelValues = new ArrayList<Float>();
			
			int j = 0;
			while(j > -1) {
				
				nX = (int) Math.floor(iX + (j * uX));
				nY = (int) Math.floor(iY + (j * uY));
				nZ = (int) Math.round(iZ + (j * uZ));
				
				/* Break out if outside the image */
				if(nX < 0 || nX > w || nY < 0 || nY > h || (nZ / vD) < 1 || nZ > d) {
					j = -2;
					break;
				}
				
				int nowX = (int) Math.floor(nX / vW);
				int nowY = (int) Math.floor(nY / vH);
				this.currentSlice = (int) Math.floor(nZ / vD);
				
				IJ.log("(nX,nY,nZ): (" + nX + "," + nY + "," + nZ + "), slice: " + currentSlice + "; nowX, nowY: " + nowX + "," + nowY + "");
				
				Float nowValue = new Float(sliceProcessors[currentSlice].getPixelValue(nowX, nowY));
				
				IJ.log("(nX,nY,nZ): (" + nX + "," + nY + "," + nZ + "); slice: " + currentSlice + "; int: " + nowValue + "");
				
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
		
		this.medianValues = new double[pixelValues.length];
		this.boundarySteps = new int[pixelValues.length];
		this.distances = new double[pixelValues.length];
		
		/* Ragged array of (double) pixel values for each vector */
		this.pxVals = new double[pixelValues.length][];
		for(int i = 0; i < pxVals.length; i++) {
			pxVals[i] = new double[pixelValues[i].size()];
			for(int j = 0; j < pixelValues[i].size(); j++) {
				pxVals[i][j] = (double) (Float) pixelValues[i].toArray()[j];
			}
			medianValues[i] = ShaftGuesser.median(pxVals[i]);
			if(medianValues[i] < 0) {
				adjustedLimit = medianValues[i] * 0.5;
			}
			else {
				adjustedLimit = medianValues[i] * 1.25;
			}
			boundarySteps[i] = boundaryLimiter(pxVals[i], adjustedLimit, (pxVals[i].length - 1), false);
			distances[i] = boundarySteps[i] * Math.sqrt((unitVectors[i][0] * unitVectors[i][0]) + (unitVectors[i][1] * unitVectors[i][1]) + (unitVectors[i][2] * unitVectors[i][2]));
		}
		
		Plot aPlot = new Plot("Pixel Values along vector " + 1 + "", "Distance (mm)", "Pixel value", xValues[0], pxVals[0]);
		aPlot.show();
		
		coOrdinates = new ArrayList<int[]>();
		this.coOrds = new int[3];
		this.meanDistance = Centroid.getCentroid(distances);
		
		/* Get list of boundary coordinates (limited by meanDistance to them:
		 * here using those distances <= the mean). */
		for(int i = 0; i < distances.length; i++) {
			if(distances[i] <= meanDistance) {
				
				fX = (int) Math.floor(iX + (i * unitVectors[i][0]));
				fY = (int) Math.floor(iY + (i * unitVectors[i][1]));
				fZ = (int) Math.floor((iZ + (i * unitVectors[i][2])) / vD);
				
				coOrds[0] = (int) Math.floor(fX / vW);
				coOrds[1] = (int) Math.floor(fY / vH);
				coOrds[2] = fZ;
				coOrdinates.add(coOrds);
			}
		}
		
		// Currently, fZ is below by 1.
		// Currently, major issue is mm jumping instead of pixel jumping.
		
		GenericDialog gd4 = new GenericDialog("Median at");
		gd4.addMessage("median: " + medianValues[0] + "");
		gd4.addMessage("distance: " + distances[0] + "");
		gd4.addMessage("fZ: " + fZ + "");
		gd4.showDialog();

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
		final double[] initialPoint = { x * vW, y * vH, z * vD };
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
		
		/* Start from the end of the vector */
		for(int i = startPoint; i >= 0; i--) {
			if(isMin && pxlValues[i] < limit) {
				boundaryPosition = i;
				/* Look 5% ahead */
				if(i - fivePc >= 0) {
					if(pxlValues[i - fivePc] > limit) {
						boundaryPosition = boundaryLimiter(pxlValues, limit, (i - fivePc), isMin);
					}
					else { break; }
				}
				else { break; }
			}
			else if(!isMin && pxlValues[i] > limit) {
				boundaryPosition = i;
				/* Look 5% ahead */
				if(i - fivePc >= 0) {
					if(pxlValues[i - fivePc] < limit) {
						boundaryPosition = boundaryLimiter(pxlValues, limit, (i - fivePc), isMin);
					}
					else { break; }
				}
				else { break; }
			}
		}
		
		return boundaryPosition;
	}
}
