package org.doube.bonej;

import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.util.ArrayList;
import java.util.Random;

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
import ij.plugin.PlugIn;
import ij.plugin.filter.Profiler;
import ij.plugin.frame.RoiManager;
import ij.process.ImageProcessor;

/**
 * Aims to drastically reduce the complexity of and time take to fit a sphere.
 * Requires one input point inside the sphere.
 * 
 * @author Nick Powell
 *
 */
public class SphereEdgeGuesser implements PlugIn, MouseListener {
	
	private ImageCanvas canvas;
	
	private int currentSlice;
	
	/** Initial (iX), current (nX) and final (fX) positions */
	private int iX, iY, iZ, nX, nY, nZ, fX, fY, fZ;
	/** Unit vector sizes */
	private double uX, uY, uZ;
	
	/** List of distances */
	
	/** Pixel values along a vector. Expands for vector length. */
	private ArrayList<Integer> vectorPixelValues;
	/** Array (known length) of ArrayLists (different lengths) 
	 * containing pixel values along each vector */
	private ArrayList[] pixelValues;
	
	private double[] sphereDim = new double[4];
	private double[] profile;
	private double[] xValues;
	/** List of total lengths */
	private double[] distances;
	private double[][] pixelValues2;
	/** User's chosen start point */
	private double[] initialPoint = new double[3];
	
	/** Random unit vectors (x, y, z) */
	private double[][] unitVectors;
	
	/** Fill with random edge points */
	private double[][] edgePoints;
	/** List of x, y, z coordinates at each point along each line */
	private double[][][] coOrdinates;
	
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
//			IJ.error("Stack required");
//			return;
//		}
		
		/* Don't throw an error if RoiManager isn't open; just open it */
//		RoiManager roiMan = RoiManager.getInstance();
//		if (roiMan == null && imp != null) {
//			IJ.run("ROI Manager...");
//			return;
//		}
		
//		ImageWindow win = imp.getWindow();
//		this.canvas = win.getCanvas();
		
//		Roi roi = imp.getRoi();
		Calibration cal = imp.getCalibration();
		ImageProcessor ip = imp.getProcessor();
		
		/** Pixel dimensions */
		final double vW = cal.pixelWidth;
		final double vH = cal.pixelHeight;
		final double vD = cal.pixelDepth;
		
//		GenericDialog gd = new GenericDialog("Options");
//		gd.showDialog();
//		if (gd.wasCanceled()) {
//			return;
//		}
		
		/* Fill initialPoint with the coordinates of the first point in RoiMan */
//		this.initialPoint = RoiMan.getRoiManPoints(imp, roiMan)[0];
		
		// remove stale MouseListeners
		MouseListener[] l = this.canvas.getMouseListeners();
		for (int n = 0; n < l.length; n++) {
			this.canvas.removeMouseListener(l[n]);
		}
		// add a new MouseListener
		this.canvas.addMouseListener(this);
		new WaitForUserDialog("Click inside the femoral head.\n" 
				+ "Then hit \'OK\'").show();
		this.canvas.removeMouseListener(this);
		
//		final double[] iP = getInitialPoint();
		
		final int iX = (int) Math.floor(initialPoint[0]);
		final int iY = (int) Math.floor(initialPoint[1]);
		final int iZ = (int) Math.floor(initialPoint[2]);
		
		/* Confirm initial point */
		GenericDialog gd = new GenericDialog("Info");
		gd.addMessage("Initial[0] " + iX + " 1 " + iY + " 2 " + iZ + "");
		gd.showDialog();
		
		ImageStack stack = imp.getStack();
		
		/* Create array of random 3D unit vectors */
		this.unitVectors = Vectors.random3D(100);
		
		/* Set up array of ImageProcessors for reference later */
		ImageProcessor[] sliceProcessors = new ImageProcessor[imp.getStackSize() + 1];
		for(int s = 1; s < sliceProcessors.length; s++) {
			sliceProcessors[s] = stack.getProcessor(s);
		}
		
		final int w = (int) Math.floor(imp.getWidth() * vW);
		final int h = (int) Math.floor(imp.getHeight() * vH);
		final int d = (int) Math.floor(imp.getStackSize() * vD);
		
		/* Array of ArrayLists */
		this.pixelValues = new ArrayList[unitVectors.length];
		
		/* Cycle through all vectors */
		for(int i = 0; i < unitVectors.length; i++) {
			
			uX = unitVectors[i][0];
			uY = unitVectors[i][1];
			uZ = unitVectors[i][2];
			
			/* Use ArrayList as we don't know how long each line (j) will be.
			 * Unfortunately, ArrayList can only be filled with Objects. */
			vectorPixelValues = new ArrayList<Integer>();
			
			int j = 0;
			while(j > -1) {
				
				nX = (int) Math.floor(iX + (j * uX));
				nY = (int) Math.floor(iY + (j * uY));
				nZ = (int) Math.floor(iZ + (j * uZ));
				
				/* Break out if outside the image */
				if(nX < 0 || nX > w || nY < 0 || nY > h || nZ < 0 || nZ > d) {
					j = -2;
					break;
				}
				
//				Integer uKK = new Integer(sliceProcessors[nZ].getPixel(nX, nY));
				// may not work
//				pixelValues2[i][j] = aL.add(uKK);
				
				vectorPixelValues.add(sliceProcessors[nZ].getPixel(nX, nY));
				j++;
			}
			
			pixelValues[i] = vectorPixelValues;
		}
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
//		/* List random points along stack edges */
//		this.edgePoints = new double[100][3];
//		Random r = new Random();
//		for(int i = 0; i < edgePoints.length; i++) {
//			// top edge, y = 0
//			if(i < edgePoints.length * 0.25) {
//				edgePoints[i][0] = r.nextInt(w);
//				edgePoints[i][1] = 0;
//				edgePoints[i][2] = r.nextInt(d);
//			}
//			// bottom edge, y = height
//			else if(i < edgePoints.length * 0.5) {
//				edgePoints[i][0] = r.nextInt(w);
//				edgePoints[i][1] = h;
//				edgePoints[i][2] = r.nextInt(d);
//			}
//			// left edge, x = 0
//			else if(i < edgePoints.length * 0.75) {
//				edgePoints[i][0] = 0;
//				edgePoints[i][1] = r.nextInt(h);
//				edgePoints[i][2] = r.nextInt(d);
//			}
//			// right edge, x = width
//			else {
//				edgePoints[i][0] = w;
//				edgePoints[i][1] = r.nextInt(h);
//				edgePoints[i][2] = r.nextInt(d);
//			}
//		}
		
		
		
		/** Total distances */
		double xD, yD, zD;
		
		this.distances = new double[edgePoints.length];
		for(int i = 0; i < edgePoints.length; i++) {
			xD = edgePoints[i][0] - iX;
			yD = edgePoints[i][1] - iY;
			zD = edgePoints[i][2] - iZ;
			
			distances[i] = Math.sqrt((xD * xD) + (yD * yD) + (zD * zD));
		}
		
		/** Interpolatable step distances - can be negative */
		double dX, dY, dZ;
		
		/* Ragged array */
		this.pixelValues = new double[edgePoints.length][];		// Array of rows
		this.coOrdinates = new double[edgePoints.length][][];
		
		for(int i = 0; i < pixelValues.length; i++) {
			
			pixelValues[i] = new double[(int) distances[i]];	// One row, nb. cast to int
			// consider rounding with: int number = Convert.ToInt32(doubleValue);
			coOrdinates[i] = new double[pixelValues[i].length][3];
			
			dX = (edgePoints[i][0] - iX) * 1/pixelValues[i].length;
			dY = (edgePoints[i][1] - iY) * 1/pixelValues[i].length;
			dZ = (edgePoints[i][2] - iZ) * 1/pixelValues[i].length;
			
			for(int j = 0; j < pixelValues[i].length; j++) {
				
				coOrdinates[i][j][0] = iX + (j * dX);
				coOrdinates[i][j][1] = iY + (j * dY);
				coOrdinates[i][j][2] = iZ + (j * dZ);
				
				IJ.log("line length: " + pixelValues[i].length + "; coOrdinates: (" + coOrdinates[i][j][0] + "," + coOrdinates[i][j][1] + ", "
						+ coOrdinates[i][j][2] + ")");
				
				// This may be messy
				if(edgePoints[i][2] - iZ < 0) {
					currentSlice = (int) Math.floor((iZ - (j * dZ)) / vD);
				}
				else {
					currentSlice = (int) Math.floor((iZ + (j * dZ)) / vD);
				}
				
				ImageProcessor sliceIP = imp.getImageStack().getProcessor(currentSlice);
				
				pixelValues[i][j] = sliceIP.getInterpolatedPixel(coOrdinates[i][j][0], coOrdinates[i][j][1]);
			}
		}
		
		/* For plotting */
		this.xValues = new double[pixelValues[1].length];
		for (int j = 0; j < xValues.length; j++) {
			xValues[j] = (double) j;
		}
		Plot omgPlot = new Plot("OMG startpoint" + iX + " " + iY + " " + iZ + " endpoint" + edgePoints[1][0] + " " + edgePoints[1][1] + " " + edgePoints[1][2] + " ", "distance", "value", xValues, pixelValues[1]);
		omgPlot.show();
		
		this.profile = new ProfilePlot(imp).getProfile();
		
//		/* For plotting */
//		this.xValues = new double[profile.length];
//		for (int j = 0; j < xValues.length; j++) {
//			xValues[j] = (double) j;
//		}
		
//		ProfilePlot pp = new ProfilePlot();
//		Profiler p = new Profiler();
//		this.profile = pp.getProfile();
//		Plot profilePlot = new Plot("Profile plot", "distance", "value", xValues, profile);
//		profilePlot.show();
		
//		pp.getStraightLineProfile(roi, cal, ip);
		
		return;
	}
	
	/**
	 * 
	 * @return
	 */
	public double[] getSphereDim() {
		return this.sphereDim;
	}
	
//	public double[] getInitialPoint() {
//		return this.initialPoint;
//	}
	public void setInitialPoint(double[] a) {
		this.initialPoint = a;
	}
	
	public void mousePressed(MouseEvent e) {
		ImagePlus imp = IJ.getImage();
		Calibration cal = imp.getCalibration();
		int x = canvas.offScreenX(e.getX());
		int y = canvas.offScreenY(e.getY());
		int z = imp.getCurrentSlice();
		final double[] initialPoint = { x * cal.pixelWidth, y * cal.pixelHeight,
				z * cal.pixelDepth };
//		IJ.log("neckPoint: (" + initialPoint[0] + "," + initialPoint[1] + ", "
//				+ initialPoint[2] + ")");
		setInitialPoint(initialPoint);
//		calculateAngles(imp, neckPoint);
//		canvas.removeMouseListener(this);
	}

	public void mouseReleased(MouseEvent e) { }
	public void mouseExited(MouseEvent e) { }
	public void mouseClicked(MouseEvent e) { }
	public void mouseEntered(MouseEvent e) { }
	public void mouseMoved(MouseEvent e) { }

}
