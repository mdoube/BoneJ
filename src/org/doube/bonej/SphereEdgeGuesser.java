package org.doube.bonej;

import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.util.Random;

import org.doube.util.ImageCheck;
import org.doube.util.RoiMan;

import ij.IJ;
import ij.ImagePlus;
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
	
	private double[] sphereDim = new double[4];
	private double[] profile;
	private double[] xValues;
	/** List of total lengths */
	private double[] distances;
	private double[][] pixelValues;
	/** User's chosen start point */
	private double[] initialPoint = new double[3];
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
		
		ImageWindow win = imp.getWindow();
		this.canvas = win.getCanvas();
		
		Roi roi = imp.getRoi();
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
		
		GenericDialog gd = new GenericDialog("Info");
		gd.addMessage("Hello " + initialPoint[2] + " bye!");
		gd.showDialog();
		
		final int w = (int) Math.floor(imp.getWidth() * vW);
		final int h = (int) Math.floor(imp.getHeight() * vH);
		final int d = (int) Math.floor(imp.getStackSize() * vD);
		
		/* List random points along stack edges */
		this.edgePoints = new double[100][3];
		Random r = new Random();
		for(int i = 0; i < edgePoints.length; i++) {
			// top edge, y = 0
			if(i < edgePoints.length * 0.25) {
				edgePoints[i][0] = r.nextInt(w);
				edgePoints[i][1] = 0;
				edgePoints[i][2] = r.nextInt(d);
			}
			// bottom edge, y = height
			else if(i < edgePoints.length * 0.5) {
				edgePoints[i][0] = r.nextInt(w);
				edgePoints[i][1] = h;
				edgePoints[i][2] = r.nextInt(d);
			}
			// left edge, x = 0
			else if(i < edgePoints.length * 0.75) {
				edgePoints[i][0] = 0;
				edgePoints[i][1] = r.nextInt(h);
				edgePoints[i][2] = r.nextInt(d);
			}
			// right edge, x = width
			else {
				edgePoints[i][0] = w;
				edgePoints[i][1] = r.nextInt(h);
				edgePoints[i][2] = r.nextInt(d);
			}
		}
		
		final double iX = initialPoint[0];
		final double iY = initialPoint[1];
		final double iZ = initialPoint[2];
		
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
		
		this.xValues = new double[pixelValues[2].length];
		for (int j = 0; j < xValues.length; j++) {
			xValues[j] = (double) j;
		}
		Plot omgPlot2 = new Plot("OMG", "distance", "value", xValues, pixelValues[2]);
		omgPlot2.show();
		
		this.xValues = new double[pixelValues[3].length];
		for (int j = 0; j < xValues.length; j++) {
			xValues[j] = (double) j;
		}
		Plot omgPlot3 = new Plot("OMG", "distance", "value", xValues, pixelValues[3]);
		omgPlot3.show();
		
		this.xValues = new double[pixelValues[4].length];
		for (int j = 0; j < xValues.length; j++) {
			xValues[j] = (double) j;
		}
		Plot omgPlot4 = new Plot("OMG", "distance", "value", xValues, pixelValues[4]);
		omgPlot4.show();
		
		this.xValues = new double[pixelValues[5].length];
		for (int j = 0; j < xValues.length; j++) {
			xValues[j] = (double) j;
		}
		Plot omgPlot5 = new Plot("OMG", "distance", "value", xValues, pixelValues[5]);
		omgPlot5.show();
		
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
