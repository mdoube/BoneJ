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
 * @author Nick Powell
 *
 */
public class SphereEdgeGuesser implements PlugIn, MouseListener {
	
	private double fill_value;
	private double get;
	private double getpx;
	private double getpxval;
	private double btdpth;
	
	private ImageProcessor[] sliceProcessors = null;
	private ImageProcessor ip = null;
	private ImageCanvas canvas;
	
	/** Pixel dimensions in 'unit's */
	private double vH, vW, vD;
	/** Initial (iX), current (nX) and final (fX) coordinates in 'unit's */
	private int iX, iY, iZ, nX, nY, nZ, fX, fY, fZ;
	/** Unit vector sizes, -1 < value < 1 */
	private double uX, uY, uZ;
	/** Slice number of current z coordinate, nZ */
	private int currentSlice;
	
	/** List of distances */
	
	/** Pixel values along a vector. Expands for vector length. */
	private ArrayList<Float> vectorPixelValues;
//	private ArrayList<Float> vectorPixelValues;
	/** Array (known length) of ArrayLists (differing lengths) 
	 * containing pixel values along each vector */
	private ArrayList[] pixelValues;
	
	private double[] sphereDim = new double[4];
	private double[] profile;
	
	/** Ragged array of pixel values for each vector */
	private double[][] pxVals;
	/** Ragged array of x-axis values for plotting pixel value profiles along vectors */
	private double[][] xValues;
	
	
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
		
		ImageWindow win = imp.getWindow();
		this.canvas = win.getCanvas();
		Calibration cal = imp.getCalibration();
		
		vW = cal.pixelWidth;
		vH = cal.pixelHeight;
		vD = cal.pixelDepth;
		
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
		GenericDialog gd = new GenericDialog("Info");
		gd.addMessage("Initial point (x,y,z): (" + iX + "," + iY + "," + iZ + ")");
		gd.addMessage("Pixel sizes (vW,vH,vD): (" + vW + "," + vH + "," + vD + ")");
		gd.addMessage("imp.getCurrentSlice" + imp.getCurrentSlice() + "");
		gd.showDialog();
		
		ImageStack stack = imp.getStack();
		
		/* Create array of random 3D unit vectors */
		this.unitVectors = Vectors.random3D(1);
		
		/* Set up array of ImageProcessors for reference later */
		this.sliceProcessors = new ImageProcessor[imp.getStackSize() + 1];
		for(int s = 1; s < sliceProcessors.length; s++) {
			sliceProcessors[s] = stack.getProcessor(s);
		}
		
//		int bitDepth = imp.getBitDepth();
//		switch (bitDepth) {
//			case 8: sliceProcessors = new ByteProcessor[imp.getStackSize() + 1]; 
//				for(int s = 1; s < sliceProcessors.length; s++) {
//					sliceProcessors[s] = (ByteProcessor) stack.getProcessor(s);
//				}
//				break;  
//			case 16: sliceProcessors = new ShortProcessor[imp.getStackSize() + 1]; 
//				for(int s = 1; s < sliceProcessors.length; s++) {
//					sliceProcessors[s] = (ShortProcessor) stack.getProcessor(s);
//				}
//				break;  
//			case 32: sliceProcessors = new FloatProcessor[imp.getStackSize() + 1]; 
//				for(int s = 1; s < sliceProcessors.length; s++) {
//					sliceProcessors[s] = (FloatProcessor) stack.getProcessor(s);
//				}
//				break;
//			case 24: sliceProcessors = new ColorProcessor[imp.getStackSize() + 1]; 
//				for(int s = 1; s < sliceProcessors.length; s++) {
//					sliceProcessors[s] = (ColorProcessor) stack.getProcessor(s);
//				}
//				break;
//		}
		
		/** Image dimensions in 'unit's */
		final int w = (int) Math.floor(imp.getWidth() * vW);
		final int h = (int) Math.floor(imp.getHeight() * vH);
		final int d = (int) Math.floor(imp.getStackSize() * vD);
		
		GenericDialog gd2 = new GenericDialog("Info");
		gd2.addMessage("w,h,d: " + w + "," + h + "," + d + "");
		gd2.showDialog();
		
		/* Array of ArrayLists */
		this.pixelValues = new ArrayList[unitVectors.length];
		
		/* Cycle through all vectors */
		for(int i = 0; i < unitVectors.length; i++) {
			
			uX = unitVectors[i][0];
			uY = unitVectors[i][1];
//			uZ = unitVectors[i][2];
			uZ = 0;
			
			/* Use ArrayList as we don't know how long each line (j) will be.
			 * Unfortunately, ArrayList can only be filled with Objects. */
			vectorPixelValues = new ArrayList<Float>();
			
			int j = 0;
			while(j > -1) {
				
				nX = (int) Math.floor(iX + (j * uX));
				nY = (int) Math.floor(iY + (j * uY));
				nZ = (int) Math.floor(iZ + (j * uZ));
				
				/* Break out if outside the image */
				if(nX < 0 || nX > w || nY < 0 || nY > h || nZ < 1 || nZ > d) {
					j = -2;
					break;
				}
				
				int nowX = (int) Math.floor(nX / vW);
				int nowY = (int) Math.floor(nY / vH);
				
				this.currentSlice = (int) Math.floor(nZ / vD);
				IJ.log("(nX,nY,nZ): (" + nX + "," + nY + "," + nZ + ", slice: " + currentSlice + "; nowX, nowY: " + nowX + "," + nowY + "");
				Float floating = new Float(sliceProcessors[currentSlice].getPixelValue(nowX, nowY));
				
				IJ.log("(nX,nY,nZ): (" + nX + "," + nY + "," + nZ + "); slice: " + currentSlice + "; int: " + floating + "");
				
				vectorPixelValues.add(floating);
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
		/* Ragged array of (double) pixel values for each vector */
		this.pxVals = new double[pixelValues.length][];
		for(int i = 0; i < pxVals.length; i++) {
			pxVals[i] = new double[pixelValues[i].size()];
			for(int j = 0; j < pixelValues[i].size(); j++) {
				pxVals[i][j] = (double) (Float) pixelValues[i].toArray()[j];
			}
		}
		
		Plot aPlot = new Plot("Pixel Values along 1st vector", "Distance (mm)", "Pixel value", xValues[0], pxVals[0]);
		aPlot.show();

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

}
