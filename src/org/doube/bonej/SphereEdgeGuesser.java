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
	
	private double[] sphereDim = new double[4];
	private double[] profile;
	/** User's chosen start point */
	private double[] initialPoint = new double[3];
	/** Fill with random edge points */
	private double[][] edgePoints;
	
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
		
		final int w = imp.getWidth();
		final int h = imp.getHeight();
		final int d = imp.getStackSize();
		
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
		
		
		this.profile = new ProfilePlot(imp).getProfile();
		
//		ProfilePlot pp = new ProfilePlot();
//		Profiler p = new Profiler();
//		this.profile = pp.getProfile();
//		Plot profilePlot = new Plot("Profile plot", "distance", "value", xVals, yVals);
//		pp.getStraightLineProfile(roi, cal, ip);
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
		this.initialPoint = a;;
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
