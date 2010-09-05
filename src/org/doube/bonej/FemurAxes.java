package org.doube.bonej;

import java.util.ArrayList;
import java.util.List;

import javax.vecmath.Color3f;
import javax.vecmath.Point3f;

import org.doube.util.ImageCheck;
import org.doube.util.ThresholdGuesser;

import customnode.CustomPointMesh;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.gui.ImageCanvas;
import ij.gui.ImageWindow;
import ij.gui.WaitForUserDialog;
import ij.measure.Calibration;
import ij.plugin.Duplicator;
import ij.plugin.PlugIn;
import ij.process.StackConverter;
import ij3d.Content;
import ij3d.Image3DUniverse;

/**
 * 
 * 
 * @author Nick Powell
 *
 */
public class FemurAxes implements PlugIn {
	
	private ImagePlus imp;
	
	private int[] shaftEndSlices;
	
	private double[] shaftCentroid;
	private double[][] shaftVector;
	
	/** Image thresholds */
	private double[] thresholds;
	private double min, max;
	
	/** Pixel dimensions in 'unit's */
	private double vD;
	/** Image dimensions in 'unit's */
	private int h, w;
	
	public void run(String arg) {
		
		if (!ImageCheck.checkEnvironment())
			return;
		this.imp = IJ.getImage();
		if (null == imp) {
			IJ.noImage();
			return;
		}
		
		Calibration cal = imp.getCalibration();
		vD = cal.pixelDepth;
		w = imp.getWidth();
		h = imp.getHeight();
		
		thresholds = ThresholdGuesser.setDefaultThreshold(imp);
		min = thresholds[0];
		max = thresholds[1];
		boolean isHUCalibrated = ImageCheck.huCalibrated(imp);
		if (isHUCalibrated) {
			min = cal.getRawValue(min);
			max = cal.getRawValue(max);
		}
		IJ.log("min: " +min+ "; max: " +max+ "");
		
		/* Shaft with visual confirmation */
		int i = -1;
		while(i < 0) {
			shaftEndSlices = getShaftEndSlices();
			shaftCentroid = getShaftCentroid(imp, shaftEndSlices, min, max);
			
			double[] proxCentroid = getSliceCentroid(imp, min, max, shaftEndSlices[0]);
			double[] distCentroid = getSliceCentroid(imp, min, max, shaftEndSlices[1]);
			
			IJ.log("prox: " + proxCentroid[2] + "");
			IJ.log("dist: " + distCentroid[2] + "");
			
			showShaftEnds(imp, proxCentroid, distCentroid);
			
			new WaitForUserDialog("Check shaft ends are reasonable, then hit \'OK\'").show();
			GenericDialog gdOK = new GenericDialog("OK?");
			gdOK.addCheckbox("These shaft ends look good", true);
			gdOK.showDialog();
			if(gdOK.getNextBoolean()) {
				i = 1;
			}
		}
		
		
		/* Condyles */
		Condyles cons = new Condyles();
//		cons.manualOptions();	// set numEllipsoids manually
		Object[] properties1 = cons.findProperties(cons.getEllipsoids(imp, 3));
		cons.showResults(properties1).show("Results");
		Object[] properties2 = cons.findProperties(cons.getEllipsoids(imp, 3));
		cons.showResults(properties2).show("Results");
		
		
		shaftVector = BicondylarAngle.regression3D(imp, shaftCentroid, shaftEndSlices[0], shaftEndSlices[1], min, max);
		if (this.shaftVector == null) { return; }
		
		
		
		
		
		return;
	}
	
//	public void manualOptions() {
//		
//		GenericDialog gd = new GenericDialog("Manual femur axes options");
//		gd.addNumericField("Take mean of", this.numEllipsoids, 0, 3, "ellipsoids");
//		gd.showDialog();
//		if(gd.wasCanceled()) { return; }
//		
//		this.min = gd.getNextNumber();
//		this.max = gd.getNextNumber();
//		
//		return;
//	}
	
	/**
	 * Use ShaftGuesser to get the end slices of the shaft.
	 * 
	 * @return shaftPosition
	 * 				int[2] containing proximal and distal slices of the shaft.
	 */
	public int[] getShaftEndSlices() {
		
		ShaftGuesser sg = new ShaftGuesser();
		sg.run(null);
		int[] shaftPosition = sg.getShaftPosition();
		
		return shaftPosition;
	}
	
	/**
	 * Use SliceGeometry to get the centroid of the desired slice
	 * 
	 * @param imp
	 * @param min
	 * @param max
	 * @param shaftSlice
	 * @return
	 */
	public double[] getSliceCentroid(ImagePlus imp, double min, double max, int shaftSlice) {
		
		SliceGeometry sgeo = new SliceGeometry();
		sgeo.setParameters(imp, shaftSlice, shaftSlice);
		sgeo.calculateCentroids(imp, min, max);
		double[][] sliceCentroids = sgeo.getSliceCentroids();
		
		double[] sliceCentroid = new double[3];
		sliceCentroid[0] = sliceCentroids[0][shaftSlice];
		sliceCentroid[1] = sliceCentroids[1][shaftSlice];
		sliceCentroid[2] = shaftSlice * this.vD;
		
		return sliceCentroid;
	}
	
	public double[] getShaftCentroid(ImagePlus imp, int[] shaftEndSlices, double min, double max) {
		
		Moments m = new Moments();
		final double[] centroid = m.getCentroid3D(imp, shaftEndSlices[0], shaftEndSlices[1], min, max, 0, 1);
		if (centroid[0] < 0) {
			IJ.error("Empty Stack", "No voxels available for calculation."
					+ "\nCheck your ROI and threshold.");
			return null;
		}
		
		return centroid;
	}
	
	private void showShaftEnds(ImagePlus imp, double[] a, double[] b) {
		
		ImagePlus line3Dimp = new Duplicator().run(imp, 1, imp.getImageStackSize());
		List<Point3f> line3Da = new ArrayList<Point3f>();
		List<Point3f> line3Db = new ArrayList<Point3f>();
		/* initialise and show the 3D universe */
		Image3DUniverse univ = new Image3DUniverse();
		univ.show();
		
		/* Draw cross confirming location of proximal slice centroid */
		Point3f start1 = new Point3f();
		start1.x = (float) (a[0] + (w * 0.2));
		start1.y = (float) a[1];
		start1.z = (float) a[2];
		line3Da.add(start1);

		Point3f end1 = new Point3f();
		end1.x = (float) (a[0] - (w * 0.2));
		end1.y = (float) a[1];
		end1.z = (float) a[2];
		line3Da.add(end1);
		
		Point3f start2 = new Point3f();
		start2.x = (float) a[0];
		start2.y = (float) (a[1] + (h * 0.2));
		start2.z = (float) a[2];
		line3Da.add(start2);

		Point3f end2 = new Point3f();
		end2.x = (float) a[0];
		end2.y = (float) (a[1] - (h * 0.2));
		end2.z = (float) a[2];
		line3Da.add(end2);
		
		/* Draw cross confirming location of distal slice centroid */
		Point3f start3 = new Point3f();
		start3.x = (float) (b[0] + (w * 0.2));
		start3.y = (float) b[1];
		start3.z = (float) b[2];
		line3Db.add(start3);

		Point3f end3 = new Point3f();
		end3.x = (float) (b[0] - (w * 0.2));
		end3.y = (float) b[1];
		end3.z = (float) b[2];
		line3Db.add(end3);
		
		Point3f start4 = new Point3f();
		start4.x = (float) b[0];
		start4.y = (float) (b[1] + (h * 0.2));
		start4.z = (float) b[2];
		line3Db.add(start4);

		Point3f end4 = new Point3f();
		end4.x = (float) b[0];
		end4.y = (float) (b[1] - (h * 0.2));
		end4.z = (float) b[2];
		line3Db.add(end4);
		
		float red = 0.0f;
		float green = 0.5f;
		float blue = 1.0f;
		Color3f aColour = new Color3f(red, green, blue);
		/* show the line */
		try {
			univ.addLineMesh(line3Da, aColour, "Proximal slice", false).setLocked(true);
			univ.addLineMesh(line3Db, aColour, "Distal slice", false).setLocked(true);
			new StackConverter(line3Dimp).convertToGray8();
			Content c = univ.addVoltex(line3Dimp);
			c.setLocked(true);
		} catch (NullPointerException npe) {
			IJ.log("3D Viewer was closed before rendering completed.");
			return;
		}
		
		return;
	}
	
//	public Object[] 

}
