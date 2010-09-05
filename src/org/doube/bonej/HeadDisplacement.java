package org.doube.bonej;

import java.awt.Color;
import java.util.ArrayList;
import java.util.List;

import javax.vecmath.Color3f;
import javax.vecmath.Point3f;

import org.doube.geometry.Trig;
import org.doube.util.ImageCheck;
import org.doube.util.ThresholdGuesser;

import customnode.CustomPointMesh;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.ImageCanvas;
import ij.gui.ImageWindow;
import ij.measure.Calibration;
import ij.measure.ResultsTable;
import ij.plugin.Duplicator;
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;
import ij.process.StackConverter;
import ij3d.Content;
import ij3d.Image3DUniverse;

/**
 * <p>
 * Head Displacement of Femur <br />
 * Find the cartesian displacement of the head centre from the centroid of the 
 * proximal diaphysis. </p>
 * 
 * <p>
 * 
 * </p>
 * 
 * @author Nick Powell
 *
 */
public class HeadDisplacement implements PlugIn {
	
	/** Contains the absolute x,y,z distance between the centre of the femoral head and the centroid of the chosen slice of the shaft */
	private double[] displacement;
	private double[] proximalShaftCentroid;
	/** Image thresholds */
	private double[] thresholds;
	private double min, max;
	/** List of slice centroids */
	private double[][] sliceCentroids;
	/** Pixel dimensions in 'unit's */
	private double vD;
	/** Image dimensions in 'unit's */
	private int h, w;
	
	public void run(String arg) {
		
		/* Copied from SphereFitter */
		if (!ImageCheck.checkEnvironment())
			return;
		final ImagePlus imp = IJ.getImage();
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
		
		/* Get the proximal slice of the diaphysis */
		ShaftGuesser sg = new ShaftGuesser();
		sg.run(arg);
		int proximalShaftSlice = sg.getShaftPosition()[0];
		
//		IJ.log("proximalShaftSlice (slice #): " + proximalShaftSlice + "; ");
		
		/* How to use SphereEdgeGuesser */
//		ImageCanvas canvas = imp.getWindow().getCanvas();
		
		SphereEdgeGuesser seg = new SphereEdgeGuesser();
		seg.getManualSettings();
		seg.getInitialPointUser(imp);
		seg.findSphere(imp, seg.runNum, seg.numVectors, seg.sd2DMult, seg.sd3DMult, seg.fitEllipsoid, seg.ignoreCentroidDirection);
		double[] meanDimensions = seg.getMeanDimensions();
		
		/* Just to prove we've got the meanDimensions */
		seg.annotateCentre(imp, meanDimensions).show();
//		IJ.log("head centre (mm): " + meanDimensions[0] + "; " + meanDimensions[1] + "; " + meanDimensions[2] + "; ");
		
		displacement = headDisplacement(imp, meanDimensions, proximalShaftSlice);
		
//		IJ.log("displacement (x,y,z,distance) (mm): " + displacement[0] + "; " + displacement[1] + "; " + displacement[2] + "; dist: " + displacement[3] + "; ");
		
		annotateImage(imp, meanDimensions, proximalShaftCentroid);
		
		ResultsTable rt = ResultsTable.getResultsTable();
		rt.incrementCounter();
		rt.addLabel("Type", "displacement");
		rt.addValue("x (mm)", displacement[0]);
		rt.addValue("y (mm)", displacement[1]);
		rt.addValue("z (mm)", displacement[2]);
		rt.addValue("distance (mm)", displacement[3]);
		rt.show("Results");
		
		return;
	}
	
	
	/**
	 * Find the displacement of the centre of the femoral head, in 'unit's.
	 * 
	 * @param imp
	 * @param headCentre
	 * @param proximalShaftSlice
	 * @return
	 */
	private double[] headDisplacement(ImagePlus imp, double[] headCentre, int proximalShaftSlice) {
		
		/* Get the centroid of the proximal shaft slice */
		SliceGeometry sgeo = new SliceGeometry();
		sgeo.setParameters(imp, proximalShaftSlice, proximalShaftSlice);
		sgeo.calculateCentroids(imp, min, max);
		this.sliceCentroids = sgeo.getSliceCentroids();
		
		this.proximalShaftCentroid = new double[3];
		proximalShaftCentroid[0] = sliceCentroids[0][proximalShaftSlice];
		proximalShaftCentroid[1] = sliceCentroids[1][proximalShaftSlice];
		proximalShaftCentroid[2] = proximalShaftSlice * this.vD;
		
		/* Confirm we've got the centroid */
//		IJ.log("Slice Centroid x, y (mm): " + "; " + sliceCentroids[0][proximalShaftSlice] + "; " + sliceCentroids[1][proximalShaftSlice] + "; ");
		
		displacement = new double[4];
		
		/* Displacement in mm */
		displacement[0] = Math.abs(proximalShaftCentroid[0] - headCentre[0]);
		displacement[1] = Math.abs(proximalShaftCentroid[1] - headCentre[1]);
		displacement[2] = Math.abs(proximalShaftCentroid[2] - headCentre[2]);
		
		displacement[3] = Trig.distance3D(headCentre, proximalShaftCentroid);
		
		return displacement;
	}
	
	/**
	 * Draw the connecting line between the proximal shaft slice centroid and the 
	 * femoral head centre.
	 * 
	 * @param imp
	 * @param a
	 * @param b
	 */
	private void annotateImage(ImagePlus imp, double[] a, double[] b) {
		
		ImagePlus line3Dimp = new Duplicator().run(imp, 1, imp.getImageStackSize());
		List<Point3f> line3D = new ArrayList<Point3f>();
		/* initialise and show the 3D universe */
		Image3DUniverse univ = new Image3DUniverse();
		univ.show();
			
		/* Draw line connecting head to proximal slice centroid */
		Point3f start1 = new Point3f();
		start1.x = (float) a[0];
		start1.y = (float) a[1];
		start1.z = (float) a[2];
		line3D.add(start1);

		Point3f end1 = new Point3f();
		end1.x = (float) b[0];
		end1.y = (float) b[1];
		end1.z = (float) b[2];
		line3D.add(end1);
		
		/* Draw cross confirming location of proximal slice centroid */
		Point3f start2 = new Point3f();
		start2.x = (float) (b[0] + (w * 0.2));
		start2.y = (float) b[1];
		start2.z = (float) b[2];
		line3D.add(start2);

		Point3f end2 = new Point3f();
		end2.x = (float) (b[0] - (w * 0.2));
		end2.y = (float) b[1];
		end2.z = (float) b[2];
		line3D.add(end2);
		
		Point3f start3 = new Point3f();
		start3.x = (float) b[0];
		start3.y = (float) (b[1] + (h * 0.2));
		start3.z = (float) b[2];
		line3D.add(start3);

		Point3f end3 = new Point3f();
		end3.x = (float) b[0];
		end3.y = (float) (b[1] - (h * 0.2));
		end3.z = (float) b[2];
		line3D.add(end3);
		
		float red = 0.0f;
		float green = 0.5f;
		float blue = 1.0f;
		Color3f aColour = new Color3f(red, green, blue);
		/* show the line */
		try {
			univ.addLineMesh(line3D, aColour, "Principal axes", false).setLocked(true);
			new StackConverter(line3Dimp).convertToGray8();
			Content c = univ.addVoltex(line3Dimp);
			c.setLocked(true);
		} catch (NullPointerException npe) {
			IJ.log("3D Viewer was closed before rendering completed.");
			return;
		}
		
		return;
	}

}
