package org.doube.bonej;

import java.util.ArrayList;
import java.util.List;

import javax.vecmath.Color3f;
import javax.vecmath.Point3f;

import org.doube.geometry.Centroid;
import org.doube.geometry.FitEllipsoid;
import org.doube.geometry.Trig;
import org.doube.util.ImageCheck;
import org.doube.util.RoiMan;

import customnode.CustomPointMesh;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.gui.WaitForUserDialog;
import ij.measure.ResultsTable;
import ij.plugin.Duplicator;
import ij.plugin.PlugIn;
import ij.plugin.frame.RoiManager;
import ij.process.StackConverter;
import ij3d.Content;
import ij3d.Image3DUniverse;

/**
 * Find condyle properties using FitEllipsoid. Takes a mean of a user-determined 
 * number of manual inputs.
 * 
 * @author Nick Powell
 *
 */
public class Condyles implements PlugIn {
	
	private ImagePlus imp;
	
	/** If we have the centroids already from a previous run, ask for them. */
	public boolean alreadyHaveCentroids = false;
	
	/** Number of ellipsoids to fit (manually), from which to take the mean dimensions */
	private int numEllipsoids = 3;
	
	/** The mean cartesian coordinates x,y,z of the centroids of all the ellipsoids found */
	double[] meanCentroid;
	
	/** The mean radii of all the ellipsoids found */
	double[] meanRadii;
	
	/** The volume of the mean ellipsoid, in 'unit's */
	double meanVolume;
	
	/** The standard deviations of the radii */
	double[] sdRadii = new double[3];
	
	/** The unit vector between the two condyles, as well as the distance */
	double[] iCV;
	
	/** The x,y,z coordinates of the midpoint between the centroids of the condyles */
	double[] midPoint;
	
	public void run(String arg) {
		
		/* Modified from SphereFitter */
		if (!ImageCheck.checkEnvironment())
			return;
		this.imp = IJ.getImage();
		if (null == imp) {
			IJ.noImage();
			return;
		}
		
		manualOptions();
//		if(this.alreadyHaveCentroids) {
//			double[]
//			
//			GenericDialog gd = new GenericDialog("Get manual centroids");
//			for(int i = 0; i < numEllipsoids; i++) {
//				gd.addNumericField("Centroid 1 x: ", 0.0, 14, 17, " ");
//			}
//			gd.showDialog();
//			if(gd.wasCanceled()) {return;}
//			for(int i = 0; i < numEllipsoids; i++) {
//				
//			}
//		}
//		
		Object[][] ellipsoids1 = getEllipsoids(imp, numEllipsoids);
		Object[] properties1 = findProperties(ellipsoids1);
		
		meanCentroid = (double[]) properties1[0];
		meanRadii = (double[]) properties1[1];
		meanVolume = 4/3 * meanRadii[0] * meanRadii[1] * meanRadii[2] * Math.PI;
		
		sdRadii[0] = Math.sqrt(meanRadii[0]);
		sdRadii[1] = Math.sqrt(meanRadii[1]);
		sdRadii[2] = Math.sqrt(meanRadii[2]);
		
//		IJ.log("Mean centroid x,y,z: " + meanCentroid[0] + "; " + meanCentroid[1] + "; " + meanCentroid[2] + "; ");
		
		Object[][] ellipsoids2 = getEllipsoids(imp, numEllipsoids);
		Object[] properties2 = findProperties(ellipsoids2);
		
		showResults(properties1).show("Results");
		showResults(properties2).show("Results");
		
		midPoint = getMidPoint((double[]) properties1[0], (double[]) properties2[0]);
		
		iCV = getUnitVector((double[]) properties1[0], (double[]) properties2[0]);
		
		annotate3D(imp, (double[]) properties1[0], (double[]) properties2[0], midPoint);
		
		return;
	}
	
	public void manualOptions() {
		
		GenericDialog gd = new GenericDialog("Manual condyle fit options");
		gd.addNumericField("Take mean of", this.numEllipsoids, 0, 3, "ellipsoids");
//		gd.addCheckbox("Already have centroids?", this.alreadyHaveCentroids);
		gd.showDialog();
		if(gd.wasCanceled()) { return; }
		
		this.numEllipsoids = (int) gd.getNextNumber();
//		this.alreadyHaveCentroids = gd.getNextBoolean();
		
		return;
	}
	
	/**
	 * Run getEllipsoids after setting numEllipsoids manually.
	 * 
	 * @param imp
	 * @return
	 */
	public Object[][] getEllipsoids(ImagePlus imp) {
		Object[][] ellipsoids = getEllipsoids(imp, this.numEllipsoids);
		return ellipsoids;
	}
	
	/**
	 * getEllipsoids without manual input.
	 * 
	 * @param imp
	 * @param numEllipsoids
	 * @return
	 */
	public Object[][] getEllipsoids(ImagePlus imp, int numEllipsoids) {
		
		ImageCheck ic = new ImageCheck();
		if (!ic.isVoxelIsotropic(imp, 0.05)) {
			if (!IJ.showMessageWithCancel("Voxel depth problem",
					"Voxels are anisotropic." + "\nWidth = "
							+ imp.getCalibration().pixelWidth + "\nHeight = "
							+ imp.getCalibration().pixelHeight + "\nDepth = "
							+ imp.getCalibration().pixelDepth
							+ "\nClick OK if voxel dimensions are correct.")) {
				imp.unlock();
				IJ.run("Properties...");
			}
		}
		if (!imp.lock()) {
			imp.lock();			// relock if unlocked to reset properties
		}
		
		Object[][] ellipsoids = new Object[numEllipsoids][5];
//		Object[] ellipsoid = new Object[5];
		
		for(int i = 0; i < numEllipsoids; i++) {
			
			IJ.run("ROI Manager...");
			
			/* Manual clicking */
			new WaitForUserDialog("Fill the ROI Manager with >9 points,\n" 
										+ "then hit \'OK\'").show();
			
			RoiManager roiMan = RoiManager.getInstance();
			double[][] points = RoiMan.getRoiManPoints(imp, roiMan);
			
//			int p = -1;
//			while (p < 0) {
//				if(points.length < 9) {
//					p = 1;
//				}
//				else {
//					new WaitForUserDialog("Fill the ROI Manager with >9 points,\n" 
//							+ "then hit \'OK\'").show();
//				}
//			}
			
			try {
				ellipsoids[i] = FitEllipsoid.yuryPetrov(points);
//				ellipsoid = FitEllipsoid.yuryPetrov(points);
			} catch (IllegalArgumentException ia) {
				IJ.showMessage(ia.getMessage());
				continue;
			} catch (RuntimeException re) {
				IJ.showMessage("Can't fit ellipsoid to points.\n"
								+ "Add more points and try again.");
				continue;
			}
			
			roiMan.close();
			
			double[] thisCentroid = (double[]) ellipsoids[i][0];
			IJ.log("Centroid x,y,z: " + thisCentroid[0] + "; " + thisCentroid[1] + "; " + thisCentroid[2] + "; ");
			
//			ellipsoids[i] = ellipsoid;
		}
		
		return ellipsoids;
	}
	
	/**
	 * Find the mean centroid and radii of several ellipsoids, as well as the 
	 * standard deviations on the radii.
	 * 
	 * @param ellipsoids
	 * @return Object[] properties, containing:
	 * 		properties[0] = centroid[x,y,z];
			properties[1] = radii[a,b,c];
			properties[2] = sdR[a,b,c];
	 */
	public Object[] findProperties(Object[][] ellipsoids) {
		
		double[][] ellipseCentroids = new double[ellipsoids.length][3];
		double[][] ellipseRadii = new double[ellipsoids.length][3];
		
		for(int num = 0; num < ellipsoids.length; num++) {
			
			ellipseCentroids[num] = (double[]) ellipsoids[num][0];
			ellipseRadii[num] = (double[]) ellipsoids[num][1];
			
		}
		
		double[] centroid = Centroid.getCentroid(ellipseCentroids);
		double[] radii = Centroid.getCentroid(ellipseRadii);
		double[] sdR = ShaftGuesser.variance(ellipseRadii);
		
		Object[] properties = new Object[3];
		properties[0] = centroid;
		properties[1] = radii;
		properties[2] = sdR;
		
		return properties;
	}
	
	/**
	 * Returns the unit vector of the line connecting two points in 3D space, as 
	 * well as the distance between them.
	 * 
	 * (Was: Get the unit vector of the line connecting the two condyle centroids.)
	 * 
	 * @param point1
	 * 				double[3]
	 * @param point2
	 * 				double[3]
	 * @return double[4] containing unitVector (x,y,z) followed by the distance.
	 */
	public double[] getUnitVector(double[] point1, double[] point2) {
		
		double distance = Trig.distance3D(point1, point2);
		double[] unitVector = new double[4];
		unitVector[0] = (point1[0] - point2[0]) / distance;
		unitVector[1] = (point1[1] - point2[1]) / distance;
		unitVector[2] = (point1[2] - point2[2]) / distance;
		
		unitVector[3] = distance;
		
		return unitVector;
	}
	
	/**
	 * Get the centroid of (mid-point between) two points [x,y,z].
	 * 
	 * @param meanCentroid1[x,y,z]
	 * @param meanCentroid2[x,y,z]
	 * @return midPoint[x,y,z]
	 */
	public double[] getMidPoint(double[] meanCentroid1, double[] meanCentroid2) {
		
		double[][] centroids = new double[2][3];
		centroids[0] = meanCentroid1;
		centroids[1] = meanCentroid2;
		double[] midPoint = Centroid.getCentroid(centroids);
		
		return midPoint;
	}
	
	/**
	 * Show the centroid, radii, condyle volume and the standard deviations 
	 * on each of the radii.
	 * 
	 * @param properties
	 * @return
	 */
	public ResultsTable showResults(Object[] properties) {
		
		meanCentroid = (double[]) properties[0];
		meanRadii = (double[]) properties[1];
		meanVolume = 4/3 * meanRadii[0] * meanRadii[1] * meanRadii[2] * Math.PI;
		
		sdRadii[0] = Math.sqrt(meanRadii[0]);
		sdRadii[1] = Math.sqrt(meanRadii[1]);
		sdRadii[2] = Math.sqrt(meanRadii[2]);
		
		ResultsTable rt = ResultsTable.getResultsTable();
		rt.incrementCounter();
		rt.addValue("Centroid_x", meanCentroid[0]);
		rt.addValue("Centroid_y", meanCentroid[1]);
		rt.addValue("Centroid_z", meanCentroid[2]);
		rt.addValue("Radius_a_(mm)", meanRadii[0]);
		rt.addValue("Radius_b_(mm)", meanRadii[1]);
		rt.addValue("Radius_c_(mm)", meanRadii[2]);
		rt.addValue("SD rad_a_(mm)", sdRadii[0]);
		rt.addValue("SD rad_b_(mm)", sdRadii[1]);
		rt.addValue("SD rad_c_(mm)", sdRadii[2]);
		rt.addValue("Volume_(mm^3)", meanVolume);
		
		return rt;
	}
	
	public double[] getMeanCentroid() {
		return this.meanCentroid;
	}
	public double[] getMeanRadii() {
		return this.meanRadii;
	}
	
	/**
	 * Draw the centroids of each ellipse and a line connecting them. Also show the 
	 * mid-point of this line.
	 * 
	 * @param imp
	 * @param a
	 * @param b
	 */
	public void annotate3D(ImagePlus imp, double[] a, double[] b, double[] midPoint) {
		
		ImagePlus con3Dimp = new Duplicator().run(imp, 1, imp.getImageStackSize());
		List<Point3f> line = new ArrayList<Point3f>();
		List<Point3f> centrePoints = new ArrayList<Point3f>();
		
		/* initialise and show the 3D universe */
		Image3DUniverse univ = new Image3DUniverse();
		univ.show();
		
		/* Draw line connecting the two points */
		Point3f start1 = new Point3f();
		start1.x = (float) a[0];
		start1.y = (float) a[1];
		start1.z = (float) a[2];
		line.add(start1);

		Point3f end1 = new Point3f();
		end1.x = (float) b[0];
		end1.y = (float) b[1];
		end1.z = (float) b[2];
		line.add(end1);
		
		/* Draw centroids */
		Point3f cent1 = new Point3f();
		cent1.x = (float) a[0];
		cent1.y = (float) a[1];
		cent1.z = (float) a[2];
		centrePoints.add(cent1);
		
		Point3f cent2 = new Point3f();
		cent2.x = (float) b[0];
		cent2.y = (float) b[1];
		cent2.z = (float) b[2];
		centrePoints.add(cent2);
		
		Point3f mid = new Point3f();
		mid.x = (float) midPoint[0];
		mid.y = (float) midPoint[1];
		mid.z = (float) midPoint[2];
		centrePoints.add(mid);
		
		/* Line properties */
		float red1 = 0.0f;
		float green1 = 0.5f;
		float blue1 = 1.0f;
		Color3f Colour1 = new Color3f(red1, green1, blue1);
		
		/* Point properties */
		CustomPointMesh points = new CustomPointMesh(centrePoints);
		points.setPointSize(5.0f);
		Color3f Colour2 = new Color3f(0.0f, 1.0f, 0.5f);
		points.setColor(Colour2);
		
		try {
			univ.addLineMesh(line, Colour1, "Inter condylar axis", false).setLocked(true);
			univ.addCustomMesh(points, "Centroid").setLocked(true);
			
			new StackConverter(con3Dimp).convertToGray8();
			Content c = univ.addVoltex(con3Dimp);
			c.setLocked(true);
		} catch (NullPointerException npe) {
			IJ.log("3D Viewer was closed before rendering completed.");
			return;
		}
		
		return;
	}

}
