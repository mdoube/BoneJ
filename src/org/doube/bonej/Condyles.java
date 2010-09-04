package org.doube.bonej;

import org.doube.geometry.Centroid;
import org.doube.geometry.FitEllipsoid;
import org.doube.util.ImageCheck;
import org.doube.util.RoiMan;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.gui.WaitForUserDialog;
import ij.plugin.PlugIn;
import ij.plugin.frame.RoiManager;

/**
 * Find condyle properties using FitEllipsoid. Takes a mean of a user-determined 
 * number of manual inputs.
 * 
 * @author Nick Powell
 *
 */
public class Condyles implements PlugIn {
	
	private ImagePlus imp;
	
	/** Number of ellipsoids to fit (manually), from which to take the mean dimensions */
	private int numEllipsoids = 3;
	
	/** The mean cartesian coordinates x,y,z of the centroids of all the ellipsoids found */
	double[] meanCentroid;
	
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
		Object[][] ellipsoids = getEllipsoids(imp, numEllipsoids);
		meanCentroid = findMeans(ellipsoids);
		
		IJ.log("Mean centroid x,y,z: " + meanCentroid[0] + "; " + meanCentroid[1] + "; " + meanCentroid[2] + "; ");
		
		return;
	}
	
	public void manualOptions() {
		
		GenericDialog gd = new GenericDialog("Manual condyle fit options");
		gd.addNumericField("Take mean of", this.numEllipsoids, 0, 3, "ellipsoids");
		gd.showDialog();
		if(gd.wasCanceled()) { return; }
		
		this.numEllipsoids = (int) gd.getNextNumber();
		
		return;
	}
	
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
	
	public double[] findMeans(Object[][] ellipsoids) {
		
		double[][] ellipseCentroids = new double[ellipsoids.length][3];
		double[][] ellipseRadii = new double[ellipsoids.length][3];
		
		for(int num = 0; num < ellipsoids.length; num++) {
			
			ellipseCentroids[num] = (double[]) ellipsoids[num][0];
			ellipseRadii[num] = (double[]) ellipsoids[num][1];
			
		}
		
		double[] centroid = Centroid.getCentroid(ellipseCentroids);
		
		return centroid;
	}
	
	public double[] getMeanCentroid() {
		return this.meanCentroid;
	}
	
	private void showAxes3D() {
		
		
		
		return;
	}

}
