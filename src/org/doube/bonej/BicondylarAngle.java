package org.doube.bonej;

import java.awt.Rectangle;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.measure.Calibration;
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;

import org.doube.geometry.Trig;
import org.doube.geometry.Vectors;
import org.doube.jama.Matrix;
import org.doube.jama.SingularValueDecomposition;
import org.doube.util.ImageCheck;
import org.doube.util.ResultInserter;
import org.doube.util.ThresholdGuesser;

/**
 * <p>
 * Bicondylar Angle<br />
 * Tool to calculate the bicondylar angle of 3D images of femora.
 * </p>
 * 
 * <p>
 * 
 * Bicondylar angle is the angle formed at the intersection between the coplanar
 * lines <b>N</b> and <b>S</b>. <br>
 * 
 * 
 * 
 * 
 * <b>S</b> is the singular value decomposition (orthogonal distance regression)
 * vector that passes through the centroid (<b>B</b>) of the bone and describes
 * the long axis of the bone.<br/>
 * <b>C</b> is the centre of a sphere fit to the femoral head.<br/>
 * 
 * <b>P</b> is the plane that contains <b>S</b> and <b>C</b>. <br />
 * <b>N</b> is the projection onto <b>P</b> of a vector originating at <b>C</b>
 * and passing through the 'middle' of the femoral neck.
 * </p>
 * 
 *@author Nick Powell
 */
public class BicondylarAngle implements PlugIn {
	
	/** For HU units (copied from Neck Shaft Angle) */
	private boolean fieldUpdated = false;
	
	/** Pixel dimensions in 'unit's */
	private double vH, vW, vD;
	/** Linear unit of measure */
	private String units;
	/** Working in units */
	private boolean inUnits;
	
	/** First and last slices of the femoral shaft (inclusive) */
	private int startSlice, endSlice;
	private int[] shaftPosition = new int[2];
	
	private double[] medCentre, latCentre;
	private double[] centroid;
	private double[][] shaftVector;
	
	
	public void run(String arg) {

		if (!ImageCheck.checkEnvironment())
			return;
		final ImagePlus imp = IJ.getImage();
		if (null == imp) {
			IJ.noImage();
			return;
		}
		
		Calibration cal = imp.getCalibration();
		
		vW = cal.pixelWidth;
		vH = cal.pixelHeight;
		vD = cal.pixelDepth;
		units = cal.getUnits();
		
		double[] thresholds = ThresholdGuesser.setDefaultThreshold(imp);
		double min = thresholds[0];
		double max = thresholds[1];
		String pixUnits;
		if (ImageCheck.huCalibrated(imp)) {
			pixUnits = "HU";
			fieldUpdated = true;
		} else
			pixUnits = "grey";
		cal = imp.getCalibration();
		
		/* Based on Neck Shaft Angle */
		GenericDialog gd = new GenericDialog("Setup");
		gd.addNumericField("Shaft Start Slice:", 1, 0);
		gd.addNumericField("Shaft End Slice:", imp.getImageStackSize(), 0);
		gd.addCheckbox("HU Calibrated", ImageCheck.huCalibrated(imp));
		gd.addNumericField("Bone Min:", min, 1, 6, pixUnits + " ");
		gd.addNumericField("Bone Max:", max, 1, 6, pixUnits + " ");
		gd
				.addMessage("Only pixels >= bone min\n"
						+ "and <= bone max are used.");
		gd.addCheckbox("Calculate curvature", true);
		gd.addCheckbox("Show annotations in 3D", true);
//		gd.addDialogListener(this);
		gd.showDialog();
		if (gd.wasCanceled()) {
			return;
		}
		final int startSlice = (int) gd.getNextNumber();
		final int endSlice = (int) gd.getNextNumber();
		boolean isHUCalibrated = gd.getNextBoolean();
		min = gd.getNextNumber();
		max = gd.getNextNumber();
		if (isHUCalibrated) {
			min = cal.getRawValue(min);
			max = cal.getRawValue(max);
		}
		
		
		
		
		ShaftGuesser sg = new ShaftGuesser();
		sg.run(arg);
		this.shaftPosition = sg.getShaftPosition();
		this.startSlice = shaftPosition[0];
		this.endSlice = shaftPosition[1];
		
		/* Centroid and regression vector of the shaft - based on NeckShaftAngle */
		Moments m = new Moments();
		final double[] centroid = m.getCentroid3D(imp, startSlice, endSlice, min, max, 0, 1);
		this.centroid = centroid;
		if (centroid[0] < 0) {
			IJ.error("Empty Stack", "No voxels available for calculation."
					+ "\nCheck your ROI and threshold.");
			return;
		}
		
		this.shaftVector = regression3D(imp, centroid, startSlice, endSlice, min, max);
		if (this.shaftVector == null) {
			return;
		}
		
		
	}
	
	/**
	 * Gives unit vector of the line joining the centroids of the condyles.
	 * Written with the centres of the condyles in mind; will work for the 
	 * distal points.
	 * 
	 * Based on Neck_Shaft_Angle.neckVector()
	 * 
	 * @param medCentre (x,y,z) location of the centroid of the ellipsoid which approximates the medial condyle.
	 * @param latCentre (x,y,z) location of the centroid of the ellipsoid which approximates the lateral condyle.
	 * @return
	 */
	private double[][] bicondylarVector(double[] medCentre, double[] latCentre) {
		// have to calculate d to make sure that bicondylarVector is a unit vector
		double d = Trig.distance3D(medCentre, latCentre);

		double[][] bicondylarVector = new double[3][1];
		bicondylarVector[0][0] = (medCentre[0] - latCentre[0]) / d;
		bicondylarVector[1][0] = (medCentre[1] - latCentre[1]) / d;
		bicondylarVector[2][0] = (medCentre[2] - latCentre[2]) / d;
		return bicondylarVector;
	}
	
	/**
	 * Gives the distance of the line joining the centroids of the condyles, 
	 * in 'unit's (use Trig.distance3D(medCentre, latCentre) if already in 'unit's).
	 * 
	 * @param medCentre
	 * @param latCentre
	 * @return
	 */
	private double interCondylarDistance(double[][] bicondylarVector) {
		
//		if(this.inUnits) {
//			double interCondylarDistance = Trig.distance3D(medCentre, latCentre);
//		}
		
		double uX = bicondylarVector[0][0] * vW;
		double uY = bicondylarVector[1][0] * vH;
		double uZ = bicondylarVector[2][0] * vD;
			
		double interCondylarDistance = Trig.distance3D(uX * vW, uY * vH, uZ * vD);
		
		return interCondylarDistance;
	}
	
	/**
	 * Gives the distance in 'unit's between the most medial and most lateral points 
	 * of the condyles.
	 * 
	 * @param medCentre
	 * @param latCentre
	 * @param medialCondyle
	 * @param lateralCondyle
	 * @return
	 */
	private double epiCondylarDistance(double[] medCentre, double[] latCentre, Object[] medialCondyle, Object[] lateralCondyle) {
		
		double[] medRadii = (double[]) medialCondyle[1];
		double[] latRadii = (double[]) lateralCondyle[1];
		
		double epiCondylarDistance = Trig.distance3D(medCentre, latCentre) + (medRadii[2] / 2) + (latRadii[2] / 2);
		
		return epiCondylarDistance;
	}
	
	
	/* Copied directly from Neck_Shaft_Angle */
//	private void calculateAngles(ImagePlus imp, double[] latCentre) {
//		double[][] bicondylarVector = bicondylarVector(medCentre, latCentre);
//		double[][] projectionPlane = getProjectionPlane(shaftVector, medCentre, this.centroid);
//		double[][] neckPlane = neckPlane(bicondylarVector, projectionPlane);
//		double[][] testVector = testVector(projectionPlane, neckPlane);
//		// P . Q = ||P|| ||Q|| cos(a) so if P and Q are unit vectors, then P.Q =
//		// cos(a)
//		Matrix PP = new Matrix(projectionPlane);
//		PP.printToIJLog("projectionPlane");
//
//		Matrix tV = new Matrix(testVector);
//		tV.printToIJLog("testVector");
//
//		Matrix sV = new Matrix(shaftVector);
//		sV.printToIJLog("shaftVector");
//
//		Matrix nV = new Matrix(bicondylarVector);
//		nV.printToIJLog("neckVector");
//
//		double cosA1 = sV.get(0, 0) * tV.get(0, 0) + sV.get(1, 0)
//				* tV.get(1, 0) + sV.get(2, 0) * tV.get(2, 0);
//		// printMatrix(cosA1, "cosA1");
//		IJ.log("cosA1: " + cosA1);
//
//		double cosA2 = nV.get(0, 0) * tV.get(0, 0) + nV.get(1, 0)
//				* tV.get(1, 0) + nV.get(2, 0) * tV.get(2, 0);
//		// printMatrix(cosA2, "cosA2");
//		IJ.log("cosA2: " + cosA2);
//
//		double bicondylarAngle = Math.acos(cosA1);
//		double neckShaftSkew = Math.acos(cosA2);
//		ResultInserter ri = ResultInserter.getInstance();
//		ri.setResultInRow(imp, "Angle (rad)", bicondylarAngle);
//		ri.setResultInRow(imp, "Skew (rad)", neckShaftSkew);
//		ri.updateTable();
//		return;
//	}
	
	/* Copied directly from Neck_Shaft_Angle */
	/**
	 * Calculate the vector associated with the projection plane from the
	 * regression vector and the vector connecting the centroid and the femoral
	 * head centre
	 * 
	 * @param shaftVector
	 *            double[][]
	 * @param headCentre
	 *            double[][]
	 * @param centroid
	 *            double[][]
	 * @return double[][] projectionPlane
	 */
	public double[][] getProjectionPlane(double[][] shaftVector,
			double[] headCentre, double[] centroid) {
		
		// have to calculate distance between points for a unit vector
		double d = Trig.distance3D(headCentre, centroid);
		double[][] cHVec = new double[3][1];
		cHVec[0][0] = (headCentre[0] - centroid[0]) / d;
		cHVec[1][0] = (headCentre[1] - centroid[1]) / d;
		cHVec[2][0] = (headCentre[2] - centroid[2]) / d;

		Matrix cH = new Matrix(cHVec);
		cH.printToIJLog("cHVec");

		// projectionPlane is the cross product of cHVec and shaftVector
		double[][] projectionPlane = Vectors.crossProduct(cHVec, shaftVector);

		d = Trig.distance3D(projectionPlane[0][0], projectionPlane[1][0], projectionPlane[2][0]);
		projectionPlane[0][0] /= d;
		projectionPlane[1][0] /= d;
		projectionPlane[2][0] /= d;

		return projectionPlane;
	}
	
	
	/**
	 * Copied directly from Neck_Shaft_Angle.regression3D().
	 * NB returns a vector as: double[3][1]
	 * 
	 * Calculate the orthogonal distance regression plane of a set of points by
	 * the covariance method and Singular Value Decomposition
	 * 
	 * @param stack
	 * @param centroid
	 * @return SingularValueDecomposition containing eigenvector and eigenvalue
	 * 
	 * @see <a
	 *      href="http://mathforum.org/library/drmath/view/63765.html">Description
	 *      on Ask Dr Math</a>
	 * 
	 */
	public static double[][] regression3D(ImagePlus imp, double[] centroid,
			int startSlice, int endSlice, double min, double max) {
		IJ.showStatus("Calculating SVD");
		ImageStack stack = imp.getImageStack();
		Rectangle r = stack.getRoi();
		final double cX = centroid[0];
		final double cY = centroid[1];
		final double cZ = centroid[2];
		Calibration cal = imp.getCalibration();
		final double vW = cal.pixelWidth;
		final double vH = cal.pixelHeight;
		final double vD = cal.pixelDepth;
		final int rX = r.x;
		final int rY = r.y;
		final int rW = r.x + r.width;
		final int rH = r.y + r.height;
		double sDxDx = 0;
		double sDyDy = 0;
		double sDzDz = 0;
		double sDxDy = 0;
		double sDxDz = 0;
		double sDyDz = 0;
		double count = 0;
		for (int z = startSlice; z <= endSlice; z++) {
			IJ.showStatus("Getting covariance matrix...");
			IJ.showProgress(z, endSlice);
			final double dz = z * vD - cZ;
			final double dzdz = dz * dz;
			final ImageProcessor ip = stack.getProcessor(z);
			for (int y = rY; y < rH; y++) {
				final double dy = y * vH - cY;
				final double dydy = dy * dy;
				final double dydz = dy * dz;
				for (int x = rX; x < rW; x++) {
					final double testPixel = (double) ip.get(x, y);
					if (testPixel >= min && testPixel <= max) {
						final double dx = x * vW - cX;
						sDxDx += dx * dx;
						sDyDy += dydy;
						sDzDz += dzdz;
						sDxDy += dx * dy;
						sDxDz += dx * dz;
						sDyDz += dydz;
						count++;
					}
				}
			}
		}
		if (count == 0){
			IJ.log("Count == 0");
			return null;
		}
		double[][] C = new double[3][3];
		C[0][0] = sDxDx;
		C[1][1] = sDyDy;
		C[2][2] = sDzDz;
		C[0][1] = sDxDy;
		C[0][2] = sDxDz;
		C[1][0] = sDxDy;
		C[1][2] = sDyDz;
		C[2][0] = sDxDz;
		C[2][1] = sDyDz;
		double invCount = 1 / count;
		Matrix covarianceMatrix = new Matrix(C).times(invCount);
		covarianceMatrix.printToIJLog("Covariance matrix");
		SingularValueDecomposition S = new SingularValueDecomposition(
				covarianceMatrix);
		Matrix leftVectors = S.getU();
		leftVectors.printToIJLog("Left vectors");
		double[][] orthogonalDistanceRegression = new double[3][1];
		orthogonalDistanceRegression[0][0] = leftVectors.get(0, 0);
		orthogonalDistanceRegression[1][0] = leftVectors.get(1, 0);
		orthogonalDistanceRegression[2][0] = leftVectors.get(2, 0);
		return orthogonalDistanceRegression;
	}/* end Regression3D */
	
}
