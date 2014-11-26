package org.bonej;

/**
 * EllipsoidFactor plugin for ImageJ
 * Copyright 2013 Michael Doube
 * 
 *This program is free software: you can redistribute it and/or modify
 *it under the terms of the GNU General Public License as published by
 *the Free Software Foundation, either version 3 of the License, or
 *(at your option) any later version.
 *
 *This program is distributed in the hope that it will be useful,
 *but WITHOUT ANY WARRANTY; without even the implied warranty of
 *MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *GNU General Public License for more details.
 *
 *You should have received a copy of the GNU General Public License
 *along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

import javax.vecmath.Color3f;
import javax.vecmath.Point3f;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.plugin.PlugIn;
import ij.process.ByteProcessor;
import ij.gui.GenericDialog;
import ij.macro.Interpreter;
import ij.measure.Calibration;
import ij3d.Image3DUniverse;

import org.doube.geometry.Trig;
import org.doube.geometry.Vectors;
import org.doube.geometry.Ellipsoid;
import org.doube.jama.Matrix;
import org.doube.skeleton.Skeletonize3D;
import org.doube.util.ArrayHelper;
import org.doube.util.ImageCheck;
import org.doube.util.Multithreader;
//import org.doube.util.ResultInserter;
import org.doube.util.UsageReporter;

import customnode.CustomPointMesh;

/**
 * <p>
 * <b>Plate_Rod</b>
 * </p>
 * <p>
 * ImageJ plugin to describe the local geometry of a binary image in an
 * oblate/prolate spheroid space. Uses Skeletonize3D to generate a 3D skeleton,
 * the points of which are used as centres for star volumes. Local geometry is
 * determined by the ratio between the first and second eigenvalues and first
 * and third eigenvalues of each star volume.
 * </p>
 * 
 * @author Michael Doube
 * 
 */
public class EllipsoidFactor implements PlugIn, Comparator<Ellipsoid> {
	private final byte foreground = (byte) 255;
	private int nVectors = 100;
	private Image3DUniverse universe = new Image3DUniverse();

	/**
	 * increment for vector searching in real units. Defaults to ~Nyquist
	 * sampling of a unit pixel
	 */
	private double vectorIncrement = 1 / 2.3;

	/**
	 * Number of skeleton points per ellipsoid. Sets the granularity of the
	 * ellipsoid fields.
	 */
	private int skipRatio = 50;
	private int contactSensitivity = 5;

	public void run(String arg) {
		if (!ImageCheck.checkEnvironment())
			return;
		ImagePlus imp = IJ.getImage();
		if (imp == null) {
			IJ.noImage();
			return;
		}
		ImageCheck ic = new ImageCheck();
		if (!ic.isBinary(imp) || !ic.isMultiSlice(imp)
				|| !ic.isVoxelIsotropic(imp, 0.001)) {
			IJ.error("8-bit binary stack with isotropic pixel spacing required.");
			return;
		}
		Calibration cal = imp.getCalibration();
		String units = cal.getUnits();
		vectorIncrement *= Math.min(cal.pixelDepth,
				Math.min(cal.pixelHeight, cal.pixelWidth));
		GenericDialog gd = new GenericDialog("Setup");
		gd.addNumericField("Sampling increment", vectorIncrement, 3, 8, units);
		gd.addNumericField("Vectors", nVectors, 0, 8, "");
		gd.addNumericField("Skeleton points per ellipsoid", skipRatio, 0);
		gd.addNumericField("Contact sensitivity", contactSensitivity, 0, 4, "");
		gd.addHelp("http://bonej.org/ef");
		gd.showDialog();
		if (!Interpreter.isBatchMode()) {
			vectorIncrement = gd.getNextNumber();
			nVectors = (int) Math.round(gd.getNextNumber());
			skipRatio = (int) Math.round(gd.getNextNumber());
			contactSensitivity = (int) Math.round(gd.getNextNumber());
		}
		if (gd.wasCanceled())
			return;

		final double[][] unitVectors = Vectors.regularVectors(nVectors);
		int[][] skeletonPoints = skeletonPoints(imp);

		IJ.log("Found " + skeletonPoints.length + " skeleton points");

		if (IJ.debugMode)
			universe.show();

		Ellipsoid[] ellipsoids = findEllipsoids(imp, skeletonPoints,
				unitVectors);

		IJ.log("Found " + ellipsoids.length + " ellipsoids");

		float[][] biggestEllipsoid = findBiggestEllipsoid(imp, ellipsoids);

		ImageStack bigStack = new ImageStack(imp.getWidth(), imp.getHeight());
		for (int i = 1; i < biggestEllipsoid.length; i++)
			bigStack.addSlice("" + i, biggestEllipsoid[i]);

		ImagePlus bigImp = new ImagePlus("", bigStack);
		bigImp.setDisplayRange(-ellipsoids.length / 2, ellipsoids.length);
		bigImp.show();

		// ResultInserter ri = ResultInserter.getInstance();
		// ri.updateTable();
		UsageReporter.reportEvent(this).send();
	}

	/**
	 * For each foreground pixel of the input image, find the ellipsoid of
	 * greatest volume
	 * 
	 * @param imp
	 * @param ellipsoids
	 * @return array containing the indexes of the biggest ellipsoids which
	 *         contain each point
	 */
	private float[][] findBiggestEllipsoid(ImagePlus imp,
			final Ellipsoid[] ellipsoids) {

		final ImageStack stack = imp.getImageStack();
		final int w = stack.getWidth();
		final int h = stack.getHeight();
		final int d = stack.getSize();

		Calibration cal = imp.getCalibration();
		final double vW = cal.pixelWidth;
		final double vH = cal.pixelHeight;
		final double vD = cal.pixelDepth;

		final float[][] biggest = new float[d + 1][w * h];

		final AtomicInteger ai = new AtomicInteger(1);
		Thread[] threads = Multithreader.newThreads();
		for (int thread = 0; thread < threads.length; thread++) {
			threads[thread] = new Thread(new Runnable() {
				public void run() {
					for (int z = ai.getAndIncrement(); z <= d; z = ai
							.getAndIncrement()) {
						byte[] slicePixels = (byte[]) stack.getPixels(z);
						float[] bigSlice = biggest[z];
						Arrays.fill(bigSlice, -ellipsoids.length);
						final double zvD = z * vD;
						for (int y = 0; y < h; y++) {
							final int offset = y * w;
							final double yvH = y * vH;
							for (int x = 0; x < w; x++) {
								if (slicePixels[offset + x] == foreground) {
									bigSlice[offset + x] = biggestEllipsoid(
											ellipsoids, x * vW, yvH, zvD);
								}
							}
						}

					}
				}
			});
		}
		Multithreader.startAndJoin(threads);

		return biggest;
	}

	/**
	 * Search the list of ellipsoids and return the index of the largest
	 * ellipsoid which contains the point x, y, z
	 * 
	 * @param ellipsoids
	 * @param x
	 * @param y
	 * @param z
	 * @return the index of the largest ellipsoid which contains this point
	 */
	private int biggestEllipsoid(Ellipsoid[] ellipsoids, double x, double y,
			double z) {
		final int l = ellipsoids.length;
		for (int i = 0; i < l; i++) {
			if (ellipsoids[i].contains(x, y, z))
				return i;
		}
		return -1;
	}

	/**
	 * Using skeleton points as seeds, propagate along each vector until a
	 * boundary is hit. Use the resulting cloud of boundary points as input into
	 * an ellipsoid fit.
	 * 
	 * @param imp
	 * @param skeletonPoints
	 * @param unitVectors
	 * @return
	 */
	private Ellipsoid[] findEllipsoids(final ImagePlus imp,
			final int[][] skeletonPoints, final double[][] unitVectors) {
		final int nPoints = skeletonPoints.length;
		final Ellipsoid[] ellipsoids = new Ellipsoid[nPoints];

		// make sure array contains null in the non-calculated elements
		Arrays.fill(ellipsoids, null);

		final AtomicInteger ai = new AtomicInteger(1);
		Thread[] threads = Multithreader.newThreads();
		for (int thread = 0; thread < threads.length; thread++) {
			threads[thread] = new Thread(new Runnable() {
				public void run() {
					for (int i = ai.getAndAdd(skipRatio); i <= nPoints; i = ai
							.getAndAdd(skipRatio)) {
						IJ.showStatus("Optimising ellipsoid " + (i + 1) + "/"
								+ nPoints);
						ellipsoids[i] = optimiseEllipsoid(imp,
								skeletonPoints[i], unitVectors);
					}
				}
			});
		}
		Multithreader.startAndJoin(threads);

		Ellipsoid[] sortedEllipsoids = ArrayHelper.removeNulls(ellipsoids);

		// Sort using this class' compare method
		Arrays.sort(sortedEllipsoids, this);
		return sortedEllipsoids;
	}

	/**
	 * given a seed point, find the ellipsoid which best fits the binarised
	 * structure
	 * 
	 * @param imp
	 * @param is
	 * @param unitVectors
	 * @return ellipsoid fitting the point cloud of boundaries lying at the end
	 *         of vectors surrounding the seed point. If ellipsoid fitting
	 *         fails, returns null
	 */
	private Ellipsoid optimiseEllipsoid(final ImagePlus imp,
			int[] skeletonPoint, double[][] unitVectors) {

		Calibration cal = imp.getCalibration();
		final double pW = cal.pixelWidth;
		final double pH = cal.pixelHeight;
		final double pD = cal.pixelDepth;

		ImageStack stack = imp.getImageStack();
		final int stackSize = stack.getSize();

		// cache slices into an array
		ByteProcessor[] ips = new ByteProcessor[stackSize + 1];
		for (int i = 1; i <= stackSize; i++) {
			ips[i] = (ByteProcessor) stack.getProcessor(i);
		}

		final int w = ips[1].getWidth();
		final int h = ips[1].getHeight();
		final int d = ips.length - 1;

		// centre point of vector field
		final double px = skeletonPoint[0] * pW;
		final double py = skeletonPoint[1] * pH;
		final double pz = skeletonPoint[2] * pD;

		// Instantiate a small spherical ellipsoid
		final double[][] orthogonalVectors = { { 1, 0, 0 }, { 0, 1, 0 },
				{ 0, 0, 1 } };

		Ellipsoid ellipsoid = new Ellipsoid(0, 0, 0, px, py, pz,
				orthogonalVectors);

		// dilate the sphere until it hits the background
		while (isContained(ellipsoid, ips, pW, pH, pD, w, h, d)) {
			ellipsoid.dilate(vectorIncrement, vectorIncrement, vectorIncrement);
		}

		// get the points of contact
		ArrayList<double[]> contactPoints = findContactPoints(ellipsoid, ips,
				pW, pH, pD, w, h, d);

		// contract the ellipsoid by one increment so all points are
		// inside the foregrounds
		// ellipsoid.contract(vectorIncrement);

		// find the mean unit vector pointing to the points of contact from the
		// centre
		double[] summedVector = new double[3];
		final double[] c = ellipsoid.getCentre();
		for (double[] p : contactPoints) {
			final double l = Trig.distance3D(p, c);
			double[] unitVector = { (p[0] - c[0]) / l, (p[1] - c[1]) / l,
					(p[2] - c[2]) / l };
			summedVector[0] += unitVector[0];
			summedVector[1] += unitVector[1];
			summedVector[2] += unitVector[2];
		}
		// put the short axis on the mean contact vector
		double[] shortAxis = new double[3];
		shortAxis[0] = summedVector[0] / contactPoints.size();
		shortAxis[1] = summedVector[1] / contactPoints.size();
		shortAxis[2] = summedVector[2] / contactPoints.size();

		shortAxis = Vectors.norm(shortAxis);

		// find an orthogonal axis
		final double[] xAxis = { 1, 0, 0 };
		double[] middleAxis = Vectors.crossProduct(shortAxis, xAxis);
		middleAxis = Vectors.norm(middleAxis);

		// find a mutually orthogonal axis by forming the cross product
		double[] longAxis = Vectors.crossProduct(shortAxis, middleAxis);
		longAxis = Vectors.norm(longAxis);

		// construct a rotation matrix
		double[][] rotation = { shortAxis, middleAxis, longAxis };

		// needs transpose because each vector is put in as row to begin with
		Matrix R = new Matrix(rotation).transpose();

		// R.printToIJLog("Rotation Matrix: det() = " + R.det());

		// rotate ellipsoid to point this way...
		ellipsoid.setRotation(R);

		// dilate other two axes until number of contact points increases
		// by contactSensitivity number of contacts

		int maxContacts = contactPoints.size() + contactSensitivity;
		while (contactPoints.size() < maxContacts) {
			ellipsoid.dilate(0, vectorIncrement, vectorIncrement);
			contactPoints = findContactPoints(ellipsoid, ips, pW, pH, pD, w, h,
					d);
		}

		// until ellipsoid is totally jammed within the structure, go through
		// cycles of contraction, wiggling, dilation
		// goal is maximal inscribed ellipsoid, maximal being defined by volume

//		double maxVolumeSoFar = ellipsoid.getVolume();
//
//		// maybe store a copy of the 'best ellipsoid so far'?
//		// needs TODO copy method in Ellipsoid.
//		double volumeThisTime = 0;
//		int noImprovementCount = 0;
//		// number of times to cycle without improvement before quitting
//		final int triedEnoughTimes = 5;
//		int totalIterations = 0;
//		while (noImprovementCount < triedEnoughTimes && totalIterations < 100) {
//
//			// contract until no contact
//			while (contactPoints.size() > 0) {
//				ellipsoid.contract(0.01);
//				contactPoints = findContactPoints(ellipsoid, ips, pW, pH, pD,
//						w, h, d);
//			}
//
//			// rotate a little bit
//			wiggle(ellipsoid);
//
//			// dilate a
//			while (contactPoints.size() < contactSensitivity) {
//				ellipsoid.dilate(vectorIncrement, 0, 0);
//				contactPoints = findContactPoints(ellipsoid, ips, pW, pH, pD,
//						w, h, d);
//			}
//
////			shrink(ellipsoid, contactPoints, ips, pW, pH, pD, w, h, d);
//
//			// dilate b
//			while (contactPoints.size() < contactSensitivity) {
//				ellipsoid.dilate(0, vectorIncrement, 0);
//				contactPoints = findContactPoints(ellipsoid, ips, pW, pH, pD,
//						w, h, d);
//			}
//
////			shrink(ellipsoid, contactPoints, ips, pW, pH, pD, w, h, d);
//
//			// dilate c
//			while (contactPoints.size() < contactSensitivity) {
//				ellipsoid.dilate(0, vectorIncrement, 0);
//				contactPoints = findContactPoints(ellipsoid, ips, pW, pH, pD,
//						w, h, d);
//			}
//
//			volumeThisTime = ellipsoid.getVolume();
//			if (volumeThisTime <= maxVolumeSoFar)
//				noImprovementCount++;
//			else
//				maxVolumeSoFar = volumeThisTime;
//			totalIterations++;
//		}
		// add them to the 3D viewer
		if (IJ.debugMode) {
			ArrayList<Point3f> contactPointsf = new ArrayList<Point3f>(
					contactPoints.size());
			for (double[] p : contactPoints) {
				Point3f point = new Point3f((float) p[0], (float) p[1],
						(float) p[2]);
				contactPointsf.add(point);
			}
			double[][] pointCloud = ellipsoid.getSurfacePoints(100);

			List<Point3f> pointList = new ArrayList<Point3f>();
			for (int p = 0; p < pointCloud.length; p++) {
				if (pointCloud[p] == null)
					continue;
				Point3f e = new Point3f();
				e.x = (float) pointCloud[p][0];
				e.y = (float) pointCloud[p][1];
				e.z = (float) pointCloud[p][2];
				pointList.add(e);
			}
			CustomPointMesh mesh = new CustomPointMesh(pointList);
			mesh.setPointSize(2.0f);
			Color3f cColour = new Color3f((float) (px / pW) / w,
					(float) (py / pH) / h, (float) (pz / pD) / d);
			mesh.setColor(cColour);

			CustomPointMesh contactPointMesh = new CustomPointMesh(
					contactPointsf);
			contactPointMesh.setPointSize(2.5f);
			Color3f invColour = new Color3f(1 - cColour.x, 1 - cColour.y,
					1 - cColour.z);
			contactPointMesh.setColor(invColour);

			try {
				universe.addCustomMesh(mesh,
						"Point cloud " + px + " " + py + " " + pz).setLocked(
						true);
				universe.addCustomMesh(contactPointMesh,
						"Contact points of " + px + " " + py + " " + pz)
						.setLocked(true);

			} catch (NullPointerException npe) {
				IJ.log("3D Viewer was closed before rendering completed.");
			}
		}
		return ellipsoid;
	}

	private void shrink(Ellipsoid ellipsoid, ArrayList<double[]> contactPoints,
			ByteProcessor[] ips, double pW, double pH, double pD, int w, int h,
			int d) {

		// contract until no contact
		while (contactPoints.size() > 0) {
			ellipsoid.contract(0.001);
			contactPoints = findContactPoints(ellipsoid, ips, pW, pH, pD, w, h,
					d);
		}
	}

	/**
	 * Rotate the ellipsoid by a small random amount
	 * 
	 * @param ellipsoid
	 */
	private void wiggle(Ellipsoid ellipsoid) {

		double b = nudge(0.1);
		double c = nudge(0.1);
		double a = Math.sqrt(1 - b * b - c * c);

		// zeroth column, should be very close to [1, 0, 0]^T (mostly x)
		double[] zerothColumn = { a, b, c };

		// form triangle in nearly xy plane
		double[] vector = { 0, 1, 0 };

		// first column, should be very close to [0, 1, 0]^T
		double[] firstColumn = Vectors.norm(Vectors.crossProduct(zerothColumn,
				vector));

		// second column, should be very close to [0, 0, 1]^T
		double[] secondColumn = Vectors.norm(Vectors.crossProduct(zerothColumn,
				firstColumn));

		double[][] rotation = { zerothColumn, firstColumn, secondColumn };

		// rotation has vectors as rows, to need to transpose
		Matrix N = new Matrix(rotation).transpose();

		// N.printToIJLog("Wiggle rotation matrix");
		ellipsoid.rotate(N);
	}

	/**
	 * generate a random number between -a and +a
	 * 
	 * @param a
	 * @return
	 */
	private double nudge(double a) {
		return Math.random() * (a + a) - a;
	}

	private ArrayList<double[]> findContactPoints(Ellipsoid ellipsoid,
			ByteProcessor[] ips, final double pW, final double pH,
			final double pD, final int w, final int h, final int d) {
		double[][] points = ellipsoid.getSurfacePoints(nVectors);

		ArrayList<double[]> contactPoints = new ArrayList<double[]>();

		for (double[] p : points) {
			final int x = (int) Math.floor(p[0] / pW);
			final int y = (int) Math.floor(p[1] / pH);
			final int z = (int) Math.floor(p[2] / pD);
			if (isOutOfBounds(x, y, z, w, h, d))
				continue;
			if ((byte) ips[z].get(x, y) != foreground)
				contactPoints.add(p);
		}
		return contactPoints;
	}

	private boolean isContained(Ellipsoid ellipsoid, ByteProcessor[] ips,
			final double pW, final double pH, final double pD, final int w,
			final int h, final int d) {
		double[][] points = ellipsoid.getSurfacePoints(nVectors);
		for (double[] p : points) {
			final int x = (int) Math.floor(p[0] / pW);
			final int y = (int) Math.floor(p[1] / pH);
			final int z = (int) Math.floor(p[2] / pD);
			if (isOutOfBounds(x, y, z, w, h, d))
				continue;
			if ((byte) ips[z].get(x, y) != foreground)
				return false;
		}
		return true;
	}

	/**
	 * return true if pixel coordinate is out of image bounds
	 * 
	 * @param x
	 * @param y
	 * @param z
	 * @param w
	 * @param h
	 * @param d
	 * @return
	 */
	private boolean isOutOfBounds(int x, int y, int z, int w, int h, int d) {
		if (x < 0 || x >= w || y < 0 || y >= h || z < 1 || z > d)
			return true;
		else
			return false;
	}

	private int[][] skeletonPoints(ImagePlus imp) {
		Skeletonize3D sk = new Skeletonize3D();
		ImagePlus skeleton = sk.getSkeleton(imp);
		final ImageStack skeletonStack = skeleton.getStack();

//		if (IJ.debugMode)
//			skeleton.show();

		final int d = imp.getStackSize();
		final int h = imp.getHeight();
		final int w = imp.getWidth();

		IJ.log("Skeleton image is " + w + " x " + h + " x " + d);

		// Bare ArrayList is not thread safe for concurrent add() operations.
		final List<int[]> list = Collections
				.synchronizedList(new ArrayList<int[]>());

		final AtomicInteger ai = new AtomicInteger(1);
		Thread[] threads = Multithreader.newThreads();
		for (int thread = 0; thread < threads.length; thread++) {
			threads[thread] = new Thread(new Runnable() {
				public void run() {
					for (int z = ai.getAndIncrement(); z <= d; z = ai
							.getAndIncrement()) {
						byte[] slicePixels = (byte[]) skeletonStack
								.getPixels(z);
						for (int y = 0; y < h; y++) {
							final int offset = y * w;
							for (int x = 0; x < w; x++) {
								if (slicePixels[offset + x] == foreground) {
									final int[] array = { x, y, z };
									list.add(array);
								}
							}
						}
					}
				}
			});
		}
		Multithreader.startAndJoin(threads);

		if (IJ.debugMode)
			IJ.log("Skeleton point ArrayList contains " + list.size()
					+ " points");

		int[][] skeletonPoints = list.toArray(new int[list.size()][]);

		return skeletonPoints;
	}

	/**
	 * Compare Ellipsoids by volume.
	 * 
	 * Sorting based on this method will result in Ellipsoids sorted in order of
	 * <b>descending</b> volume.
	 * 
	 */
	public int compare(Ellipsoid o1, Ellipsoid o2) {
		return Double.compare(o2.getVolume(), o1.getVolume());
	}

}
