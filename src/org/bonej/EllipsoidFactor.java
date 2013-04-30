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
import java.util.Comparator;
import java.util.List;

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

import org.doube.geometry.FitEllipsoid;
import org.doube.geometry.Vectors;
import org.doube.geometry.Ellipsoid;
import org.doube.skeleton.Skeletonize3D;
import org.doube.util.ArrayHelper;
import org.doube.util.ImageCheck;
import org.doube.util.ResultInserter;
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
	private int nVectors = 1000;
	private Image3DUniverse universe = new Image3DUniverse();
	
	/**
	 * increment for vector searching in pixel units. Defaults to ~Nyquist
	 * sampling
	 */
	private double vectorIncrement = 1 / 2.3;

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
		GenericDialog gd = new GenericDialog("Setup");
		gd.addNumericField("Sampling increment", vectorIncrement, 3, 8, units);
		gd.addNumericField("Vectors", nVectors, 0, 8, "");
		gd.addHelp("http://bonej.org/ef");
		gd.showDialog();
		if (!Interpreter.isBatchMode()) {
			vectorIncrement = gd.getNextNumber();
			nVectors = (int) Math.round(gd.getNextNumber());
		}
		if (gd.wasCanceled())
			return;
		
		final double[][] unitVectors = Vectors.regularVectors(nVectors);
		int[][] skeletonPoints = skeletonPoints(imp);

		IJ.log("Found " + skeletonPoints.length + " skeleton points");

		Ellipsoid[] ellipsoids = findEllipsoids(imp, skeletonPoints,
				unitVectors);

		IJ.log("Found " + ellipsoids.length + " ellipsoids");

//		Image3DUniverse univ = new Image3DUniverse();
		
		for (Ellipsoid e : ellipsoids) {
			IJ.log("" + e.getVolume());
			displayEllipsoid(universe, e, 1000);
		}
		
		universe.show();

		float[][] biggestEllipsoid = findBiggestEllipsoid(imp, ellipsoids);

		ImageStack bigStack = new ImageStack(imp.getWidth(), imp.getHeight());
		for (int i = 1; i < biggestEllipsoid.length; i++)
			bigStack.addSlice("" + i, biggestEllipsoid[i]);
		
		ImagePlus bigImp = new ImagePlus("", bigStack);
		bigImp.setDisplayRange(-ellipsoids.length/2, ellipsoids.length);
		bigImp.show();

		ResultInserter ri = ResultInserter.getInstance();
		ri.updateTable();
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
	private float[][] findBiggestEllipsoid(ImagePlus imp, Ellipsoid[] ellipsoids) {

		ImageStack stack = imp.getImageStack();
		final int w = stack.getWidth();
		final int h = stack.getHeight();
		final int d = stack.getSize();

		Calibration cal = imp.getCalibration();
		final double vW = cal.pixelWidth;
		final double vH = cal.pixelHeight;
		final double vD = cal.pixelDepth;

		float[][] biggest = new float[d + 1][w * h];

		for (int z = 1; z <= d; z++) {
			byte[] slicePixels = (byte[]) stack.getPixels(z);
			float[] bigSlice = biggest[z];
			// -1 means background, 0 will be the biggest ellipsoid
			Arrays.fill(bigSlice, -ellipsoids.length);
			for (int y = 0; y < h; y++) {
				int offset = y * w;
				for (int x = 0; x < w; x++) {
					if (slicePixels[offset + x] == foreground) {
						bigSlice[offset + x] = biggestEllipsoid(ellipsoids, x
								* vW, y * vH, z * vD);
					}
				}
			}
		}

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
	private Ellipsoid[] findEllipsoids(ImagePlus imp, int[][] skeletonPoints,
			double[][] unitVectors) {
		final int nPoints = skeletonPoints.length;
		Ellipsoid[] ellipsoids = new Ellipsoid[nPoints];

		//TODO tweak number of skipped skeleton points
		for (int i = 0; i < nPoints; i += 20) {
			ellipsoids[i] = optimiseEllipsoid(imp, skeletonPoints[i],
					unitVectors);
		}
		
		ellipsoids = ArrayHelper.removeNulls(ellipsoids);
		
		// Sort using this class' compare method
		Arrays.sort(ellipsoids, this);
		return ellipsoids;
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

		// cache slices into an array
		ByteProcessor[] ips = new ByteProcessor[stack.getSize() + 1];
		for (int i = 1; i <= stack.getSize(); i++) {
			ips[i] = (ByteProcessor) stack.getProcessor(i);
		}

		final int w = ips[1].getWidth();
		final int h = ips[1].getHeight();
		final int d = ips.length - 1;

		// centre point of vector field
		final int px = skeletonPoint[0];
		final int py = skeletonPoint[1];
		final int pz = skeletonPoint[2];

		double[][] pointCloud = new double[nVectors][3];

		// for each vector
		vectorLoop:
		for (int v = 0; v < nVectors; v++) {
			final double[] vector = unitVectors[v];
			// search the direction
			final double dx = vector[0];
			final double dy = vector[1];
			final double dz = vector[2];

			byte pixel = foreground;
			double l = 0;
			while (pixel == foreground) {
				l += vectorIncrement;
				// pixel is defined as minimal corner of box
				final int x = (int) Math.floor(px + l * dx);
				final int y = (int) Math.floor(py + l * dy);
				final int z = (int) Math.floor(pz + l * dz);

				// set points to null when their vector goes out of bounds
				if (isOutOfBounds(x, y, z, w, h, d)) {
					pointCloud[v] = null;
					continue vectorLoop;
				} else
					pixel = (byte) ips[z].get(x, y);
			}

			// point is in real units, simply pixel location times pixel spacing
			double[] point = { ((int) Math.floor(px + l * dx)) * pW,
					((int) Math.floor(py + l * dy)) * pH,
					((int) Math.floor(pz + l * dz)) * pD };

			pointCloud[v] = point;
		}

		List<Point3f> pointList = new ArrayList<Point3f>();
		for (int p = 0; p < nVectors; p++) {
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
		Color3f cColour = new Color3f((float)px/w, (float)py/h, (float)pz/d);
		mesh.setColor(cColour);
		try {
			universe.addCustomMesh(mesh, "Point cloud "+px+" "+py+" "+pz).setLocked(false);
		} catch (NullPointerException npe) {
			IJ.log("3D Viewer was closed before rendering completed.");
		}
		
		
		try {
			Ellipsoid ellipsoid = FitEllipsoid.inertia(ArrayHelper.removeNulls(pointCloud));
			return ellipsoid;
		} catch (Exception e) {
			IJ.log("Couldn't fit ellipsoid: "+e.getMessage());
			return null;
		}
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
		ImageStack skeletonStack = skeleton.getStack();
		
		skeleton.show();

		final int d = imp.getStackSize();
		final int h = imp.getHeight();
		final int w = imp.getWidth();

		IJ.log("Skeleton image is " + w + " x " + h + " x " + d);

		ArrayList<int[]> list = new ArrayList<int[]>();

		for (int z = 1; z <= d; z++) {
			byte[] slicePixels = (byte[]) skeletonStack.getPixels(z);
			for (int y = 0; y < h; y++) {
				int offset = y * w;
				for (int x = 0; x < w; x++) {
					if (slicePixels[offset + x] == foreground) {
						int[] array = { x, y, z };
						list.add(array);
					}
				}
			}
		}

		IJ.log("Skeleton point ArrayList contains " + list.size() + " points");

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
	
	private void displayEllipsoid(Image3DUniverse univ, Ellipsoid ellipsoid, int nPoints){
		
		final double[][] points = ellipsoid.getSurfacePoints(nPoints);
		
		List<Point3f> pointList = new ArrayList<Point3f>();
		for (int p = 0; p < nPoints; p++) {
			Point3f e = new Point3f();
			e.x = (float) points[p][0];
			e.y = (float) points[p][1];
			e.z = (float) points[p][2];
			pointList.add(e);
		}
		CustomPointMesh mesh = new CustomPointMesh(pointList);
		mesh.setPointSize(2.0f);
		Color3f cColour = new Color3f(0.0f, 0.5f, 1.0f);
		mesh.setColor(cColour);
		try {
			univ.addCustomMesh(mesh, "Ellipsoid "+ellipsoid.getVolume()).setLocked(true);
		} catch (NullPointerException npe) {
			IJ.log("3D Viewer was closed before rendering completed.");
			return;
		}
		return;
		
	}
	
}
