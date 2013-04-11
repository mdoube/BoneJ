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

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.plugin.PlugIn;
import ij.process.ByteProcessor;
import ij.gui.GenericDialog;
import ij.macro.Interpreter;
import ij.measure.Calibration;

import org.doube.geometry.FitEllipsoid;
import org.doube.geometry.Vectors;
import org.doube.geometry.Ellipsoid;
import org.doube.skeleton.Skeletonize3D;
import org.doube.util.ImageCheck;
import org.doube.util.ResultInserter;
import org.doube.util.UsageReporter;

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
/*
 * TODO Summarise plateness / rodness, e.g. ratio of sums of middle and biggest
 * Eigenvalues plot ev2/ev1 on 3d skeleton
 */

public class EllipsoidFactor implements PlugIn, Comparator<Ellipsoid> {
	private int nVectors = 1000;
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
		if (!ic.isBinary(imp) || !ic.isMultiSlice(imp)) {
			IJ.error("8-bit binary stack required.");
			return;
		}
		Calibration cal = imp.getCalibration();
		final double vD = cal.pixelDepth;
		final double vH = cal.pixelHeight;
		final double vW = cal.pixelWidth;
		String units = cal.getUnits();
		vectorIncrement = Math.max(vH, Math.max(vW, vD));
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
		Ellipsoid[] ellipsoids = findEllipsoids(imp, skeletonPoints,
				unitVectors);

		ResultInserter ri = ResultInserter.getInstance();
		ri.updateTable();
		UsageReporter.reportEvent(this).send();
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

		for (int i = 0; i < nPoints; i++) {
			ellipsoids[i] = optimiseEllipsoid(imp, skeletonPoints[i],
					unitVectors);
		}

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

		final int w = ips[0].getWidth();
		final int h = ips[1].getHeight();
		final int d = ips.length - 1;

		// centre point of vector field
		final int px = skeletonPoint[0];
		final int py = skeletonPoint[1];
		final int pz = skeletonPoint[2];

		double[][] pointCloud = new double[nVectors][3];

		// for each vector
		for (int v = 0; v < nVectors; v++) {
			final double[] vector = unitVectors[v];
			// search the direction
			final double dx = vector[0];
			final double dy = vector[1];
			final double dz = vector[2];

			// TODO check that foreground actually is (byte) 255
			final byte foreground = (byte) 255;
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
					break;
				} else
					pixel = (byte) ips[z].get(x, y);
			}
			pointCloud[v][0] = (int) Math.floor(px + l * dx) * pW;
			pointCloud[v][1] = (int) Math.floor(py + l * dy) * pH;
			pointCloud[v][2] = (int) Math.floor(pz + l * dz) * pD;
		}

		try {
			Ellipsoid ellipsoid = FitEllipsoid.fitTo(pointCloud);
			return ellipsoid;
		} catch (Exception e) {
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
		ImageStack skeletonStack = sk.getSkeleton(imp).getStack();

		final int d = imp.getStackSize();
		final int h = imp.getHeight();
		final int w = imp.getWidth();

		ArrayList<int[]> list = new ArrayList<int[]>();

		for (int z = 1; z <= d; z++) {
			byte[] slicePixels = (byte[]) skeletonStack.getPixels(z);
			for (int y = 0; y < h; y++) {
				int offset = y * w;
				for (int x = 0; x < w; x++) {
					if (slicePixels[offset + x] < 0) {
						int[] array = { x, y, z };
						list.add(array);
					}
				}
			}
		}

		int[][] skeletonPoints = list.toArray(new int[list.size()][]);

		return skeletonPoints;
	}

	/**
	 * Compare Ellipsoids by volume. Ordering could be reversed by swapping o1
	 * and o2
	 * 
	 */
	public int compare(Ellipsoid o1, Ellipsoid o2) {
		return Double.compare(o1.getVolume(), o2.getVolume());
	}
}
