/**
 *Anisotropy_ plugin for ImageJ
 *Copyright 2009 Michael Doube 
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

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ImageProcessor;
import ij.plugin.PlugIn;
import ij.gui.*;
import ij.measure.Calibration; //import ij.measure.ResultsTable;

// for 3D plotting of coordinates
import javax.vecmath.Point3f;
import javax.vecmath.Color3f;
import customnode.CustomPointMesh;

import java.util.Enumeration;
import java.util.List;
import java.util.ArrayList;
import java.util.Vector;

import ij3d.Image3DUniverse;
import ij3d.Content;

import org.doube.bonej.ImageCheck;
import org.doube.bonej.ResultInserter;
import org.doube.jama.*;

/**
 * <p>
 * Calculate degree of anisotropy using the mean intercept length method
 * </p>
 * 
 * @see <p>
 *      Odgaard A (1997) Three-dimensional methods for quantification of
 *      cancellous bone architecture. Bone 20: 315-28. <a
 *      href="http://dx.doi.org/10.1016/S8756-3282(97)00007-0"
 *      >doi:10.1016/S8756-3282(97)00007-0</a>
 *      </p>
 * @author Michael Doube
 * 
 */
// TODO implement the star volume method
// TODO implement the autocorrelation function (ACF) method
// TODO multithread
// TODO split off anisotropy algorithms into classes in org.doube.bonej
// and call them from the main Anisotropy_ class
// TODO run to stable result.
public class Anisotropy_ implements PlugIn {

	/** Graph of results */
	private ImagePlus plotImage;

	ImageProcessor ip;

	private double[][] coOrdinates;

	public void run(String arg) {
		if (!ImageCheck.checkIJVersion()) {
			return;
		}
		ImagePlus imp = IJ.getImage();
		if (null == imp) {
			IJ.noImage();
			return;
		}
		ImageCheck ic = new ImageCheck();
		if (!ic.isBinary(imp)) {
			IJ.error("8-bit binary (black and white only) image required.");
			return;
		}
		if (!ic.isMultiSlice(imp)) {
			IJ.error("Stack required");
			return;
		}
		/*
		 * if (!showDialog()) { return; }
		 */
		Calibration cal = imp.getCalibration();
		final double vW = cal.pixelWidth;
		final double vH = cal.pixelHeight;
		final double vD = cal.pixelDepth;
		final int w = imp.getWidth();
		final int h = imp.getHeight();
		final int d = imp.getStackSize();
		final double vectorSampling = Math.max(vW, Math.max(vH, vD)) * 2.3;
		final double radius = Math.min(h * vH, Math.min(d * vD, w * vW)) / 4;

		GenericDialog gd = new GenericDialog("Setup");
		gd.addCheckbox("Auto Mode", true);
		// number of random vectors in vector field
		gd.addNumericField("Vectors", 50000, 0, 6, "vectors");
		// number of randomly-positioned vector fields
		gd.addNumericField("Min_Spheres", 100, 0, 5, "");
		gd.addNumericField("Max_Spheres", 2000, 0, 5, "");
		gd.addCheckbox("Show_Plot", true);
		gd.addCheckbox("3D_Result", false);
		gd.showDialog();
		if (gd.wasCanceled()) {
			return;
		}
		final boolean doAutoMode = gd.getNextBoolean();
		final int nVectors = (int) gd.getNextNumber();
		final int minSpheres = (int) gd.getNextNumber();
		final int maxSpheres = (int) gd.getNextNumber();
		final boolean doPlot = gd.getNextBoolean();
		final boolean do3DResult = gd.getNextBoolean();

		double anisotropy = Double.NaN;
		if (doAutoMode)
			anisotropy = runToStableResult(imp, minSpheres, maxSpheres,
					nVectors, radius, vectorSampling, doPlot);
		else
			anisotropy = runToStableResult(imp, minSpheres, minSpheres,
					nVectors, radius, vectorSampling, doPlot);

		ResultInserter ri = ResultInserter.getInstance();
		ri.setResultInRow(imp, "Anisotropy", anisotropy);
		ri.updateTable();

		if (do3DResult) {
			plotPoints3D(coOrdinates, "Intercept Lengths");
		}
	}

	/**
	 * 
	 * @param imp
	 * @param minSpheres
	 * @param nVectors
	 * @param radius
	 * @param vectorSampling
	 * @return
	 * @deprecated
	 */
	@SuppressWarnings("unused")
	private double runOnce(ImagePlus imp, int minSpheres, int nVectors,
			double radius, double vectorSampling) {
		final int nSpheres = minSpheres;
		double[][] centroidList = gridCalculator(imp, nSpheres, radius);
		double[][] vectorList = randomVectors(nVectors);
		double[] sumInterceptCounts = new double[nVectors];
		for (int s = 0; s < nSpheres; s++) {
			double[] centroid = new double[3];
			centroid[0] = centroidList[s][0];
			centroid[1] = centroidList[s][1];
			centroid[2] = centroidList[s][2];
			IJ.showStatus("Counting intercepts at site " + (s + 1) + "/"
					+ nSpheres);
			double[] interceptCounts = countIntercepts(imp, centroid,
					vectorList, nVectors, radius, vectorSampling);
			for (int i = 0; i < nVectors; i++) {
				sumInterceptCounts[i] += interceptCounts[i];
			}
		}
		double[] meanInterceptLengths = new double[nVectors];
		for (int v = 0; v < nVectors; v++) {
			// in case by chance there has been no intercept, we want to
			// avoid a divide by 0 error
			if (sumInterceptCounts[v] == 0) {
				sumInterceptCounts[v] = 1;
			}
			// MIL = total vector length / number of intercepts
			meanInterceptLengths[v] = radius * (double) nSpheres
					/ sumInterceptCounts[v];
		}
		coOrdinates = new double[nVectors][3];
		for (int v = 0; v < nVectors; v++) {
			coOrdinates[v][0] = meanInterceptLengths[v] * vectorList[v][0];
			coOrdinates[v][1] = meanInterceptLengths[v] * vectorList[v][1];
			coOrdinates[v][2] = meanInterceptLengths[v] * vectorList[v][2];
		}
		// the ellipsoid is cool but a bit tricky to work with, sticking with
		// PCA for now.
		// double[][] ellipsoid = fitEllipsoid(coOrdinates);
		EigenvalueDecomposition E = principalComponents(coOrdinates);
		Matrix eigenValues = E.getD();
		Matrix eigenVectors = E.getV();
		IJ.log("eigenValues:");
		printMatrix(eigenValues);
		IJ.log("eigenVectors");
		printMatrix(eigenVectors);
		double[][] eVal = eigenValues.getArrayCopy();
		// double[][] eVec = eigenVectors.getArrayCopy();
		double anisotropy = 1 - eVal[0][0] / eVal[2][2];
		return anisotropy;
	}

	private double runToStableResult(ImagePlus imp, int minSpheres,
			int maxSpheres, int nVectors, double radius, double vectorSampling,
			boolean doPlot) {
		final int minIterations = minSpheres;
		final int maxIterations = maxSpheres;
		final double[][] vectorList = randomVectors(nVectors);
		double variance = Double.MAX_VALUE;
		final double tolerance = 0.0001;
		double anisotropy = Double.NaN;
		double[][] centroidList = new double[1][3];
		double[] centroid = new double[3];
		double[] interceptCounts = new double[nVectors];
		double[] sumInterceptCounts = new double[nVectors];
		double previous = 2; // Anisotropy cannot be greater than 1, so 2 gives
		// a very high variance
		coOrdinates = new double[nVectors][3];
		if (doPlot)
			createGraph();
		Vector<Double> anisotropyHistory = new Vector<Double>();
		int s = 0;
		while (s < minIterations
				|| (s >= minIterations && s < maxIterations && variance > tolerance)) {
			s++;
			// return a single centroid within the bounds
			centroidList = gridCalculator(imp, 1, radius);
			// count intercepts at centroid
			centroid[0] = centroidList[0][0];
			centroid[1] = centroidList[0][1];
			centroid[2] = centroidList[0][2];
			IJ.showStatus("Counting intercepts at site " + s
					+ ", anisotropy = " + IJ.d2s(anisotropy, 5));
			interceptCounts = countIntercepts(imp, centroid, vectorList,
					nVectors, radius, vectorSampling);

			// add intercepts to vectors
			for (int i = 0; i < nVectors; i++) {
				sumInterceptCounts[i] += interceptCounts[i];
			}

			// work out the current mean intercept length
			double[] meanInterceptLengths = new double[nVectors];
			for (int v = 0; v < nVectors; v++) {
				if (sumInterceptCounts[v] == 0) {
					sumInterceptCounts[v] = 1;
				}
				// MIL = total vector length / number of intercepts
				meanInterceptLengths[v] = radius * (double) s
						/ sumInterceptCounts[v];
			}
			// work out coordinates of vector cloud
			for (int v = 0; v < nVectors; v++) {
				coOrdinates[v][0] = meanInterceptLengths[v] * vectorList[v][0];
				coOrdinates[v][1] = meanInterceptLengths[v] * vectorList[v][1];
				coOrdinates[v][2] = meanInterceptLengths[v] * vectorList[v][2];
			}
			// calculate principal components
			EigenvalueDecomposition E = principalComponents(coOrdinates);
			Matrix eigenValues = E.getD();
			double[][] eVal = eigenValues.getArrayCopy();
			anisotropy = 1 - eVal[0][0] / eVal[2][2];
			anisotropyHistory.add(anisotropy);
			if (doPlot)
				updateGraph(anisotropyHistory);
			variance = Math.abs(previous - anisotropy);
			previous = anisotropy;
		}
		return anisotropy;
	}

	private void createGraph() {
		double[] x = { 0 }, y = { 0 };
		Plot plot = new Plot("Anisotropy", "Repeats", "Anisotropy", x, y);
		plot.setLimits(0, 1, 0, 1);
		ImageProcessor plotIp = plot.getProcessor();
		plotImage = new ImagePlus("Anisotropy", plotIp);
		plotImage.setProcessor(null, plotIp);
		plotImage.show();
		return;
	}

	/**
	 * Generate an array of randomly-oriented 3D unit vectors
	 * 
	 * @param nVectors
	 *            number of vectors to generate
	 * @return 2D array (nVectors x 3) containing unit vectors
	 */
	private double[][] randomVectors(int nVectors) {
		double[][] randomVectors = new double[nVectors][3];
		for (int n = 0; n < nVectors; n++) {
			randomVectors[n][2] = 2 * Math.random() - 1;
			final double rho = Math.sqrt(1 - randomVectors[n][2]
					* randomVectors[n][2]);
			final double phi = Math.PI * (2 * Math.random() - 1);
			randomVectors[n][0] = rho * Math.cos(phi);
			randomVectors[n][1] = rho * Math.sin(phi);
		}
		return randomVectors;
	} /* end randomVectors */

	/*-------------------------------------------------------------------*/
	/**
	 * Generate a set of centroids padded from the image edges
	 * 
	 * @param imp
	 *            ImagePlus
	 * @param nCentroids
	 *            Number of centroids to generate
	 * @param radius
	 *            amount of padding between stack edges and centroid field
	 * @return nCentroids x 3 array of 3D coordinates
	 */
	private double[][] gridCalculator(final ImagePlus imp,
			final int nCentroids, final double radius) {
		Calibration cal = imp.getCalibration();
		final double vW = cal.pixelWidth;
		final double vH = cal.pixelHeight;
		final double vD = cal.pixelDepth;
		final int w = imp.getWidth();
		final int h = imp.getHeight();
		final int d = imp.getStackSize();
		final double stackWidth = vW * w;
		final double stackHeight = vH * h;
		final double stackDepth = vD * d;
		// strategy: n random coordinates within bounding box (easy, no bias.)
		double[][] gridCentroids = new double[nCentroids][3];
		for (int n = 0; n < nCentroids; n++) {
			gridCentroids[n][0] = Math.random()
					* (stackWidth - 2 * radius - 2 * vW) + radius;
			gridCentroids[n][1] = Math.random()
					* (stackHeight - 2 * radius - 2 * vH) + radius;
			gridCentroids[n][2] = Math.random()
					* (stackDepth - 2 * radius - 2 * vD) + radius;
		}
		// alternative: n regularly-spaced coordinates fitting
		// within bounding box
		// allegedly more efficient but could collide with periodic data

		return gridCentroids;
	} /* end gridCalculator */

	/*------------------------------------------------------*/
	/**
	 * <p>
	 * Counts the number of intercepts between each vector in a set of vectors
	 * and a binary 3D image.
	 * </p>
	 * <p>
	 * Vector length and intercept count should be summed over all vectors with
	 * same direction and sum of length divided by sum of count to give mean
	 * intercept length for each vector.
	 * </p>
	 * 
	 * 
	 * @param centroid
	 *            3-element array containing calibrated 3D centroid location
	 * @param vectorList
	 *            array containing unit vectors
	 * @param nVectors
	 *            number of vectors in each set
	 * @param radius
	 *            length of vectors
	 * @param vectorSampling
	 *            distance between tests along each vector
	 * @return 1D array containing a count of intercepts for each vector
	 */
	private double[] countIntercepts(ImagePlus imp, double[] centroid,
			double[][] vectorList, int nVectors, double radius,
			double vectorSampling) {
		Calibration cal = imp.getCalibration();
		final double vW = cal.pixelWidth;
		final double vH = cal.pixelHeight;
		final double vD = cal.pixelDepth;

		final double cX = centroid[0];
		final double cY = centroid[1];
		final double cZ = centroid[2];

		final int width = imp.getWidth();

		boolean lastPos, thisPos;

		// create a work array containing pixels +- 1 radius from centroid
		final int w = (int) Math.round(radius / vW);
		final int h = (int) Math.round(radius / vH);
		final int d = (int) Math.round(radius / vD);
		byte[] workArray = new byte[(2 * w + 1) * (2 * h + 1) * (2 * d + 1)];

		final int startCol = (int) Math.round(cX / vW) - w;
		final int endCol = (int) Math.round(cX / vW) + w;
		final int startRow = (int) Math.round(cY / vH) - h;
		final int endRow = (int) Math.round(cY / vH) + h;
		final int startSlice = (int) Math.round(cZ / vD) - d;
		final int endSlice = (int) Math.round(cZ / vD) + d;

		final ImageStack stack = imp.getImageStack();
		if (startSlice < 0 || startSlice >= stack.getSize())
			IJ.log("!!!startSlice = " + startSlice);
		if (endSlice < 0 || endSlice >= stack.getSize())
			IJ.log("!!!endSlice = " + endSlice);
		if (startRow < 0 || startRow >= stack.getHeight())
			IJ.log("!!!startRow = " + startRow);
		if (endRow < 0 || endRow >= stack.getHeight())
			IJ.log("!!!endRow = " + endRow);
		if (startCol < 0 || startCol >= stack.getWidth())
			IJ.log("!!!startCol = " + startCol);
		if (endCol < 0 || endCol >= stack.getWidth())
			IJ.log("!!!endCol = " + endCol);

		// fill the work array
		int i = 0;
		for (int s = startSlice; s <= endSlice; s++) {
			final byte[] slicePixels = (byte[]) stack.getPixels(s + 1);
			for (int r = startRow; r <= endRow; r++) {
				final int index = width * r;
				for (int c = startCol; c <= endCol; c++) {
					workArray[i] = slicePixels[index + c];
					i++;
				}
			}
		}

		// centroid position in workArray is at (w+1, h+1, d+1), subtract one
		// for starting at 0.
		final int a = (2 * w + 1);
		final int b = a * (2 * h + 1);

		final int wVz = b * d;
		final int wVy = a * h;
		final int wVx = w;
		final int centroidIndex = wVz + wVy + wVx;

		// store an intercept count for each vector
		double[] interceptCounts = new double[nVectors];

		final double radVw = -radius / vW;
		final double radVh = -radius / vH;
		final double radVd = -radius / vD;

		// loop through all vectors
		for (int v = 0; v < nVectors; v++) {
			double nIntercepts = 0;
			final double vX = vectorList[v][0];
			final double vY = vectorList[v][1];
			final double vZ = vectorList[v][2];

			// start at negative end of vector
			final int xS = (int) Math.round(radVw * vX);
			final int yS = (int) Math.round(radVh * vY);
			final int zS = (int) Math.round(radVd * vZ);

			final int startIndex = centroidIndex + b * zS + a * yS + xS;

			if (xS + w < 0 || xS + w >= (2 * w + 1) || yS + h < 0
					|| yS + h >= (2 * h + 1) || zS + d < 0
					|| zS + d >= (2 * d + 1)) {
				IJ.log("!!!!!(xS, yS, zS) (" + (xS + w) + ", " + (yS + h)
						+ ", " + (zS + d) + ")!!!!!");
				IJ.log("!!!!!startIndex is " + startIndex + "!!!!!");
			}

			if (workArray[startIndex] == 0) {
				lastPos = true;
			} else {
				lastPos = false;
			}

			final double vXvW = vX / vW;
			final double vYvH = vY / vH;
			final double vZvD = vZ / vD;

			for (double pos = -radius; pos <= radius; pos += vectorSampling) {
				// find the index of the voxel that the sample falls within
				// offset from centroid
				final int x = (int) Math.round(pos * vXvW);
				final int y = (int) Math.round(pos * vYvH);
				final int z = (int) Math.round(pos * vZvD);
				final int testIndex = centroidIndex + b * z
						+ a * y + x;
				// determine if the voxel is thresholded or not
				if (workArray[testIndex] == 0) {
					thisPos = true;
				} else {
					thisPos = false;
				}
				// if this pos is not equal to last pos then an interface is
				// counted
				if (thisPos != lastPos) {
					nIntercepts++;
				}
				// then before incrementing the for loop, set lastPos to thisPos
				lastPos = thisPos;
			}
			interceptCounts[v] = nIntercepts;
		}

		return interceptCounts;
	}/* end meanInterceptLengths */

	/*--------------------------------------------------------------------------*/
	/**
	 * Calculate the eigenvectors and eigenvalues of a set of points by the
	 * covariance method and eigendecomposition.
	 * 
	 * @param coOrdinates
	 *            n x 3 array centred on (0,0,0)
	 * @return EigenvalueDecomposition containing eigenvectors and eigenvalues
	 * 
	 */
	private EigenvalueDecomposition principalComponents(double[][] coOrdinates) {
		IJ.showStatus("Calculating eigenvalues");
		double sumX = 0, sumY = 0, sumZ = 0;
		for (int n = 0; n < coOrdinates.length; n++) {
			sumX += coOrdinates[n][0];
			sumY += coOrdinates[n][1];
			sumZ += coOrdinates[n][2];
		}
		double centX = sumX / coOrdinates.length;
		double centY = sumY / coOrdinates.length;
		double centZ = sumZ / coOrdinates.length;

		double[][] C = new double[3][3];
		double count = 0;
		for (int n = 0; n < coOrdinates.length; n++) {
			double x = coOrdinates[n][0] - centX;
			double y = coOrdinates[n][1] - centY;
			double z = coOrdinates[n][2] - centZ;
			C[0][0] += x * x;
			C[1][1] += y * y;
			C[2][2] += z * z;
			C[0][1] += x * y;
			C[0][2] += x * z;
			C[1][0] += x * y;
			C[1][2] += y * z;
			C[2][0] += x * z;
			C[2][1] += y * z;
			count += 1;
		}
		double invCount = 1 / count;
		Matrix covarianceMatrix = new Matrix(C).times(invCount);
		EigenvalueDecomposition E = new EigenvalueDecomposition(
				covarianceMatrix);
		return E;
	}/* end PrincipalComponents */

	/*---------------------------------------------------------------------*/
	/**
	 * Calculate the best-fit ellipsoid by least squares
	 * 
	 * @param coOrdinates
	 *            n x 3 array containing 3D coordinates
	 * @return 10 x 1 array containing constants for the ellipsoid equation
	 *         <i>ax</i><sup>2</sup> + <i>by</i><sup>2</sup> +
	 *         <i>cz</i><sup>2</sup> + 2<i>fyz</i> + 2<i>gxz</i> + 2<i>hxy</i> +
	 *         2<i>px</i> + 2<i>qy</i> + 2<i>rz</i> + <i>d</i> = 0
	 * 
	 * @see <p>
	 *      Qingde Li, Griffiths J (2004) Least squares ellipsoid specific
	 *      fitting. Geometric Modeling and Processing, 2004. Proceedings. pp.
	 *      335-340. <a
	 *      href="http://dx.doi.org/10.1109/GMAP.2004.1290055">doi:10.1109
	 *      /GMAP.2004.1290055</a>
	 *      </p>
	 */
	@SuppressWarnings("unused")
	private double[][] fitEllipsoid(double[][] coOrdinates) {

		IJ.showStatus("Fitting ellipsoid");
		double[][] Darray = new double[coOrdinates.length][10];
		for (int n = 0; n < coOrdinates.length; n++) {
			// populate each column of D with ten-element Xi
			double xi = coOrdinates[n][0];
			double yi = coOrdinates[n][1];
			double zi = coOrdinates[n][2];
			Darray[n][0] = xi * xi;
			Darray[n][1] = yi * yi;
			Darray[n][2] = zi * zi;
			Darray[n][3] = yi * zi * 2;
			Darray[n][4] = xi * zi * 2;
			Darray[n][5] = xi * yi * 2;
			Darray[n][6] = xi * 2;
			Darray[n][7] = yi * 2;
			Darray[n][8] = zi * 2;
			Darray[n][9] = 1;
		}
		Matrix D = new Matrix(Darray);
		Matrix S = D.times(D.transpose());
		// Define 6x6 Matrix C1 (Eq. 7)
		double[][] C1array = { { -1, 1, 1, 0, 0, 0 }, { 1, -1, 1, 0, 0, 0 },
				{ 1, 1, -1, 0, 0, 0 }, { 0, 0, 0, -4, 0, 0 },
				{ 0, 0, 0, 0, -4, 0 }, { 0, 0, 0, 0, 0, -4 } };
		Matrix C1 = new Matrix(C1array);
		// Define S11, S12, S22 from S (Eq. 11)
		Matrix S11 = S.getMatrix(0, 5, 0, 5);
		Matrix S12 = S.getMatrix(0, 5, 6, 9);
		Matrix S22 = S.getMatrix(6, 9, 6, 9);
		// Eq. 15
		EigenvalueDecomposition E = new EigenvalueDecomposition(C1.inverse()
				.times(
						S11.minus(S12.times(S22.inverse()
								.times(S12.transpose())))));
		Matrix eigenValues = E.getD();
		Matrix eigenVectors = E.getV();
		Matrix v1 = new Matrix(6, 1);
		double[][] EigenValueMatrix = eigenValues.getArray();
		double posVal = -999999999;
		for (int p = 0; p < EigenValueMatrix.length; p++) {
			if (EigenValueMatrix[p][p] > posVal) {
				posVal = EigenValueMatrix[p][p];
				v1 = eigenVectors.getMatrix(0, 5, p, p);
			}
		}
		Matrix v2 = S22.inverse().times(S12.transpose()).times(v1).times(-1);
		Matrix v = new Matrix(10, 1);
		int[] c = { 0 };
		v.setMatrix(0, 5, c, v1);
		v.setMatrix(6, 9, c, v2);
		double[][] ellipsoid = v.getArrayCopy();
		IJ.log("Ellipse equation: " + ellipsoid[0][0] + " x^2 + "
				+ ellipsoid[1][0] + " y^2 + " + ellipsoid[2][0] + " z^2 + "
				+ ellipsoid[3][0] + " 2yz + " + ellipsoid[4][0] + " 2xz + "
				+ ellipsoid[5][0] + " 2xy + " + ellipsoid[6][0] + " 2x + "
				+ ellipsoid[7][0] + " 2y + " + ellipsoid[8][0] + " 2z + "
				+ ellipsoid[9][0] + " = 0");

		return ellipsoid;
	} /* end fitEllipsoid */

	private void printMatrix(Matrix matrix) {
		int nCols = matrix.getColumnDimension();
		int nRows = matrix.getRowDimension();
		double[][] eVal = matrix.getArrayCopy();
		for (int r = 0; r < nRows; r++) {
			String row = "||";
			for (int c = 0; c < nCols; c++) {
				row = row + eVal[r][c] + "|";
			}
			row = row + "|";
			IJ.log(row);
		}
	}

	private void updateGraph(Vector<Double> anisotropyHistory) {
		double[] yVariables = new double[anisotropyHistory.size()];
		double[] xVariables = new double[anisotropyHistory.size()];
		Enumeration<Double> e = anisotropyHistory.elements();
		int i = 0;
		while (e.hasMoreElements()) {
			yVariables[i] = e.nextElement();
			xVariables[i] = (double) i;
			i++;
		}
		Plot plot = new Plot("Anisotropy", "Number of repeats", "Anisotropy",
				xVariables, yVariables);
		plot.addPoints(xVariables, yVariables, Plot.X);
		plot.setLimits(0, anisotropyHistory.size(), 0, 1);
		ImageProcessor plotIp = plot.getProcessor();
		plotImage.setProcessor(null, plotIp);
	}

	/*-------------------------------------------------------------*/
	/**
	 * Plot a set of 3D coordinates in Benjamin Schmidt's ImageJ 3D Viewer
	 * 
	 * @param coOrdinates
	 *            float[][] n x 3 array of 3D (x,y,z) coordinates
	 * @param name
	 *            String name of the dataset
	 * 
	 */
	private void plotPoints3D(double[][] coOrdinates, String name) {
		// Create a CustomMesh from the coordinates
		List<Point3f> mesh = new ArrayList<Point3f>();
		for (int i = 0; i < coOrdinates.length; i++) {
			mesh.add(new Point3f((float) coOrdinates[i][0],
					(float) coOrdinates[i][1], (float) coOrdinates[i][2]));
		}
		// add the other ends of the vectors
		List<Point3f> mesh2 = new ArrayList<Point3f>();
		for (int i = 0; i < coOrdinates.length; i++) {
			mesh2.add(new Point3f(-(float) coOrdinates[i][0],
					-(float) coOrdinates[i][1], -(float) coOrdinates[i][2]));
		}

		CustomPointMesh cm = new CustomPointMesh(mesh);
		CustomPointMesh cm2 = new CustomPointMesh(mesh2);

		// Create a universe and show it
		Image3DUniverse univ = new Image3DUniverse();
		univ.show();

		// Add the mesh
		Content c = univ.addCustomMesh(cm, "heads");
		Content c2 = univ.addCustomMesh(cm2, "tails");
		Color3f green = new Color3f(0.0f, 0.5f, 0.0f);
		Color3f red = new Color3f(0.5f, 0.0f, 0.0f);
		c.getColor();
		c2.getColor();
		cm.setColor(green);
		cm2.setColor(red);
		cm.setPointSize(1);

		// Have a look at the source code of CustomPointMesh
		// for changing point size and anti-aliasing
	}/* end plotPoints3D */
}