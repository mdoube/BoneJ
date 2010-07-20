package org.doube.bonej;

/**
 *Anisotropy_ plugin for ImageJ
 *Copyright 2009 2010 Michael Doube 
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
import ij.macro.Interpreter;
import ij.measure.Calibration;

// for 3D plotting of coordinates
import javax.vecmath.Point3f;
import javax.vecmath.Color3f;
import customnode.CustomPointMesh;

import java.util.Enumeration;
import java.util.List;
import java.util.ListIterator;
import java.util.ArrayList;
import java.util.Vector;

import ij3d.Image3DUniverse;
import ij3d.Content;

import org.doube.geometry.FitEllipsoid;
import org.doube.jama.Matrix;
import org.doube.jama.EigenvalueDecomposition;
import org.doube.util.ImageCheck;
import org.doube.util.ResultInserter;

/**
 * <p>
 * Calculate degree of anisotropy using the mean intercept length method.
 * Samples vectors in a sphere; the location of the sphere is random within the
 * stack. New sampling spheres are generated until the coefficient of variation
 * of the anisotropy result reduces to a user-set threshold or a maximum number
 * of spheres has been sampled.
 * </p>
 * 
 * @see <p>
 *      Harrigan TP, Mann RW (1984) Characterization of microstructural
 *      anisotropy in orthotropic materials using a second rank tensor. J Mater
 *      Sci 19: 761-767. <a
 *      href="http://dx.doi.org/10.1007/BF00540446">doi:10.1007/BF00540446</a>.
 *      </p>
 * @author Michael Doube
 * 
 */
public class Anisotropy implements PlugIn {

	public void run(String arg) {
		if (!ImageCheck.checkEnvironment()) {
			return;
		}
		final ImagePlus imp = IJ.getImage();
		if (null == imp) {
			IJ.noImage();
			return;
		}
		ImageCheck ic = new ImageCheck();
		if (!ic.isBinary(imp)) {
			IJ.error("8-bit binary (black and white only) image required.");
			return;
		}
		if (!ic.isMultiSlice(imp) || imp.getStackSize() < 5) {
			IJ.error("Stack with at least 5 slices required");
			return;
		}

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
		gd.addNumericField("Tolerance", 0.005, 4, 6, "");
		gd.addCheckbox("Show_Plot", true);
		gd.addCheckbox("3D_Result", false);
		gd.addCheckbox("Align to fabric tensor", false);
		gd.addHelp("http://bonej.org/anisotropy");
		gd.showDialog();
		if (gd.wasCanceled()) {
			return;
		}
		final boolean doAutoMode = gd.getNextBoolean();
		final int nVectors = (int) gd.getNextNumber();
		final int minSpheres = (int) gd.getNextNumber();
		final int maxSpheres = (int) gd.getNextNumber();
		final double tolerance = gd.getNextNumber();
		final boolean doPlot = gd.getNextBoolean();
		final boolean do3DResult = gd.getNextBoolean();
		final boolean doAlign = gd.getNextBoolean();

		Object[] result = new Object[3];
		if (doAutoMode)
			result = runToStableResult(imp, minSpheres, maxSpheres, nVectors,
					radius, vectorSampling, tolerance, doPlot);
		else
			result = runToStableResult(imp, minSpheres, minSpheres, nVectors,
					radius, vectorSampling, tolerance, doPlot);

		double[] anisotropy = (double[]) result[0];
		double[][] coOrdinates = (double[][]) result[1];

		ResultInserter ri = ResultInserter.getInstance();
		ri.setResultInRow(imp, "DA", anisotropy[0]);
		ri.updateTable();

		if (do3DResult) {
			plotPoints3D(coOrdinates, "Intercept Lengths");
		}

		if (doAlign) {
			EigenvalueDecomposition E = (EigenvalueDecomposition) result[2];
			Moments m = new Moments();
			ImagePlus alignedImp = m.alignImage(imp, E.getV(), false, 1, d,
					128, 255, 0, 1);
			alignedImp.show();
		}

		return;
	}

	/**
	 * Calculate degree of anisotropy for a binary stack, running until a stable
	 * result is achieved, or the maximum number of iterations occurs.
	 * 
	 * @param imp
	 *            ImagePlus input. A binary stack is required.
	 * @param minSpheres
	 *            minimum number of iterations
	 * @param maxSpheres
	 *            maximum number of iterations
	 * @param nVectors
	 *            number of vectors in the sampling sphere
	 * @param radius
	 *            radius of the sampling sphere
	 * @param vectorSampling
	 *            distance between samples along the sampling vector. Set to 2.3
	 *            * maximum voxel dimension for a safe and efficient sampling
	 *            increment
	 * @param tolerance
	 *            coefficient of variation of results at which we accept result
	 *            is stable
	 * @param doPlot
	 *            set to true if you want to see a plot of anisotropy versus
	 *            number of repeats, updated in real time
	 * @return Object array containing degree of anisotropy, coordinates of rose
	 *         plot and Eigenvalue decomposition (fabric tensor)
	 */
	public Object[] runToStableResult(ImagePlus imp, int minSpheres,
			int maxSpheres, int nVectors, double radius, double vectorSampling,
			double tolerance, boolean doPlot) {
		final int minIterations = minSpheres;
		final int maxIterations = maxSpheres;
		final double[][] vectorList = randomVectors(nVectors);
		double variance = Double.NaN;
		double anisotropy = Double.NaN;
		double[][] centroidList = new double[1][3];
		double[] centroid = new double[3];
		double[] interceptCounts = new double[nVectors];
		double[] sumInterceptCounts = new double[nVectors];
		double[][] coOrdinates = new double[nVectors][3];
		ImagePlus plotImage = new ImagePlus();
		// never show a plot if this method is being called from a batch mode
		// macro
		if (Interpreter.isBatchMode())
			doPlot = false;
		if (doPlot) {
			plotImage = createGraph(imp.getTitle());
			plotImage.show();
		}
		Vector<Double> anisotropyHistory = new Vector<Double>();
		Vector<Double> errorHistory = new Vector<Double>();
		double[][] emptyArray = new double[3][3];
		Matrix emptyMatrix = new Matrix(emptyArray);
		EigenvalueDecomposition E = new EigenvalueDecomposition(emptyMatrix);
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
					+ ", anisotropy = " + IJ.d2s(anisotropy, 5) + ", CV = "
					+ IJ.d2s(variance, 3));
			interceptCounts = countIntercepts(imp, centroid, vectorList,
					nVectors, radius, vectorSampling);

			// add intercepts to vectors
			for (int i = 0; i < nVectors; i++) {
				sumInterceptCounts[i] += interceptCounts[i];
			}

			// work out the current mean intercept length
			double[] meanInterceptLengths = new double[nVectors];
			final double probeLength = radius * (double) s;
			for (int v = 0; v < nVectors; v++) {
				if (sumInterceptCounts[v] == 0) {
					sumInterceptCounts[v] = 1;
				}
				// MIL = total vector length / number of intercepts
				meanInterceptLengths[v] = probeLength / sumInterceptCounts[v];
			}

			// work out coordinates of vector cloud
			for (int v = 0; v < nVectors; v++) {
				final double milV = meanInterceptLengths[v];
				coOrdinates[v][0] = milV * vectorList[v][0];
				coOrdinates[v][1] = milV * vectorList[v][1];
				coOrdinates[v][2] = milV * vectorList[v][2];
			}
			Object[] result = harriganMann(coOrdinates);
			anisotropy = ((double[]) result[0])[0];
			anisotropyHistory.add(anisotropy);
			E = (EigenvalueDecomposition) result[1];

			variance = getVariance(anisotropyHistory, minIterations);

			if (variance + anisotropy > 1 || anisotropy - variance < 0) {
				variance = Math.max(Math.min(1 - anisotropy, anisotropy),
						tolerance);
			}

			errorHistory.add(variance);
			if (doPlot)
				updateGraph(plotImage, anisotropyHistory, errorHistory);
		}
		double[] da = { anisotropy };
		Object[] result = { da, coOrdinates, E };
		return result;
	}

	/**
	 * Create a graph for plotting anisotropy results
	 * 
	 * @param id
	 *            Identification string for Anisotropy plot
	 * @return ImagePlus for drawing a plot on
	 */
	private ImagePlus createGraph(String id) {
		double[] x = { 0 }, y = { 0 };
		Plot plot = new Plot("Anisotropy", "Repeats", "Anisotropy", x, y);
		plot.setLimits(0, 1, 0, 1);
		ImageProcessor plotIp = plot.getProcessor();
		ImagePlus plotImage = new ImagePlus("Anisotropy of " + id, plotIp);
		plotImage.setProcessor(null, plotIp);
		return plotImage;
	}

	/**
	 * 
	 * Calculate coefficient of variation of last n results
	 * 
	 * @param anisotropyHistory
	 *            list of anisotropy results, one result per iteration
	 * 
	 * @return coefficient of variation, which is standard deviation / mean.
	 */

	private double getVariance(Vector<Double> anisotropyHistory, final int n) {
		ListIterator<Double> iter = anisotropyHistory
				.listIterator(anisotropyHistory.size());

		double sum = 0;

		double sumSquares = 0;

		int count = 0;

		while (iter.hasPrevious()) {
			final double value = iter.previous();
			sum += value;
			count++;
			if (count >= n)
				break;
		}

		final double mean = sum / n;
		ListIterator<Double> itr = anisotropyHistory
				.listIterator(anisotropyHistory.size());

		count = 0;
		while (itr.hasPrevious()) {
			final double value = itr.previous();
			final double a = value - mean;
			sumSquares += a * a;
			count++;
			if (count >= n)
				break;
		}

		double stDev = Math.sqrt(sumSquares / n);
		double coeffVariation = stDev / mean;
		return coeffVariation;
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
			final double z = 2 * Math.random() - 1;
			final double rho = Math.sqrt(1 - z * z);
			final double phi = Math.PI * (2 * Math.random() - 1);
			randomVectors[n][0] = rho * Math.cos(phi);
			randomVectors[n][1] = rho * Math.sin(phi);
			randomVectors[n][2] = z;
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

		// loop through all vectors
		// start multithreading here - each thread samples a set of vectors
		int nThreads = Runtime.getRuntime().availableProcessors();
		InterceptThread[] it = new InterceptThread[nThreads];
		for (int thread = 0; thread < nThreads; thread++) {
			it[thread] = new InterceptThread(thread, nThreads, nVectors,
					centroidIndex, a, b, radius, vW, vH, vD, vectorSampling,
					interceptCounts, vectorList, workArray);
			it[thread].start();
		}
		try {
			for (int thread = 0; thread < nThreads; thread++) {
				it[thread].join();
			}
		} catch (InterruptedException ie) {
			IJ.error("A thread was interrupted.");
		}
		return interceptCounts;
	}/* end meanInterceptLengths */

	class InterceptThread extends Thread {
		int thread, nThreads, nVectors, centroidIndex, a, b;
		double radius, vW, vH, vD, vectorSampling;
		double[] interceptCounts;
		double[][] vectorList;
		byte[] workArray;

		public InterceptThread(int thread, int nThreads, int nVectors,
				int centroidIndex, int a, int b, double radius, double vW,
				double vH, double vD, double vectorSampling,
				double[] interceptCounts, double[][] vectorList,
				byte[] workArray) {
			this.thread = thread;
			this.nThreads = nThreads;
			this.nVectors = nVectors;
			this.centroidIndex = centroidIndex;
			this.a = a;
			this.b = b;
			this.radius = radius;
			this.vW = vW;
			this.vH = vH;
			this.vD = vD;
			this.vectorSampling = vectorSampling;
			this.interceptCounts = interceptCounts;
			this.vectorList = vectorList;
			this.workArray = workArray;
		}

		public void run() {
			final double radVw = -radius * vW;
			final double radVh = -radius * vH;
			final double radVd = -radius * vD;
			for (int v = this.thread; v < this.nVectors; v += this.nThreads) {
				double nIntercepts = 0;
				final double vX = vectorList[v][0];
				final double vY = vectorList[v][1];
				final double vZ = vectorList[v][2];

				// start at negative end of vector
				final int xS = (int) Math.round(radVw * vX);
				final int yS = (int) Math.round(radVh * vY);
				final int zS = (int) Math.round(radVd * vZ);

				final int startIndex = centroidIndex + b * zS + a * yS + xS;
				boolean lastPos, thisPos;
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
					final int testIndex = centroidIndex + b * z + a * y + x;
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
					// then before incrementing the for loop, set lastPos to
					// thisPos
					lastPos = thisPos;
				}
				interceptCounts[v] = nIntercepts;
			}
		}
	}

	/*--------------------------------------------------------------------------*/
	/**
	 * Draw on plotImage the data in anisotropyHistory with error bars from
	 * errorHistory
	 * 
	 * @param plotImage
	 *            the graph image
	 * @param anisotropyHistory
	 *            all anisotropy results, 1 for each iteration
	 * @param errorHistory
	 *            all error results, 1 for each iteration
	 */
	private void updateGraph(ImagePlus plotImage,
			Vector<Double> anisotropyHistory, Vector<Double> errorHistory) {
		double[] yVariables = new double[anisotropyHistory.size()];
		double[] xVariables = new double[anisotropyHistory.size()];
		Enumeration<Double> e = anisotropyHistory.elements();
		int i = 0;
		while (e.hasMoreElements()) {
			yVariables[i] = e.nextElement();
			xVariables[i] = (double) i;
			i++;
		}
		Enumeration<Double> f = errorHistory.elements();
		i = 0;
		double[] errorBars = new double[errorHistory.size()];
		while (f.hasMoreElements()) {
			errorBars[i] = f.nextElement();
			i++;
		}
		Plot plot = new Plot("Anisotropy", "Number of repeats", "Anisotropy",
				xVariables, yVariables);
		plot.addPoints(xVariables, yVariables, Plot.X);
		plot.addErrorBars(errorBars);
		plot.setLimits(0, anisotropyHistory.size(), 0, 1);
		ImageProcessor plotIp = plot.getProcessor();
		plotImage.setProcessor(null, plotIp);
		return;
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
		final int nPoints = coOrdinates.length;
		// Create a CustomMesh from the coordinates
		List<Point3f> mesh = new ArrayList<Point3f>();
		for (int i = 0; i < nPoints; i++) {
			mesh.add(new Point3f((float) coOrdinates[i][0],
					(float) coOrdinates[i][1], (float) coOrdinates[i][2]));
		}
		// add the other ends of the vectors
		List<Point3f> mesh2 = new ArrayList<Point3f>();
		for (int i = 0; i < nPoints; i++) {
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

	/**
	 * Calculate degree of anisotropy according to Harrigan and Mann's
	 * ellipsoidal tensor method
	 * 
	 * @param coOrdinates
	 * @return Object[] containing degree of anisotropy and the
	 *         eigendecomposition
	 */
	private Object[] harriganMann(double[][] coOrdinates) {
		Object[] ellipsoid = new Object[6];
		double da = 0;
		try {
			ellipsoid = FitEllipsoid.yuryPetrov(coOrdinates);
		} catch (RuntimeException re) {
			da = Math.random();
		}
		double[] coEf = (double[]) ellipsoid[3];
		double[][] tensor = { { coEf[0], coEf[3], coEf[4] },
				{ coEf[3], coEf[1], coEf[5] }, { coEf[4], coEf[5], coEf[2] } };
		Matrix M = new Matrix(tensor);
		EigenvalueDecomposition E = M.eig();
		Matrix EigenVal = E.getD();
		double[] diag = EigenVal.diag().getColumnPackedCopy();
		da = 1 - diag[0] / diag[2];
		if (da > 1)
			da = 1;
		else if (da < 0)
			da = 0;
		double[] anisotropy = { da };
		Object[] result = { anisotropy, E };
		return result;
	}

}