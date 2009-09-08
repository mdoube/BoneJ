package org.doube.bonej;

/**
 * FitCircle Java class for fitting circles to 2D coordinate data
 * 
 * Copyright 2009 Michael Doube 
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

import java.util.Arrays;

import ij.IJ;

import Jama.Matrix;
import Jama.EigenvalueDecomposition;
import Jama.SingularValueDecomposition;

/**
 * Methods for fitting circles to coordinates
 * 
 * @author Michael Doube, ported from Nikolai Chernov's MATLAB scripts
 * @see <p>
 *      Al-Sharadqha & Chernov (2009) <a
 *      href="http://dx.doi.org/10.1214/09-EJS419"> Error analysis for circle
 *      fitting algorithms</a>. Electronic Journal of Statistics 3, pp. 886-911<br/>
 *      <a href="http://www.math.uab.edu/~chernov/cl/MATLABcircle.html"
 *      >http://www.math.uab.edu/~chernov/cl/MATLABcircle.html</a>
 *      </p>
 * 
 */
public class FitCircle {

    /**
     * Kåsa fit
     * 
     * 
     */
    public double[] kasaFit(double[][] points) {
	int nPoints = points.length;
	double[][] z = new double[nPoints][1];
	double[][] xy1 = new double[nPoints][3];
	for (int n = 0; n < nPoints; n++) {
	    z[n][0] = points[n][0] * points[n][0] + points[n][1] * points[n][1];
	    xy1[n][0] = points[n][0];
	    xy1[n][1] = points[n][1];
	    xy1[n][2] = 1;
	}

	Matrix XY1 = new Matrix(xy1);
	Matrix Z = new Matrix(z);
	Matrix P = (XY1.inverse()).times(Z); // no direct left divide in Jama!
	double[] centreRadius = new double[3];

	double p0 = P.get(0, 0);
	double p1 = P.get(1, 0);
	double p2 = P.get(2, 0);
	centreRadius[0] = p0 / 2;
	centreRadius[1] = p1 / 2;
	centreRadius[2] = Math.sqrt((p0 * p0 + p1 * p1) / 4 + p2);
	return centreRadius;
    }

    /**
     * Pratt method (Newton style)
     * 
     * Can't handle noise
     * 
     * @param points
     * @return
     */
    public double[] prattNewton(double[][] points) {
	int nPoints = points.length;
	double[] centroid = getCentroid(points);
	double Mxx = 0, Myy = 0, Mxy = 0, Mxz = 0, Myz = 0, Mzz = 0;

	for (int i = 0; i < nPoints; i++) {
	    double Xi = points[i][0] - centroid[0];
	    double Yi = points[i][1] - centroid[1];
	    double Zi = Xi * Xi + Yi * Yi;
	    Mxy += Xi * Yi;
	    Mxx += Xi * Xi;
	    Myy += Yi * Yi;
	    Mxz += Xi * Zi;
	    Myz += Yi * Zi;
	    Mzz += Zi * Zi;
	}
	Mxx /= nPoints;
	Myy /= nPoints;
	Mxy /= nPoints;
	Mxz /= nPoints;
	Myz /= nPoints;
	Mzz /= nPoints;

	double Mz = Mxx + Myy;
	double Cov_xy = Mxx * Myy - Mxy * Mxy;
	double Mxz2 = Mxz * Mxz;
	double Myz2 = Myz * Myz;

	double A2 = 4 * Cov_xy - 3 * Mz * Mz - Mzz;
	double A1 = Mzz * Mz + 4 * Cov_xy * Mz - Mxz2 - Myz2 - Mz * Mz * Mz;
	double A0 = Mxz2 * Myy + Myz2 * Mxx - Mzz * Cov_xy - 2 * Mxz * Myz
		* Mxy + Mz * Mz * Cov_xy;
	double A22 = A2 + A2;

	double epsilon = 1e-12;
	double ynew = 1e+20;
	int IterMax = 20;
	double xnew = 0;
	for (int iter = 0; iter < IterMax; iter++) {
	    double yold = ynew;
	    ynew = A0 + xnew * (A1 + xnew * (A2 + 4 * xnew * xnew));
	    if (Math.abs(ynew) > Math.abs(yold)) {
		IJ.log("Newton-Pratt goes wrong direction: |ynew| > |yold|");
		xnew = 0;
		break;
	    }
	    double Dy = A1 + xnew * (A22 + 16 * xnew * xnew);
	    double xold = xnew;
	    xnew = xold - ynew / Dy;
	    if (Math.abs((xnew - xold) / xnew) < epsilon) {
		break;
	    }
	    if (iter >= IterMax) {
		IJ.log("Newton-Pratt will not converge");
		xnew = 0;
	    }
	    if (xnew < 0.) {
		IJ.log("Newton-Pratt negative root:  x= " + xnew);
		xnew = 0;
	    }
	}
	double det = xnew * xnew - xnew * Mz + Cov_xy;
	double x = Mxz * (Myy - xnew) - Myz * Mxy / (det * 2);
	double y = Myz * (Mxx - xnew) - Mxz * Mxy / (det * 2);
	double r = Math.sqrt(x * x + y * y + Mz + 2 * xnew);

	double[] centreRadius = { x + centroid[0], y + centroid[1], r };
	return centreRadius;
    }

    /**
     * Pratt method (SVD style)
     * 
     * 
     * @param double[n][2] containing n (<i>x</i>, <i>y</i>) coordinates
     * @return double[] containing (<i>x</i>, <i>y</i>) centre and radius
     */
    public double[] prattSVD(double[][] points) {
	int nPoints = points.length;
	double[] centroid = getCentroid(points);
	double[][] xyXY1 = new double[nPoints][4];
	for (int i = 0; i < nPoints; i++) {
	    double x = points[i][0] - centroid[0];
	    double y = points[i][1] - centroid[1];
	    xyXY1[i][0] = x * x + y * y;
	    xyXY1[i][1] = x;
	    xyXY1[i][2] = y;
	    xyXY1[i][3] = 1;
	}

	Matrix XYXY1 = new Matrix(xyXY1);
	SingularValueDecomposition svd = new SingularValueDecomposition(XYXY1);
	Matrix S = svd.getS();
	Matrix V = svd.getV();
	Matrix A;
	Matrix W;
	if (S.get(3, 3) / S.get(0, 0) < 1e-12) {
	    A = V.getMatrix(0, V.getRowDimension() - 1, 3, 3);
	    IJ.log("Pratt singular case");
	} else {
	    W = V.times(S);
	    double[][] bInv = { { 0, 0, 0, -0.5 }, { 0, 1, 0, 0 },
		    { 0, 0, 1, 0 }, { -0.5, 0, 0, 0 } };
	    Matrix Binv = new Matrix(bInv);
	    EigenvalueDecomposition ed = new EigenvalueDecomposition((W
		    .transpose()).times(Binv.times(W)));
	    Matrix D = ed.getD();
	    Matrix E = ed.getV();
	    int col = getNthSmallestCol(D, 2);
	    A = E.getMatrix(0, E.getRowDimension() - 1, col, col);
	    for (int i = 0; i < 4; i++) {
		S.set(i, i, 1 / S.get(i, i));
	    }
	    A = V.times(S.times(A));
	}
	double a0 = A.get(0, 0);
	double a1 = A.get(1, 0);
	double a2 = A.get(2, 0);
	double a3 = A.get(3, 0);
	double[] centreRadius = new double[3];
	centreRadius[0] = -(a1 / a0) / 2 + centroid[0];
	centreRadius[1] = -(a2 / a0) / 2 + centroid[1];
	centreRadius[2] = (Math.sqrt(a1 * a1 + a2 * a2 - 4 * a0 * a3) / Math
		.abs(a0)) / 2;
	return centreRadius;
    }

    /**
     * Taubin method (Newton Style)
     * 
     * @param double[n][2] containing n (<i>x</i>, <i>y</i>) coordinates
     * @return double[] containing (<i>x</i>, <i>y</i>) centre and radius
     */
    public double[] taubinNewton(double[][] points) {
	int nPoints = points.length;
	double[] centroid = getCentroid(points);
	double Mxx = 0, Myy = 0, Mxy = 0, Mxz = 0, Myz = 0, Mzz = 0;
	for (int i = 0; i < nPoints; i++) {
	    double Xi = points[i][0] - centroid[0];
	    double Yi = points[i][1] - centroid[1];
	    double Zi = Xi * Xi + Yi * Yi;
	    Mxy = Mxy + Xi * Yi;
	    Mxx = Mxx + Xi * Xi;
	    Myy = Myy + Yi * Yi;
	    Mxz = Mxz + Xi * Zi;
	    Myz = Myz + Yi * Zi;
	    Mzz = Mzz + Zi * Zi;

	}
	Mxx = Mxx / nPoints;
	Myy = Myy / nPoints;
	Mxy = Mxy / nPoints;
	Mxz = Mxz / nPoints;
	Myz = Myz / nPoints;
	Mzz = Mzz / nPoints;

	double Mz = Mxx + Myy;
	double Cov_xy = Mxx * Myy - Mxy * Mxy;
	double A3 = 4 * Mz;
	double A2 = -3 * Mz * Mz - Mzz;
	double A1 = Mzz * Mz + 4 * Cov_xy * Mz - Mxz * Mxz - Myz * Myz - Mz
		* Mz * Mz;
	double A0 = Mxz * Mxz * Myy + Myz * Myz * Mxx - Mzz * Cov_xy - 2 * Mxz
		* Myz * Mxy + Mz * Mz * Cov_xy;
	double A22 = A2 + A2;
	double A33 = A3 + A3 + A3;

	double xnew = 0;
	double ynew = 1e+20;
	double epsilon = 1e-12;
	double iterMax = 20;

	for (int iter = 0; iter < iterMax; iter++) {
	    double yold = ynew;
	    ynew = A0 + xnew * (A1 + xnew * (A2 + xnew * A3));
	    if (Math.abs(ynew) > Math.abs(yold)) {
		IJ.log("Newton-Taubin goes wrong direction: |ynew| > |yold|");
		xnew = 0;
		break;
	    }
	    double Dy = A1 + xnew * (A22 + xnew * A33);
	    double xold = xnew;
	    xnew = xold - ynew / Dy;
	    if (Math.abs((xnew - xold) / xnew) < epsilon) {
		break;
	    }
	    if (iter >= iterMax) {
		IJ.log("Newton-Taubin will not converge");
		xnew = 0;
	    }
	    if (xnew < 0.) {
		IJ.log("Newton-Taubin negative root: x = " + xnew);
		xnew = 0;
	    }
	}
	double[] centreRadius = new double[3];
	double det = xnew * xnew - xnew * Mz + Cov_xy;
	double x = (Mxz * (Myy - xnew) - Myz * Mxy) / (det / 2);
	double y = (Myz * (Mxx - xnew) - Mxz * Mxy) / (det / 2);
	centreRadius[0] = x + centroid[0];
	centreRadius[1] = y + centroid[1];
	centreRadius[2] = Math.sqrt(x * x + y * y + Mz);

	return centreRadius;
    }

    /**
     * Taubin method, SVD version
     * 
     * @param points
     * @return double[] containing (<i>x</i>, <i>y</i>) centre and radius
     */
    public double[] taubinSVD(double[][] points) {
	int nPoints = points.length;
	double[] centroid = getCentroid(points);

	double[][] zxy = new double[nPoints][3];
	double[] z = new double[nPoints];
	double sumZ = 0;
	for (int n = 0; n < nPoints; n++) {
	    double x = points[n][0] - centroid[0];
	    double y = points[n][1] - centroid[1];
	    zxy[n][1] = x;
	    zxy[n][2] = y;
	    z[n] = x * x + y * y;
	    sumZ += z[n];
	}
	double meanZ = sumZ / nPoints;
	double sqrtMeanZ2 = 2 * Math.sqrt(meanZ);
	for (int n = 0; n < nPoints; n++) {
	    zxy[n][0] = (z[n] - meanZ) / sqrtMeanZ2;
	}

	Matrix ZXY = new Matrix(zxy);
	SingularValueDecomposition svd = new SingularValueDecomposition(ZXY);
	Matrix V = svd.getV();
	Matrix A = V.getMatrix(0, V.getRowDimension() - 1, 2, 2);
	double a1 = A.get(0, 0) / sqrtMeanZ2;
	A.set(0, 0, a1);
	double[][] a = A.getArray();
	double[] b = new double[a.length + 1];
	for (int i = 0; i < a.length; i++) {
	    b[i] = a[i][0];
	}
	b[b.length - 1] = -meanZ * a1;
	double[] centreRadius = new double[3];

	centreRadius[0] = b[1] / (a1 * 2) + centroid[0];
	centreRadius[1] = b[2] / (a1 * 2) + centroid[1];
	centreRadius[2] = Math
		.sqrt(b[1] * b[1] + b[2] * b[2] - 4 * b[0] * b[3])
		/ (Math.abs(b[0]) * 2);
	return centreRadius;
    }

    /**
     * Chernov's non-biased Hyper algebraic method. Simple version.
     * 
     * @see <p>
     *      <a
     *      href="http://www.math.uab.edu/~chernov/cl/HyperSVD.m">http://www.math
     *      .uab.edu/~chernov/cl/HyperSVD.m</a>
     *      </p>
     * 
     * @param points
     *            double[n][2] containing n (<i>x</i>, <i>y</i>) coordinates
     * @return 3-element double[] containing (<i>x</i>, <i>y</i>) centre and
     *         circle radius
     */
    public double[] hyperSimple(double[][] points) {
	int nPoints = points.length;

	double[][] zxy1 = new double[nPoints][4];
	double xSum = 0, ySum = 0, zSum = 0;
	for (int n = 0; n < nPoints; n++) {
	    double x = points[n][0];
	    double y = points[n][1];
	    double z = x * x + y * y;
	    zxy1[n][0] = z;
	    zxy1[n][1] = x;
	    zxy1[n][2] = y;
	    zxy1[n][3] = 1;

	    xSum += x;
	    ySum += y;
	    zSum += z;

	}
	double s0 = zSum / nPoints;
	double s1 = xSum / nPoints;
	double s2 = ySum / nPoints;

	Matrix ZXY1 = new Matrix(zxy1);
	Matrix M = (ZXY1.transpose()).times(ZXY1);

	double[][] n = { { 8 * s0, 4 * s1, 4 * s2, 2 }, { 4 * s1, 1, 0, 0 },
		{ 4 * s2, 0, 1, 0 }, { 2, 0, 0, 0 } };
	Matrix N = new Matrix(n);
	Matrix NM = (N.inverse()).times(M);
	
	EigenvalueDecomposition ED = new EigenvalueDecomposition(NM);
	Matrix E = ED.getV();
	Matrix D = ED.getD();

	int col = getNthSmallestCol(D, 2);

	Matrix A = E.getMatrix(0, E.getRowDimension() - 1, col, col);

	double[] centreRadius = new double[3];

	centreRadius[0] = -1 * (A.get(1, 0) / A.get(0, 0)) / 2;
	centreRadius[1] = -1 * (A.get(2, 0) / A.get(0, 0)) / 2;

	double[][] a = A.getArray();
	centreRadius[2] = Math.sqrt(a[1][0] * a[1][0] + a[2][0] * a[2][0] - 4
		* a[0][0] * a[3][0])
		/ Math.abs(a[0][0]) / 2;

	return centreRadius;
    }

    /**
     * Chernov's non-biased Hyper algebraic method. Stability optimised version.
     * 
     * @see <p>
     *      <a
     *      href="http://www.math.uab.edu/~chernov/cl/HyperSVD.m">http://www.math
     *      .uab.edu/~chernov/cl/HyperSVD.m</a>
     *      </p>
     * 
     * @param points
     *            double[n][2] containing n (<i>x</i>, <i>y</i>) coordinates
     * @return 3-element double[] containing (<i>x</i>, <i>y</i>) centre and
     *         circle radius
     */
    public double[] hyperStable(double[][] points) {
	int nPoints = points.length;
	double[] centroid = getCentroid(points);

	double sumZ = 0;
	double[][] zxy1 = new double[nPoints][4];
	// centre data and assign vector values
	for (int n = 0; n < nPoints; n++) {
	    double x = points[n][0] - centroid[0];
	    double y = points[n][1] - centroid[1];
	    double z = x * x + y * y;
	    sumZ += z;
	    zxy1[n][0] = z;
	    zxy1[n][1] = x;
	    zxy1[n][2] = y;
	    zxy1[n][3] = 1;
	}

	Matrix ZXY1 = new Matrix(zxy1);
	SingularValueDecomposition svd = new SingularValueDecomposition(ZXY1);

	Matrix S = svd.getS();
	Matrix V = svd.getV();

	Matrix A;

	// singular case
	if (S.get(3, 3) / S.get(0, 0) < 1e-12) {
	    A = V.getMatrix(0, V.getRowDimension() - 1, 3, 3);
	} else {
	    // regular case
	    Matrix Y = V.times(S.times(V.transpose()));

	    double[][] bInv = { { 0, 0, 0, 0.5 }, { 0, 1, 0, 0 },
		    { 0, 0, 1, 0 }, { 0.5, 0, 0, -2 * sumZ / nPoints } };
	    Matrix Binv = new Matrix(bInv);
	    EigenvalueDecomposition ED = new EigenvalueDecomposition((Y
		    .transpose()).times(Binv.times(Y)));
	    Matrix D = ED.getD(); // eigenvalues
	    Matrix E = ED.getV(); // eigenvectors

	    int col = getNthSmallestCol(D, 2);

	    A = E.getMatrix(0, E.getRowDimension() - 1, col, col);

	    for (int i = 0; i < 4; i++) {
		S.set(i, i, 1 / S.get(i, i));
	    }
	    A = V.times(S.times((V.transpose()).times(A)));
	}

	double a0 = A.get(0, 0);
	double a1 = A.get(1, 0);
	double a2 = A.get(2, 0);
	double a3 = A.get(3, 0);

	double[] centreRadius = new double[3];
	centreRadius[0] = -(a1 / a0) / 2 + centroid[0];
	centreRadius[1] = -(a2 / a0) / 2 + centroid[1];
	centreRadius[2] = (Math.sqrt(a1 * a1 + a2 * a2 - 4 * a0 * a3) / Math
		.abs(a0)) / 2;
	return centreRadius;
    }

    private double[] getCentroid(double[][] points) {
	double[] centroid = new double[2];
	double sumX = 0;
	double sumY = 0;
	int nPoints = points.length;

	for (int n = 0; n < nPoints; n++) {
	    sumX += points[n][0];
	    sumY += points[n][1];
	}

	centroid[0] = sumX / nPoints;
	centroid[1] = sumY / nPoints;

	return centroid;
    }

    public void printMatrix(Matrix matrix, String title) {
	IJ.log(title);
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
	return;
    }

    /**
     * Generate coordinates of a circle
     * 
     * @param x
     *            x coordinate of centre
     * @param y
     *            y coordinate of centre
     * @param r
     *            radius of circle
     * @param n
     *            Number of coordinates
     * @param noise
     *            Add noise of intensity 'noise'
     * 
     * @return
     */
    public double[][] getTestCircle(double x, double y, double r, int n,
	    double noise) {
	double[][] testCircle = new double[n][2];

	for (int i = 0; i < n; i++) {
	    double theta = i * 2 * Math.PI / n;
	    testCircle[i][0] = r * (1 + noise * (Math.random() - 0.5))
		    * Math.sin(theta) + x;
	    testCircle[i][1] = r * (1 + noise * (Math.random() - 0.5))
		    * Math.cos(theta) + y;
	    IJ.log("testCircle[n] is (" + testCircle[i][0] + ", "
		    + testCircle[i][1] + ")");
	}

	return testCircle;
    }

    /**
     * Return the column in Matrix D that contains the nth smallest diagonal
     * value
     * 
     * @param D
     *            the matrix to search
     * @param n
     *            the order of the diagonal value to find (1 = smallest, 2 =
     *            second smallest, etc.)
     * @return column index of the nth smallest diagonal in D
     */
    private int getNthSmallestCol(Matrix D, int n) {
	double[] diagD = new double[D.getColumnDimension()];
	double[] dDtest = new double[D.getColumnDimension()];
	for (int i = 0; i < D.getColumnDimension(); i++) {
	    diagD[i] = D.get(i, i);
	    dDtest[i] = D.get(i, i);
	}
	Arrays.sort(diagD);
	double sSvalue = diagD[n - 1];
	int col = Arrays.binarySearch(dDtest, sSvalue);

	if (diagD[0] > 0) {
	    IJ.error("Error: the smallest e-value is positive...");
	}
	if (diagD[1] < 0) {
	    IJ.error("Error: the second smallest e-value is negative...");
	}
	return col;
    }

}
