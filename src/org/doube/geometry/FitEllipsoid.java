package org.doube.geometry;

/**
 *  FitEllipsoid Copyright 2009 2010 Michael Doube
 * 
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

import ij.IJ;

import org.doube.jama.EigenvalueDecomposition;
import org.doube.jama.Matrix;

/**
 * Ellipsoid fitting methods. Both rely on eigenvalue decomposition, which fails
 * if the input matrix is singular. It is worth enclosing calls to yuryPetrov
 * and liGriffiths in a try{FitEllipsoid.yuryPetrov} catch(RuntimeException
 * re){} to gracefully handle cases where the methods cannot find a best-fit
 * ellipsoid.
 * 
 * @author Michael Doube
 * 
 */
public class FitEllipsoid {
	/**
	 * Calculate the best-fit ellipsoid by least squares. Currently broken;
	 * doesn't handle large numbers of points and often returns a singular
	 * matrix for EVD.
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
	 * @deprecated until debugged
	 */
	public static double[] liGriffiths(double[][] coOrdinates, int maxPoints) {
		final int nPoints = coOrdinates.length;
		IJ.showStatus("Fitting ellipsoid");
		double[][] coOrd = new double[maxPoints][3];
		if (nPoints > maxPoints) {
			// randomly subsample
			for (int n = 0; n < maxPoints; n++) {
				final int point = (int) Math.floor(Math.random() * nPoints);
				coOrd[n][0] = coOrdinates[point][0];
			}
		}
		double[][] Darray = new double[maxPoints][10];
		for (int n = 0; n < maxPoints; n++) {
			// populate each column of D with ten-element Xi
			final double xi = coOrd[n][0];
			final double yi = coOrd[n][1];
			final double zi = coOrd[n][2];
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
		Matrix C1inv = C1.inverse();
		Matrix S12T = S12.transpose();
		Matrix S22inv = S22.inverse();
		Matrix Eq15 = C1inv.times(S11.minus(S12.times(S22inv.times(S12T))));
		// However, the matrix (S11 - S12*S22^-1*S12^T ) can be singular in some
		// cases. In this situations, the corresponding u1 can be replaced
		// with the eigenvector associated with the largest eigenvalue.

		EigenvalueDecomposition E = new EigenvalueDecomposition(Eq15);
		Matrix eigenValues = E.getD();
		eigenValues.printToIJLog("eigenValues");
		Matrix eigenVectors = E.getV();
		eigenVectors.printToIJLog("eigenVectors");
		Matrix v1 = new Matrix(6, 1);
		double[][] EigenValueMatrix = eigenValues.getArray();
		double posVal = 1 - Double.MAX_VALUE;
		final int evl = EigenValueMatrix.length;
		for (int p = 0; p < evl; p++) {
			if (EigenValueMatrix[p][p] > posVal) {
				posVal = EigenValueMatrix[p][p];
				v1 = eigenVectors.getMatrix(0, 5, p, p);
			}
		}
		v1.printToIJLog("v1");
		Matrix v2 = S22.inverse().times(S12.transpose()).times(v1).times(-1);
		v2.printToIJLog("v2");
		Matrix v = new Matrix(10, 1);
		int[] c = { 0 };
		v.setMatrix(0, 5, c, v1);
		v.setMatrix(6, 9, c, v2);
		double[] ellipsoid = v.getColumnPackedCopy();
		IJ.log("Ellipsoid equation: " + ellipsoid[0] + " x² + " + ellipsoid[1]
				+ " y² + " + ellipsoid[2] + " z² + " + ellipsoid[3] + " 2yz + "
				+ ellipsoid[4] + " 2xz + " + ellipsoid[5] + " 2xy + "
				+ ellipsoid[6] + " 2x + " + ellipsoid[7] + " 2y + "
				+ ellipsoid[8] + " 2z + " + ellipsoid[9] + " = 0");
		return ellipsoid;
	}

	/**
	 * <p>
	 * Ellipsoid fitting method by Yury Petrov.<br />
	 * Fits an ellipsoid in the form <i>Ax</i><sup>2</sup> +
	 * <i>By</i><sup>2</sup> + <i>Cz</i><sup>2</sup> + 2<i>Dxy</i> + 2<i>Exz</i>
	 * + 2<i>Fyz</i> + 2<i>Gx</i> + 2<i>Hy</i> + 2<i>Iz</i> = 1 <br />
	 * To an n * 3 array of coordinates.
	 * </p>
	 * 
	 * @see <p>
	 *      <a href=
	 *      "http://www.mathworks.com/matlabcentral/fileexchange/24693-ellipsoid-fit"
	 *      >MATLAB script</a>
	 *      </p>
	 * 
	 * @param coOrdinates
	 *            array[n][3] where n > 8
	 * @return Object[] array containing the centre, radii, eigenvectors of the
	 *         axes, the 9 variables of the ellipsoid equation and the EVD
	 * @throws IllegalArgumentException
	 *             if number of coordinates is less than 9
	 */
	public static Object[] yuryPetrov(double[][] coOrdinates) {

		final int nPoints = coOrdinates.length;
		if (nPoints < 9) {
			throw new IllegalArgumentException(
					"Too few points; need at least 9 to calculate a unique ellipsoid");
		}

		double[][] d = new double[nPoints][9];
		for (int i = 0; i < nPoints; i++) {
			final double x = coOrdinates[i][0];
			final double y = coOrdinates[i][1];
			final double z = coOrdinates[i][2];
			d[i][0] = x * x;
			d[i][1] = y * y;
			d[i][2] = z * z;
			d[i][3] = 2 * x * y;
			d[i][4] = 2 * x * z;
			d[i][5] = 2 * y * z;
			d[i][6] = 2 * x;
			d[i][7] = 2 * y;
			d[i][8] = 2 * z;
		}

		Matrix D = new Matrix(d);
		Matrix ones = Matrix.ones(nPoints, 1);
		Matrix V = ((D.transpose().times(D)).inverse()).times(D.transpose()
				.times(ones));
		double[] v = V.getColumnPackedCopy();

		double[][] a = { { v[0], v[3], v[4], v[6] },
				{ v[3], v[1], v[5], v[7] }, { v[4], v[5], v[2], v[8] },
				{ v[6], v[7], v[8], -1 }, };
		Matrix A = new Matrix(a);
		Matrix C = (A.getMatrix(0, 2, 0, 2).times(-1).inverse()).times(V
				.getMatrix(6, 8, 0, 0));
		Matrix T = Matrix.eye(4);
		T.setMatrix(3, 3, 0, 2, C.transpose());
		Matrix R = T.times(A.times(T.transpose()));
		double r33 = R.get(3, 3);
		Matrix R02 = R.getMatrix(0, 2, 0, 2);
		EigenvalueDecomposition E = new EigenvalueDecomposition(R02.times(-1
				/ r33));
		Matrix eVal = E.getD();
		Matrix eVec = E.getV();
		Matrix diagonal = eVal.diag();
		final int nEvals = diagonal.getRowDimension();
		double[] radii = new double[nEvals];
		for (int i = 0; i < nEvals; i++) {
			radii[i] = Math.sqrt(1 / diagonal.get(i, 0));
		}
		double[] centre = C.getColumnPackedCopy();
		double[][] eigenVectors = eVec.getArrayCopy();
		double[] equation = v;
		Object[] ellipsoid = { centre, radii, eigenVectors, equation, E };
		return ellipsoid;
	}

	/**
	 * Return points on an ellipsoid with optional noise. Point density is not
	 * uniform, becoming more dense at the poles.
	 * 
	 * @param a
	 *            First axis length
	 * @param b
	 *            Second axis length
	 * @param c
	 *            Third axis length
	 * @param angle
	 *            angle of axis (rad)
	 * @param xCentre
	 *            x coordinate of centre
	 * @param yCentre
	 *            y coordinate of centre
	 * @param zCentre
	 *            z coordinate of centre
	 * @param noise
	 *            Intensity of noise to add to the points
	 * @param nPoints
	 *            number of points to generate
	 * @param random
	 *            if true, use a random grid to generate points, else use a
	 *            regular grid
	 * @return array of (x,y,z) coordinates
	 */
	public static double[][] testEllipsoid(double a, double b, double c,
			double angle, double xCentre, double yCentre, double zCentre,
			double noise, int nPoints, boolean random) {

		final int n = (int) Math.floor(-3 / 4 + Math.sqrt(1 + 8 * nPoints) / 4);
		final int h = 2 * n + 2;
		final int w = n + 1;
		double[][] s = new double[h][w];
		double[][] t = new double[h][w];
		double value = -Math.PI / 2;
		// Random points
		if (random) {
			for (int j = 0; j < w; j++) {
				for (int i = 0; i < h; i++) {
					s[i][j] = value + Math.random() * 2 * Math.PI;
				}
			}
			for (int i = 0; i < h; i++) {
				for (int j = 0; j < w; j++) {
					t[i][j] = value + Math.random() * 2 * Math.PI;
				}
			}
			// Regular points
		} else {
			final double increment = Math.PI / (n - 1);

			for (int j = 0; j < w; j++) {
				for (int i = 0; i < h; i++) {
					s[i][j] = value;
				}
				value += increment;
			}
			value = -Math.PI / 2;
			for (int i = 0; i < h; i++) {
				for (int j = 0; j < w; j++) {
					t[i][j] = value;
				}
				value += increment;
			}

		}
		double[][] x = new double[h][w];
		double[][] y = new double[h][w];
		double[][] z = new double[h][w];
		for (int i = 0; i < h; i++) {
			for (int j = 0; j < w; j++) {
				x[i][j] = a * Math.cos(s[i][j]) * Math.cos(t[i][j]);
				y[i][j] = b * Math.cos(s[i][j]) * Math.sin(t[i][j]);
				z[i][j] = c * Math.sin(s[i][j]);
			}
		}
		double[][] xt = new double[h][w];
		double[][] yt = new double[h][w];
		for (int i = 0; i < h; i++) {
			for (int j = 0; j < w; j++) {
				xt[i][j] = x[i][j] * Math.cos(angle) - y[i][j]
						* Math.sin(angle);
				yt[i][j] = x[i][j] * Math.sin(angle) + y[i][j]
						* Math.cos(angle);
			}
		}
		for (int i = 0; i < h; i++) {
			for (int j = 0; j < w; j++) {
				x[i][j] = xt[i][j] + xCentre;
				y[i][j] = yt[i][j] + yCentre;
				z[i][j] = z[i][j] + zCentre;
			}
		}
		for (int i = 0; i < h; i++) {
			for (int j = 0; j < w; j++) {
				x[i][j] = x[i][j] + Math.random() * noise;
				y[i][j] = y[i][j] + Math.random() * noise;
				z[i][j] = z[i][j] + Math.random() * noise;
			}
		}
		double[][] ellipsoidPoints = new double[w * h][3];
		int p = 0;
		for (int i = 0; i < h; i++) {
			for (int j = 0; j < w; j++) {
				ellipsoidPoints[p][0] = x[i][j];
				ellipsoidPoints[p][1] = y[i][j];
				ellipsoidPoints[p][2] = z[i][j];
				p++;
			}
		}
		return ellipsoidPoints;
	}

	/**
	 * Return normal unit vectors at points on an ellipsoid given the radii of
	 * an ellipsoid and points on the ellipsoid. Assumes an ellipsoid centred on
	 * (0,0,0) and with no rotation
	 * 
	 * @param a
	 *            First axis length
	 * @param b
	 *            Second axis length
	 * @param c
	 *            Third axis length
	 * @param points
	 *            points on the ellipsoid, e.g. result of testElipsoid()
	 * @return array of (x,y,z) unit vectors
	 */
	public static double[][] testNormals(double a, double b, double c,
			double[][] points) {
		final int nPoints = points.length;
		double[][] ellipsoidNormals = new double[nPoints][3];
		final double p = 2 / (a * a);
		final double q = 2 / (b * b);
		final double r = 2 / (c * c);
		for (int i = 0; i < nPoints; i++) {
			ellipsoidNormals[i][0] = p * points[i][0];
			ellipsoidNormals[i][1] = q * points[i][1];
			ellipsoidNormals[i][2] = r * points[i][2];
		}

		for (int i = 0; i < nPoints; i++) {
			final double x = ellipsoidNormals[i][0];
			final double y = ellipsoidNormals[i][1];
			final double z = ellipsoidNormals[i][2];
			final double length = Trig.distance3D(x, y, z);
			ellipsoidNormals[i][0] /= length;
			ellipsoidNormals[i][1] /= length;
			ellipsoidNormals[i][2] /= length;
		}

		return ellipsoidNormals;
	}
}
