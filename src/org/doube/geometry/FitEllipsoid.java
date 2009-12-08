package org.doube.geometry;

import ij.IJ;

import org.doube.jama.EigenvalueDecomposition;
import org.doube.jama.Matrix;

public class FitEllipsoid {
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
	public double[][] fitEllipsoid(double[][] coOrdinates) {

		IJ.showStatus("Fitting ellipsoid");
		double[][] Darray = new double[coOrdinates.length][10];
		for (int n = 0; n < coOrdinates.length; n++) {
			// populate each column of D with ten-element Xi
			final double xi = coOrdinates[n][0];
			final double yi = coOrdinates[n][1];
			final double zi = coOrdinates[n][2];
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
		final int evl = EigenValueMatrix.length;
		for (int p = 0; p < evl; p++) {
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
}
