package org.doube.geometry;

import org.doube.jama.EigenvalueDecomposition;
import org.doube.jama.Matrix;

public class EigenCalculator {

	/**
	 * Calculate the eigenvectors and eigenvalues of a set of points by the
	 * covariance method and eigendecomposition.
	 * 
	 * @param coOrdinates
	 *            n x 3 array centred on (0,0,0)
	 * @return EigenvalueDecomposition containing eigenvectors and eigenvalues
	 * 
	 */
	public static EigenvalueDecomposition principalComponents(
			double[][] coOrdinates) {
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
	}
}
