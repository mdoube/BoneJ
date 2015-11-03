package org.doube.geometry;

import org.doube.jama.EigenvalueDecomposition;
import org.doube.jama.Matrix;

public class EigenCalculator {

	/**
	 * Calculate the eigenvectors and eigenvalues of a set of points by the
	 * covariance method and eigendecomposition.
	 *
	 * @param coOrdinates
	 *            n x 3 array
	 * @return EigenvalueDecomposition containing eigenvectors and eigenvalues
	 *
	 */
	public static EigenvalueDecomposition principalComponents(final double[][] coOrdinates) {
		final int nPoints = coOrdinates.length;
		double sumX = 0, sumY = 0, sumZ = 0;
		// calculate the centroid of the points
		for (int n = 0; n < nPoints; n++) {
			sumX += coOrdinates[n][0];
			sumY += coOrdinates[n][1];
			sumZ += coOrdinates[n][2];
		}
		// centroid is the mean (x, y, z) position
		final double centX = sumX / nPoints;
		final double centY = sumY / nPoints;
		final double centZ = sumZ / nPoints;

		// construct the covariance matrix
		final double[][] C = new double[3][3];
		double count = 0;
		for (int n = 0; n < nPoints; n++) {
			// translate so that centroid is at (0,0,0)
			final double x = coOrdinates[n][0] - centX;
			final double y = coOrdinates[n][1] - centY;
			final double z = coOrdinates[n][2] - centZ;
			C[0][0] += x * x;
			C[1][1] += y * y;
			C[2][2] += z * z;
			final double xy = x * y;
			final double xz = x * z;
			final double yz = y * z;
			C[0][1] += xy;
			C[0][2] += xz;
			C[1][0] += xy;
			C[1][2] += yz;
			C[2][0] += xz;
			C[2][1] += yz;
			count += 1;
		}
		final double invCount = 1 / count;
		final Matrix covarianceMatrix = new Matrix(C).times(invCount);
		final EigenvalueDecomposition E = new EigenvalueDecomposition(covarianceMatrix);
		return E;
	}
}
