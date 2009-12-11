package org.doube.geometry;

import org.doube.jama.Matrix;
import org.doube.jama.SingularValueDecomposition;

public class FitPlane {
	/**
	 * Find the best fit plane using the Orthogonal Distance Regression.
	 * 
	 * @param points
	 *            List of 3D coordinates in double[n][3] format
	 * @return double[6] containing plane vector in first 3 elements and
	 *         centroid in last 3 elements
	 */
	public static double[] fitPlane(double[][] points) {
		final int nPoints = points.length;
		final double[] centroid = Centroid.getCentroid(points);
		double sDxDx = 0;
		double sDyDy = 0;
		double sDzDz = 0;
		double sDxDy = 0;
		double sDxDz = 0;
		double sDyDz = 0;
		for (int p = 0; p < nPoints; p++) {
			final double dx = points[p][0] - centroid[0];
			final double dy = points[p][1] - centroid[1];
			final double dz = points[p][2] - centroid[2];
			sDxDx += dx * dx;
			sDyDy += dy * dy;
			sDzDz += dz * dz;
			sDxDy += dx * dy;
			sDxDz += dx * dz;
			sDyDz += dy * dz;
		}
		double[][] C = new double[3][3];
		C[0][0] = sDxDx;
		C[1][1] = sDyDy;
		C[2][2] = sDzDz;
		C[0][1] = sDxDy;
		C[0][2] = sDxDz;
		C[1][0] = sDxDy;
		C[1][2] = sDyDz;
		C[2][0] = sDxDz;
		C[2][1] = sDyDz;
		final double invCount = 1 / nPoints;
		Matrix covarianceMatrix = new Matrix(C).times(invCount);
		covarianceMatrix.printToIJLog("Covariance matrix");
		SingularValueDecomposition S = new SingularValueDecomposition(
				covarianceMatrix);
		Matrix leftVectors = S.getU();
		double[] plane = new double[6];
		// first 3 elements are the vector perpendicular to the plane
		plane[0] = leftVectors.get(0, 0);
		plane[1] = leftVectors.get(1, 0);
		plane[2] = leftVectors.get(2, 0);
		// last 3 elements are the centroid
		plane[3] = centroid[0];
		plane[4] = centroid[1];
		plane[5] = centroid[2];
		return plane;
	}
}
