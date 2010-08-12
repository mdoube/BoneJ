package org.doube.geometry;

public class Centroid {

	/**
	 * Find the centroid of an array in double[n][i] format,
	 * where n = number of points and i = number of dimensions
	 * 
	 * @param points
	 * @return array containing centroid in i dimensions
	 */
	public static double[] getCentroid(double[][] points) {
		final int nDimensions = points[0].length;

		switch (nDimensions) {
		case 1:
			return getCentroid1D(points);
		case 2:
			return getCentroid2D(points);
		case 3:
			return getCentroid3D(points);
		default:
			return getCentroidND(points);
		}
	}

	/**
	 * Find the centroid of a set of points in double[n][1] format
	 * 
	 * @param points
	 * @return
	 */
	private static double[] getCentroid1D(double[][] points) {
		double[] centroid = new double[1];
		double sumX = 0;
		final int nPoints = points.length;

		for (int n = 0; n < nPoints; n++) {
			sumX += points[n][0];
		}

		centroid[0] = sumX / nPoints;

		return centroid;
	}

	/**
	 * Find the centroid of a set of points in double[n][2] format
	 * 
	 * @param points
	 * @return
	 */
	private static double[] getCentroid2D(double[][] points) {
		double[] centroid = new double[2];
		double sumX = 0;
		double sumY = 0;
		final int nPoints = points.length;

		for (int n = 0; n < nPoints; n++) {
			sumX += points[n][0];
			sumY += points[n][1];
		}

		centroid[0] = sumX / nPoints;
		centroid[1] = sumY / nPoints;

		return centroid;
	}

	/**
	 * Find the centroid of a set of points in double[n][3] format
	 * 
	 * @param points
	 * @return
	 */
	private static double[] getCentroid3D(double[][] points) {
		double[] centroid = new double[3];
		double sumX = 0;
		double sumY = 0;
		double sumZ = 0;
		final int nPoints = points.length;

		for (int n = 0; n < nPoints; n++) {
			sumX += points[n][0];
			sumY += points[n][1];
			sumZ += points[n][2];
		}

		centroid[0] = sumX / nPoints;
		centroid[1] = sumY / nPoints;
		centroid[2] = sumZ / nPoints;

		return centroid;
	}

	/**
	 * Find the centroid of a set of points in double[n][i] format
	 * 
	 * @param points
	 * @return
	 */
	private static double[] getCentroidND(double[][] points) {
		final int nPoints = points.length;
		final int nDimensions = points[0].length;
		double[] centroid = new double[nDimensions];
		double[] sums = new double[nDimensions];

		for (int n = 0; n < nPoints; n++) {
			for (int i = 0; i < nDimensions; i++) {
				sums[i] += points[n][i];
			}
		}

		for (int i = 0; i < nDimensions; i++) {
			centroid[i] = sums[i] / nPoints;
		}

		return centroid;
	}
	
	/**
	 * Return the centroid of a 1D array, which is its mean value
	 * 
	 * @param points
	 * @return the mean value of the points
	 */
	public static double getCentroid(double[] points){
		final int nPoints = points.length;
		double sum = 0;
		for (int n = 0; n < nPoints; n++){
			sum += points[n];
		}
		return sum / nPoints;
	}
	
	/**
	 * Return the centroid of a 1D int[] array, which is its mean value (rounded down)
	 * 
	 * @param points
	 * @return the mean value of the points
	 */
	public static int getCentroid(int[] points){
		final int nPoints = points.length;
		double sum = 0;
		for (int n = 0; n < nPoints; n++){
			sum += points[n];
		}
		return (int) Math.round(sum / nPoints);
	}
}
