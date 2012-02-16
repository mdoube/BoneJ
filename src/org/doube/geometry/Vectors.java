package org.doube.geometry;

import javax.vecmath.Point3f;

public class Vectors {
	/**
	 * Calculate the cross product of 3 Point3f's, which describe two vectors
	 * joined at the tails. Can be used to find the plane / surface normal of a
	 * triangle. Half of its magnitude is the area of the triangle.
	 * 
	 * @param point0
	 *            both vectors' tails
	 * @param point1
	 *            vector 1's head
	 * @param point2
	 *            vector 2's head
	 * @return cross product vector
	 */
	public static Point3f crossProduct(Point3f point0, Point3f point1,
			Point3f point2) {
		final double x1 = point1.x - point0.x;
		final double y1 = point1.y - point0.y;
		final double z1 = point1.z - point0.z;
		final double x2 = point2.x - point0.x;
		final double y2 = point2.y - point0.y;
		final double z2 = point2.z - point0.z;
		Point3f crossVector = new Point3f();
		crossVector.x = (float) (y1 * z2 - z1 * y2);
		crossVector.y = (float) (z1 * x2 - x1 * z2);
		crossVector.z = (float) (x1 * y2 - y1 * x2);
		return crossVector;
	}

	/**
	 * Calculate the cross product of 3 double[3]'s, which describe two vectors
	 * joined at the tails. Can be used to find the plane / surface normal of a
	 * triangle. Half of its magnitude is the area of the triangle.
	 * 
	 * @param point0
	 *            both vectors' tails
	 * @param point1
	 *            vector 1's head
	 * @param point2
	 *            vector 2's head
	 * 
	 * @return cross product vector
	 */
	public static double[] crossProduct(double[] point0, double[] point1,
			double[] point2) {
		final double x1 = point1[0] - point0[0];
		final double y1 = point1[1] - point0[1];
		final double z1 = point1[2] - point0[2];

		final double x2 = point2[0] - point0[0];
		final double y2 = point2[1] - point0[1];
		final double z2 = point2[2] - point0[2];

		return crossProduct(x1, y1, z1, x2, y2, z2);
	}

	/**
	 * Calculate the cross product of two vectors (x1, y1, z1) and (x2, y2, z2)
	 * 
	 * @param x1
	 * @param y1
	 * @param z1
	 * @param x2
	 * @param y2
	 * @param z2
	 * @return cross product in {x, y, z} format
	 */
	public static double[] crossProduct(double x1, double y1, double z1,
			double x2, double y2, double z2) {
		final double x = y1 * z2 - z1 * y2;
		final double y = z1 * x2 - x1 * z2;
		final double z = x1 * y2 - y1 * x2;
		double[] crossVector = { x, y, z };
		return crossVector;
	}

	/**
	 * Calculate the cross product of 2 column vectors, both in double[3][1]
	 * format
	 * 
	 * @param a
	 *            first vector
	 * @param b
	 *            second vector
	 * @return resulting vector in double[3][1] format
	 */
	public static double[][] crossProduct(double[][] a, double[][] b) {
		double[][] c = new double[3][1];
		c[0][0] = a[1][0] * b[2][0] - a[2][0] * b[1][0];
		c[1][0] = a[2][0] * b[0][0] - a[0][0] * b[2][0];
		c[2][0] = a[0][0] * b[1][0] - a[1][0] * b[0][0];
		return c;
	}

	/**
	 * Calculate the cross product of 2 vectors, both in double[3] format
	 * 
	 * @param a
	 *            first vector
	 * @param b
	 *            second vector
	 * @return resulting vector in double[3] format
	 */
	public static double[] crossProduct(double[] a, double[] b) {
		return crossProduct(a[0], a[1], a[2], b[0], b[1], b[2]);
	}

	/**
	 * Generate an array of randomly-oriented 3D unit vectors
	 * 
	 * @param nVectors
	 *            number of vectors to generate
	 * @return 2D array (nVectors x 3) containing unit vectors
	 */
	public static double[][] randomVectors(int nVectors) {
		double[][] randomVectors = new double[nVectors][3];
		for (int n = 0; n < nVectors; n++) {
			final double z = 2 * Math.random() - 1;
			double rho = Math.sqrt(1 - z * z);
			double phi = Math.PI * (2 * Math.random() - 1);
			randomVectors[n][0] = rho * Math.cos(phi);
			randomVectors[n][1] = rho * Math.sin(phi);
			randomVectors[n][2] = z;
		}
		return randomVectors;
	}

}
