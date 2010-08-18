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
		double[] p0 = new double[3];
		double[] p1 = new double[3];
		double[] p2 = new double[3];
		p0[0] = point0.x;
		p0[1] = point0.y;
		p0[2] = point0.z;
		p1[0] = point1.x;
		p1[1] = point1.y;
		p1[2] = point1.z;
		p2[0] = point2.x;
		p2[1] = point2.y;
		p2[2] = point2.z;
		double[] cV = crossProduct(p0, p1, p2);
		Point3f crossVector = new Point3f();
		crossVector.x = (float) cV[0];
		crossVector.y = (float) cV[1];
		crossVector.z = (float) cV[2];
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

		double[] crossVector = new double[3];
		crossVector[0] = y1 * z2 - z1 * y2;
		crossVector[1] = z1 * x2 - x1 * z2;
		crossVector[2] = x1 * y2 - y1 * x2;

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
		double[] c = new double[3];
		c[0] = a[1] * b[2] - a[2] * b[1];
		c[1] = a[2] * b[0] - a[0] * b[2];
		c[2] = a[0] * b[1] - a[1] * b[0];
		return c;
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
