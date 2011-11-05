package org.doube.geometry;

import static org.junit.Assert.*;

import org.junit.Test;

public class FitPlaneTest {

	@Test
	public void testFitPlane() {
		final double sqrt3 = 1/Math.sqrt(3);
		final double[][] points = planePoints(1,2,3, sqrt3, sqrt3, sqrt3, 1, 10000000);
		double[] plane = {sqrt3, sqrt3, sqrt3, 1, 2, 3};
		assertArrayEquals(plane, FitPlane.fitPlane(points), 1e-3);
	}
	
	/**
	 * Generate n random points on a plane given centroid (x1, y1, z1) and normal
	 * unit vector (a, b, c)
	 * 
	 * @param x1
	 * @param y1
	 * @param z1
	 * @param a
	 * @param b
	 * @param c
	 * @param n
	 * @return
	 */
	private double[][] planePoints(double x1, double y1, double z1, double a,
			double b, double c, double range, int n) {
		final double d = -(a * x1 + b * y1 + c * z1);
		double[][] points = new double[n][3];
		for (int i = 0; i < n; i++) {
			final double x = (Math.random() - 0.5) * range * 2;
			final double y = (Math.random() - 0.5) * range * 2;
			// ax + by + cz + d = 0, solve for z
			final double z = -(a * x + b * y + d) / c;
			points[i][0] = x;
			points[i][1] = y;
			points[i][2] = z;
		}
		return points;
	}
	
}
