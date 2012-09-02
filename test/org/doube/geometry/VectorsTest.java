package org.doube.geometry;

import static org.junit.Assert.*;

import javax.vecmath.Point3f;

import org.junit.Test;

public class VectorsTest {
	Point3f p0 = new Point3f(1.0f, 2.0f, 3.0f);
	Point3f p1 = new Point3f(6.0f, 5.0f, 4.0f);
	Point3f p2 = new Point3f(9.0f, 8.0f, 7.0f);

	double[] point0 = { 1, 2, 3 };
	double[] point1 = { 6, 5, 4 };
	double[] point2 = { 9, 8, 7 };

	double[][] pa0 = { { 5 }, { 3 }, { 1 } };
	double[][] pa1 = { { 8 }, { 6 }, { 4 } };

	double[] a = { 5, 3, 1 };
	double[] b = { 8, 6, 4 };

	@Test
	public void testCrossProductPoint3fPoint3fPoint3f() {
		Point3f result = Vectors.crossProduct(p0, p1, p2);
		Point3f expected = new Point3f(6.0f, -12.0f, 6.0f);
		assertEquals(expected, result);
	}

	@Test
	public void testCrossProductDoubleArrayDoubleArrayDoubleArray() {
		double[] result = Vectors.crossProduct(point0, point1, point2);
		double[] expected = { 6, -12, 6 };
		assertArrayEquals(expected, result, 1e-12);
	}

	@Test
	public void testCrossProductDoubleArrayArrayDoubleArrayArray() {
		double[][] result = Vectors.crossProduct(pa0, pa1);
		double[][] expected = { { 6 }, { -12 }, { 6 } };
		assertEquals(expected[0][0], result[0][0], 1e-12);
		assertEquals(expected[1][0], result[1][0], 1e-12);
		assertEquals(expected[2][0], result[2][0], 1e-12);
	}

	@Test
	public void testCrossProductDoubleArrayDoubleArray() {
		double[] result = Vectors.crossProduct(a, b);
		double[] expected = { 6, -12, 6 };
		assertArrayEquals(expected, result, 1e-12);
	}

	@Test
	public void testRandomVectors() {
		final int n = 1000;
		double[][] v = Vectors.randomVectors(n);
		// check that vectors are unit vectors
		for (int i = 0; i < n; i++) {
			final double x = v[i][0];
			final double y = v[i][1];
			final double z = v[i][2];
			final double length = Math.sqrt(x * x + y * y + z * z);
			assertEquals(1, length, 1e-9);
		}
	}
	
	@Test
	public void testRegularVectors() {
		final int n = 1000;
		double[][] v = Vectors.regularVectors(n);
		// check that vectors are unit vectors
		for (int i = 0; i < n; i++) {
			final double x = v[i][0];
			final double y = v[i][1];
			final double z = v[i][2];
			final double length = Math.sqrt(x * x + y * y + z * z);
			assertEquals(1, length, 1e-9);
		}
	}
}
