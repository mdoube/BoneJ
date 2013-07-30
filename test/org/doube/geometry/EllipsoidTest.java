package org.doube.geometry;

import static org.junit.Assert.*;

import org.junit.Test;

public class EllipsoidTest {

	/** Ellipsoid of radii = 1 and centred on origin */
	Ellipsoid unitSphere = new Ellipsoid(1, 1, 1, 0, 0, 0, new double[][] {
			{ 0, 0, 1 }, { 0, 1, 0 }, { 1, 0, 0 } });

	/** Ellipsoid with radii of 1, 2, and 3, centred on (1,2,3) */
	Ellipsoid oneTwoThree = FitEllipsoid.fitTo(FitEllipsoid.testEllipsoid(1, 2,
			3, 0, 1, 2, 3, 0, 10000, true));

	@Test
	public void testGetVolume() {
		assertEquals(Math.PI * 4 / 3, unitSphere.getVolume(), 1E-9);
	}

	@Test
	public void testGetRadii() {
		assertArrayEquals(new double[] { 1, 1, 1 }, unitSphere.getRadii(), 1E-9);

		assertArrayEquals(new double[] { 3, 2, 1 }, oneTwoThree.getRadii(),
				1E-9);
	}

	@Test
	public void testGetEquation() {

		assertArrayEquals(new double[] { 1, 1, 1, 0, 0, 0, 0, 0, 0 },
				FitEllipsoid.fitTo(unitSphere.getSurfacePoints(10000))
						.getEquation(), 1E-9);

		assertArrayEquals(new double[] { -0.5, -0.125, -1.0 / 18.0, 0, 0, 0,
				0.5, 0.25, 1.0 / 6.0 }, oneTwoThree.getEquation(), 1E-9);
	}

	@Test
	public void testContains() {
		double[][] points = unitSphere.getSurfacePoints(1000);
		Ellipsoid e = FitEllipsoid.fitTo(points);
		double[][] vectors = Vectors.randomVectors(1000);

		// inside
		for (double[] v : vectors) {
			double rand = Math.random();
			double x = v[0] * rand;
			double y = v[1] * rand;
			double z = v[2] * rand;
			assertTrue("Testing (" + x + ", " + y + ", " + z + ")",
					e.contains(x, y, z));
		}

		// outside
		for (double[] v : vectors) {
			double rand = Math.random();
			double x = v[0] / rand;
			double y = v[1] / rand;
			double z = v[2] / rand;
			assertTrue("Testing (" + x + ", " + y + ", " + z + ")",
					!e.contains(x, y, z));
		}

		assertTrue(
				"Testing centroid (1, 2, 3) with solution = "
						+ oneTwoThree.solve(1, 2, 3),
				oneTwoThree.contains(1, 2, 3));

		assertTrue("Testing wildly outside value with solution = "
				+ oneTwoThree.solve(23, 213, 132),
				!oneTwoThree.contains(23, 213, 132));

	}

	@Test
	public void testSolve() {
		double[][] points = unitSphere.getSurfacePoints(1000);
		Ellipsoid e = FitEllipsoid.fitTo(points);
		for (double[] p : points) {
			assertEquals(1.0, e.solve(p[0], p[1], p[2]), 1E-9);
		}

		points = oneTwoThree.getSurfacePoints(1000);
		e = FitEllipsoid.fitTo(points);
		for (double[] p : points) {
			assertEquals(1.0, e.solve(p[0], p[1], p[2]), 1E-9);
		}

	}

	@Test
	public void testGetMajorRadius() {
		assertEquals(1, unitSphere.getMajorRadius(), 1E-9);

		assertEquals(3, oneTwoThree.getMajorRadius(), 1E-9);
	}

	@Test
	public void testGetCentre() {
		assertArrayEquals(new double[] { 0, 0, 0 }, unitSphere.getCentre(),
				1E-9);

		assertArrayEquals(new double[] { 1, 2, 3 }, oneTwoThree.getCentre(),
				1E-9);
	}

	@Test
	public void testGetSurfacePoints() {
		double[][] points = unitSphere.getSurfacePoints(10000);
		for (double[] p : points) {
			assertEquals(1, Trig.distance3D(p), 1E-9);
		}

		points = oneTwoThree.getSurfacePoints(10000);
		for (double[] p : points) {
			assertEquals(1, oneTwoThree.solve(p[0], p[1], p[2]), 1E-9);
		}

	}

	@Test
	public void testDilate() {
		unitSphere.dilate(1);
		double[][] points = unitSphere.getSurfacePoints(10000);
		for (double[] p : points) {
			assertEquals(2, Trig.distance3D(p), 1E-9);
		}
	}

	@Test
	public void testContract() {
		unitSphere.contract(1.5);
		double[][] points = unitSphere.getSurfacePoints(10000);
		for (double[] p : points) {
			assertEquals(0.5, Trig.distance3D(p), 1E-9);
		}
	}

}
