package org.doube.geometry;

import static org.junit.Assert.*;

import org.doube.bonej.EllipsoidFitter;
import org.junit.Test;

public class EllipsoidTest {

	/** Ellipsoid of radii = 1 and centred on origin */
	Ellipsoid unitSphere = new Ellipsoid(1, 1, 1, 0, 0, 0, new double[][] {
			{ 0, 0, 1 }, { 0, 1, 0 }, { 1, 0, 0 } });

	/** Ellipsoid of radii = 1 and centred on (17, 31, 71) */
	Ellipsoid unitSphereTrans = FitEllipsoid.fitTo(FitEllipsoid.testEllipsoid(
			1, 1, 1, 0, 17, 31, 71, 0, 1000, true));

	/** Ellipsoid with radii of 1, 2, and 3, centred on (1,2,3) */
	Ellipsoid oneTwoThree = FitEllipsoid.fitTo(FitEllipsoid.testEllipsoid(1, 2,
			3, 0, 1, 2, 3, 0, 10000, true));

	/** Ellipsoid with radii of 3, 5, 17, centred on (0, 0, 0) */
	Ellipsoid threeFiveSeventeen = FitEllipsoid.fitTo(FitEllipsoid
			.testEllipsoid(3, 5, 17, 0, 0, 0, 0, 0, 10000, true));

	@Test
	public void testGetVolume() {
		assertEquals(Math.PI * 4 / 3, unitSphere.getVolume(), 1E-9);
		assertEquals(1 * 2 * 3 * Math.PI * 4 / 3, oneTwoThree.getVolume(), 1E-9);
		assertEquals(3 * 5 * 17 * Math.PI * 4 / 3,
				threeFiveSeventeen.getVolume(), 1E-9);
	}

	@Test
	public void testGetRadii() {
		assertArrayEquals(new double[] { 1, 1, 1 }, unitSphere.getRadii(), 1E-9);

		assertArrayEquals(new double[] { 3, 2, 1 }, oneTwoThree.getRadii(),
				1E-9);

		assertArrayEquals(new double[] { 17, 5, 3 },
				threeFiveSeventeen.getRadii(), 1E-9);
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
	public void testIntercepts() {
		Ellipsoid e = FitEllipsoid.fitTo(unitSphere.getSurfacePoints(1000));
		assertArrayEquals(new double[] { 1, -1, 1, -1, 1, -1 }, e.intercepts(),
				1E-9);
	}

	@Test
	public void testContains() {
		double[][] points = unitSphere.getSurfacePoints(1000);

		// centroid
		assertTrue(unitSphere.contains(0, 0, 0));

		// inside
		for (double[] v : points) {
			double rand = Math.random();
			double x = v[0] * rand;
			double y = v[1] * rand;
			double z = v[2] * rand;
			assertTrue("Testing (" + x + ", " + y + ", " + z + ")",
					unitSphere.contains(x, y, z));
		}

		// outside
		for (double[] v : points) {
			double rand = Math.random();
			double x = v[0] / rand;
			double y = v[1] / rand;
			double z = v[2] / rand;
			assertTrue("Testing (" + x + ", " + y + ", " + z + ")",
					!unitSphere.contains(x, y, z));
		}

		assertTrue("Testing centroid (1, 2, 3)", oneTwoThree.contains(1, 2, 3));

		// random vectors greater than major radius
		double[][] vectors = Vectors.randomVectors(1000);
		for (double[] v : vectors) {
			double x = v[0] * 18;
			double y = v[1] * 18;
			double z = v[2] * 18;
			assertTrue(!threeFiveSeventeen.contains(x, y, z));
		}

		// random vectors smaller than the minor radius
		for (double[] v : vectors) {
			double x = v[0] * 2;
			double y = v[1] * 2;
			double z = v[2] * 2;
			assertTrue(threeFiveSeventeen.contains(x, y, z));
		}

		double[][] pointsU = unitSphereTrans.getSurfacePoints(1000);

		// random vectors translated
		for (double[] v : pointsU) {
			double rand = Math.random();
			double x = 17 + (v[0] - 17) * rand;
			double y = 31 + (v[1] - 31) * rand;
			double z = 71 + (v[2] - 71) * rand;
			assertTrue(unitSphereTrans.contains(x, y, z));
		}

		for (double[] v : pointsU) {
			double rand = Math.random();
			double x = 17 + (v[0] - 17) / rand;
			double y = 31 + (v[1] - 31) / rand;
			double z = 71 + (v[2] - 71) / rand;
			assertTrue(!unitSphereTrans.contains(x, y, z));
		}

		double[][] points2 = threeFiveSeventeen.getSurfacePoints(1000);

		// inside
		for (double[] p : points2) {
			// contract by random fraction
			final double rand = Math.random();
			double x = p[0] * rand;
			double y = p[1] * rand;
			double z = p[2] * rand;
			System.out.println("Testing (" + x + ", " + y + ", " + z + ")");
			assertTrue(threeFiveSeventeen.contains(x, y, z));
		}

		// outside
		for (double[] p : points2) {
			// dilate by random fraction
			final double rand = Math.random();
			double x = p[0] / rand;
			double y = p[1] / rand;
			double z = p[2] / rand;
			System.out.println("Testing (" + x + ", " + y + ", " + z + ")");
			assertTrue(!threeFiveSeventeen.contains(x, y, z));
		}
	}

	@Test
	public void testSolve() {
		double[][] points = unitSphere.getSurfacePoints(1000);
		Ellipsoid e = FitEllipsoid.fitTo(points);
		for (double[] p : points) {
			assertEquals(1.0, e.solve(p[0], p[1], p[2]), 1E-9);
		}

		points = oneTwoThree.getSurfacePoints(1000);
		for (double[] p : points) {
			assertEquals(1.0, oneTwoThree.solve(p[0], p[1], p[2]), 1E-9);
		}

		points = threeFiveSeventeen.getSurfacePoints(1000);
		for (double[] p : points) {
			assertEquals(1.0, threeFiveSeventeen.solve(p[0], p[1], p[2]), 1E-9);
		}

	}

	@Test
	public void testGetMajorRadius() {
		assertEquals(1, unitSphere.getMajorRadius(), 1E-9);
		assertEquals(3, oneTwoThree.getMajorRadius(), 1E-9);
		assertEquals(17, threeFiveSeventeen.getMajorRadius(), 1E-9);
	}

	@Test
	public void testGetCentre() {
		assertArrayEquals(new double[] { 0, 0, 0 }, unitSphere.getCentre(),
				1E-9);
		assertArrayEquals(new double[] { 1, 2, 3 }, oneTwoThree.getCentre(),
				1E-9);
		assertArrayEquals(new double[] { 0, 0, 0 },
				threeFiveSeventeen.getCentre(), 1E-9);
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

		points = threeFiveSeventeen.getSurfacePoints(10000);
		for (double[] p : points) {
			assertEquals(1, threeFiveSeventeen.solve(p[0], p[1], p[2]), 1E-9);
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
