package org.doube.geometry;

import static org.junit.Assert.*;

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

	/** Ellipsoid with radii of 17, 5, 3, centred on (0, 0, 0) */
	Ellipsoid seventeenFiveThree = FitEllipsoid.fitTo(FitEllipsoid
			.testEllipsoid(17, 5, 3, 0, 0, 0, 0, 0, 10000, true));

	/** Ellipsoid rotated a bit */
	Ellipsoid rotated = FitEllipsoid.fitTo(FitEllipsoid.testEllipsoid(7, 13,
			17, Math.PI / 4.32, 0, 0, 0, 0, 1000, true));

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

		assertArrayEquals(new double[] { 17, 13, 7 }, rotated.getRadii(), 1E-9);
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

		double[][] points4 = rotated.getSurfacePoints(1000);

		// inside
		for (double[] p : points4) {
			// contract by random fraction
			final double rand = Math.random();
			double x = p[0] * rand;
			double y = p[1] * rand;
			double z = p[2] * rand;
			System.out.println("Testing inside (" + x + ", " + y + ", " + z
					+ ") with random scaling = " + rand);
			assertTrue(rotated.contains(x, y, z));
		}

		// outside
		for (double[] p : points4) {
			// dilate by random fraction
			final double rand = Math.random() * 0.5;
			double x = p[0] / rand;
			double y = p[1] / rand;
			double z = p[2] / rand;
			System.out.println("Testing outside (" + x + ", " + y + ", " + z
					+ ") with random scaling = " + 1 / rand);
			assertTrue(!rotated.contains(x, y, z));
		}

		double[][] points3 = seventeenFiveThree.getSurfacePoints(1000);

		// inside
		for (double[] p : points3) {
			// contract by random fraction
			final double rand = Math.random();
			double x = p[0] * rand;
			double y = p[1] * rand;
			double z = p[2] * rand;
			System.out.println("Testing inside (" + x + ", " + y + ", " + z
					+ ") with random scaling = " + rand);
			assertTrue(seventeenFiveThree.contains(x, y, z));
		}

		// outside
		for (double[] p : points3) {
			// dilate by random fraction
			final double rand = Math.random() * 0.1;
			double x = p[0] / rand;
			double y = p[1] / rand;
			double z = p[2] / rand;
			System.out.println("Testing outside (" + x + ", " + y + ", " + z
					+ ") with random scaling = " + 1 / rand);
			assertTrue(!seventeenFiveThree.contains(x, y, z));
		}

		double[][] points2 = threeFiveSeventeen.getSurfacePoints(1000);

		// inside
		for (double[] p : points2) {
			// contract by random fraction
			final double rand = Math.random();
			double x = p[0] * rand;
			double y = p[1] * rand;
			double z = p[2] * rand;
			System.out.println("Testing inside (" + x + ", " + y + ", " + z
					+ ")");
			assertTrue(threeFiveSeventeen.contains(x, y, z));
		}

		// outside
		for (double[] p : points2) {
			// dilate by random fraction
			final double rand = Math.random() / 3;
			double x = p[0] / rand;
			double y = p[1] / rand;
			double z = p[2] / rand;
			System.out.println("Testing (" + x + ", " + y + ", " + z + ")");
			assertTrue(!threeFiveSeventeen.contains(x, y, z));
		}
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
		unitSphere.contract(0.015);
		double[][] points = unitSphere.getSurfacePoints(10000);
		for (double[] p : points) {
			assertEquals(0.985, Trig.distance3D(p), 1E-9);
		}
	}

	@Test
	public void testGetSortedRadii() {
		final int t = 1000;
		Ellipsoid e = unitSphere.copy();
		double[] r = new double[3];
		for (int i = 0; i < t; i++) {
			e.dilate(Math.random(), Math.random(), Math.random());
			r = e.getSortedRadii();
			assertTrue(r[0] < r[1]);
			assertTrue(r[1] < r[2]);
		}
	}

	@Test
	public void testGetEquation() {

		for (int i = 1; i < 1000; i++) {
			final double q = i * 0.1;
			final double a = q;
			final double b = Math.pow(q, 1.1);
			final double c = Math.pow(q, 1.5);
			final double x = i;
			final double y = i * 0.1;
			final double z = i;

			Object[] fit = FitEllipsoid.yuryPetrov(FitEllipsoid.testEllipsoid(
					a, b, c, Math.PI / 4.32, x, y, z, 0, 1000, true));
			Ellipsoid ellipsoid = new Ellipsoid(fit);
			double[] eq = ellipsoid.getEquation();
			double[] fitEq = (double[]) fit[3];
			assertArrayEquals(fitEq, eq, 1E-3);
		}
	}
}
