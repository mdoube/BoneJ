package org.doube.geometry;

import static org.junit.Assert.assertArrayEquals;

import org.junit.Test;

public class FitEllipseTest {
	final static double semiMinor = 0.861310905613;
	final static double semiMajor = 2.902552357932;
	final static double theta = 1.925109462859;
	final static double x = 1.0;
	final static double y = -3.5;
	final double[] expected = { x, y, semiMinor, semiMajor, theta };
	final static int n = 100;
	final static double[][] points = FitEllipse.testEllipse(semiMinor, semiMajor, theta, x, y, 0.0, n);
	final static double[][] noisyPoints = FitEllipse.testEllipse(semiMinor, semiMajor, theta, x, y, 0.01, n);

	@Test
	public void testDirect() {
		assertArrayEquals(expected, FitEllipse.varToDimensions(FitEllipse.direct(points)), 1e-10);
		assertArrayEquals(expected, FitEllipse.varToDimensions(FitEllipse.direct(noisyPoints)), 1e-2);
	}

	// @Test
	// public void testTaubin() {
	// assertArrayEquals(expected,
	// FitEllipse.varToDimensions(FitEllipse.taubin(points)), 1e-10);
	// assertArrayEquals(expected,
	// FitEllipse.varToDimensions(FitEllipse.taubin(noisyPoints)), 1e-10);
	// }

	@Test
	public void testTestEllipse() {
		// ellipse with minor and major axes, centred on (0,0)
		double[][] p = FitEllipse.testEllipse(2.3, 5.7, 0.0, 0.0, 0.0, 0.0, 4);
		final double[][] expectedPoints = { { 2.3, 0.0 }, { 0.0, 5.7 }, { -2.3, 0.0 }, { 0.0, -5.7 } };
		for (int i = 0; i < 4; i++) {
			assertArrayEquals(expectedPoints[i], p[i], 1e-10);
		}

		// off-centre ellipse
		p = FitEllipse.testEllipse(2.3, 5.7, 0.0, 8.1, -7.2, 0.0, 4);
		final double[][] expectedPoints2 = { { 10.4, -7.2 }, { 8.1, -1.5 }, { 5.8, -7.2 }, { 8.1, -12.9 } };
		for (int i = 0; i < 4; i++) {
			assertArrayEquals(expectedPoints2[i], p[i], 1e-10);
		}

	}

	@Test
	public void testVarToDimensions() {
		// a, b, c, d, f and g
		final double[] p = { 9.0, 6.0, 2.0, 3.0, 8.0, 5.0 };
		final double[] d = { 1.0, -3.5, 2.902552357932, 0.861310905613, 1.925109462859 };
		assertArrayEquals(d, FitEllipse.varToDimensions(p), 1e-10);
	}

}
