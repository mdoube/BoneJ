package org.doube.geometry;

import static org.junit.Assert.assertArrayEquals;

import org.junit.Test;

public class FitSphereTest {

	@Test
	public void testFitSphere() {
		for (double r = 1; r < 5; r += Math.random()) {
			for (double x = -2; x < 2; x += Math.random()) {
				for (double y = -2; y < 2; y += Math.random()) {
					for (double z = -2; z < 2; z += Math.random()) {
						for (double theta = -Math.PI; theta < Math.PI; theta += Math.random()) {
							final double[][] points = FitEllipsoid.testEllipsoid(r, r, r, theta, x, y, z, 0.001, 10,
									true);
							final double[] sphere = FitSphere.fitSphere(points);
							final double[] expected = { x, y, z, r };
							assertArrayEquals(expected, sphere, 1e-2);
						}
					}
				}
			}
		}
	}

}
