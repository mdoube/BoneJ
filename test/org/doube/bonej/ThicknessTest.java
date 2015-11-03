package org.doube.bonej;

import static org.junit.Assert.assertEquals;

import org.doube.geometry.TestDataMaker;
import org.doube.util.StackStats;
import org.junit.Test;

import ij.ImagePlus;

public class ThicknessTest {

	@Test
	public void testGetLocalThicknessRod() {
		final Thickness th = new Thickness();
		for (int d = 1; d < 25; d += 1) {
			final ImagePlus rod = TestDataMaker.rod(d * 100, d);
			final ImagePlus imp = th.getLocalThickness(rod, false);
			final double[] stats = StackStats.meanStdDev(imp);
			System.out.print(d + ", " + stats[0] + "\n");
			assertEquals(d, stats[0], 1.5);
		}
	}

	@Test
	public void testGetLocalThicknessSphere() {
		final Thickness th = new Thickness();
		for (int r = 2; r < 25; r++) {
			final ImagePlus sphere = TestDataMaker.sphere(r);
			final ImagePlus imp = th.getLocalThickness(sphere, false);
			final double[] stats = StackStats.meanStdDev(imp);
			final double regression = r * 1.9441872882 - 1.218936;
			System.out.print(r * 2 + ", " + stats[0] + "\n");
			assertEquals(regression, stats[0], regression * 0.1);
		}
	}

	@Test
	public void testGetLocalThicknessBrick() {
		final Thickness th = new Thickness();
		for (int t = 1; t < 21; t++) {
			final ImagePlus brick = TestDataMaker.brick(128, 128, t);
			final ImagePlus imp = th.getLocalThickness(brick, false);
			final double[] stats = StackStats.meanStdDev(imp);
			int expected = t;
			// pixelation and *2 (radius to diameter conversion) weirdness
			if (t % 2 != 0)
				expected++;
			System.out.print(t + ", " + stats[0] + "\n");
			assertEquals(expected, stats[0], expected * 0.05);
		}
	}

}
