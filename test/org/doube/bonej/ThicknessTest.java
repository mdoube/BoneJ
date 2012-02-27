package org.doube.bonej;

import static org.junit.Assert.*;
import ij.ImagePlus;

import org.doube.geometry.TestDataMaker;
import org.doube.util.StackStats;
import org.junit.Test;

public class ThicknessTest {

	@Test
	public void testGetLocalThicknessRod() {
		Thickness th = new Thickness();
		for (int d = 1; d < 25; d += 1) {
			ImagePlus rod = TestDataMaker.rod(d * 100, d);
			ImagePlus imp = th.getLocalThickness(rod, false);
			double[] stats = StackStats.meanStdDev(imp);
			assertEquals(d, stats[0], 1.5);
		}
	}

	@Test
	public void testGetLocalThicknessSphere() {
		Thickness th = new Thickness();
		for (int r = 2; r < 25; r++) {
			ImagePlus sphere = TestDataMaker.sphere(r);
			ImagePlus imp = th.getLocalThickness(sphere, false);
			double[] stats = StackStats.meanStdDev(imp);
			double regression = r * 1.9441872882 - 1.218936;
			assertEquals(regression, stats[0], regression * 0.1);
		}
	}

	@Test
	public void testGetLocalThicknessBrick() {
		Thickness th = new Thickness();
		for (int t = 1; t < 21; t++) {
			ImagePlus brick = TestDataMaker.brick(128, 128, t);
			ImagePlus imp = th.getLocalThickness(brick, false);
			double[] stats = StackStats.meanStdDev(imp);
			int expected = t;
			//pixelation and *2 (radius to diameter conversion) weirdness
			if (t%2 != 0)
				expected++;
			System.out.print(t + ", " + stats[0] + "\n");
			assertEquals(expected, stats[0], expected * 0.05);
		}
	}

}
