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
		for (int r = 2; r < 25; r += 1) {
			ImagePlus sphere = TestDataMaker.sphere(r);
			ImagePlus imp = th.getLocalThickness(sphere, false);
			double[] stats = StackStats.meanStdDev(imp);
			System.out.print(r+", "+stats[0]+"\n");
			double regression = r*1.9441872882 - 1.218936;
			assertEquals(regression, stats[0], regression * 0.1);
		}
	}

}
