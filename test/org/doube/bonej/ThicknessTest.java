package org.doube.bonej;

import static org.junit.Assert.*;
import ij.ImagePlus;

import org.doube.geometry.TestDataMaker;
import org.doube.util.StackStats;
import org.junit.Test;

public class ThicknessTest {

	@Test
	public void testGetLocalThickness() {
		Thickness th = new Thickness();
		for (int d = 1; d < 25; d += 1) {
			ImagePlus rod = TestDataMaker.rod(d * 100, d);
			ImagePlus imp = th.getLocalThickness(rod, false);
			double[] stats = StackStats.meanStdDev(imp);
			assertEquals(d, stats[0], 1.5);
		}
	}

}
