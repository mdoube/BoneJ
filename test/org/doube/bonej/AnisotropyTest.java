package org.doube.bonej;

import static org.junit.Assert.*;
import ij.ImagePlus;

import org.doube.geometry.TestDataMaker;
import org.junit.Test;

public class AnisotropyTest {

	private final Anisotropy anisotropy = new Anisotropy();

	@Test
	public void testRunToStableResultIsotropy() {
		ImagePlus imp = TestDataMaker.binaryNoise(256, 256, 256, 0.25);
		Object[] result = anisotropy.runToStableResult(imp, 100, 2000, 50000,
				256 / 4, 2.3, 0.005, true);
		double da = ((double[])result[0])[0];
		assertEquals(0, da, 1e-2);
	}
	
	@Test
	public void testRunToStableResultAnisotropy() {
		ImagePlus imp = TestDataMaker.plates(256, 256, 256, 8);
		imp.show();
		Object[] result = anisotropy.runToStableResult(imp, 100, 2000, 50000,
				256 / 4, 2.3, 0.005, true);
		double da = ((double[])result[0])[0];
		assertEquals(1, da, 1e-12);
	}

	@Test
	public void testCalculateSingleSphere() {
		fail("Not yet implemented"); // TODO
	}

}
