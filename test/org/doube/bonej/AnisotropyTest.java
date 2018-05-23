package org.doube.bonej;

import static org.junit.Assert.assertEquals;

import org.doube.geometry.TestDataMaker;
import org.junit.Test;

import ij.ImagePlus;

import java.util.Random;

public class AnisotropyTest {

	private final Anisotropy anisotropy = new Anisotropy();

	@Test
	public void testRunToStableResultIsotropy() {
		final Random random = new Random(1234);
		final ImagePlus imp = TestDataMaker.binaryNoise(256, 256, 256, 0.25, random);
		final Object[] result = anisotropy.runToStableResult(imp, 100, 2000, 50000, 256 / 4, 2.3, 0.005, true);
		final double da = ((double[]) result[0])[0];
		assertEquals(0, da, 1e-2);
	}

	@Test
	public void testRunToStableResultAnisotropy() {
		final ImagePlus imp = TestDataMaker.plates(256, 256, 256, 8);
		final Object[] result = anisotropy.runToStableResult(imp, 100, 2000, 50000, 256 / 4, 2.3, 0.005, true);
		final double da = ((double[]) result[0])[0];
		assertEquals(1, da, 1e-12);
	}

	@Test
	public void testCalculateSingleSphereIsotropy() {
		final Random random = new Random(12345);
		final ImagePlus imp = TestDataMaker.binaryNoise(256, 256, 256, 0.25, random);
		final double[] centroid = { 128, 128, 128 };
		final Object[] result = anisotropy.calculateSingleSphere(imp, centroid, 127, 2.3, 50000, false);
		final double da = ((double[]) result[0])[0];
		assertEquals(0, da, 1e-1);
	}

	@Test
	public void testCalculateSingleSphereAnisotropy() {
		final ImagePlus imp = TestDataMaker.plates(256, 256, 256, 8);
		final double[] centroid = { 128, 128, 128 };
		final Object[] result = anisotropy.calculateSingleSphere(imp, centroid, 127, 2.3, 50000, false);
		final double da = ((double[]) result[0])[0];
		assertEquals(1, da, 1e-2);
	}

}
