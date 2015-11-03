package org.doube.bonej;

import static org.junit.Assert.assertEquals;

import org.doube.geometry.TestDataMaker;
import org.junit.Test;

import ij.ImagePlus;
import ij.measure.Calibration;

public class ConnectivityTest {

	private final Connectivity conn = new Connectivity();

	@Test
	public void testGetConnDensity() {
		final ImagePlus imp = TestDataMaker.boxFrame(32, 64, 128);
		final Calibration cal = imp.getCalibration();
		cal.pixelDepth = 0.2;
		cal.pixelHeight = 0.2;
		cal.pixelWidth = 0.2;
		// boxFrame adds 32 pixels of padding around the box
		final double stackVolume = (32 + 64) * (64 + 64) * (128 + 64) * (0.2 * 0.2 * 0.2);
		final double sumEuler = conn.getSumEuler(imp);
		final double deltaChi = conn.getDeltaChi(imp, sumEuler);
		final double connectivity = conn.getConnectivity(deltaChi);
		final double connD = conn.getConnDensity(imp, connectivity);
		assertEquals(5 / stackVolume, connD, 1e-12);
	}

	@Test
	public void testGetSumEulerCrossedCircle() {
		for (int size = 16; size < 1024; size *= 2) {
			final ImagePlus imp = TestDataMaker.crossedCircle(size);
			final double sumEuler = conn.getSumEuler(imp);
			assertEquals(-3, sumEuler, 1e-12);
		}
	}

	@Test
	public void testGetSumEulerBoxFrame() {
		for (int size = 16; size < 256; size *= 2) {
			final ImagePlus imp = TestDataMaker.boxFrame(size, size, size);
			final double sumEuler = conn.getSumEuler(imp);
			assertEquals(-4, sumEuler, 1e-12);
		}
	}
}
