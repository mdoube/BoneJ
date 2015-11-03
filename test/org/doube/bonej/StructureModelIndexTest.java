package org.doube.bonej;

import static org.junit.Assert.assertEquals;

import org.doube.geometry.TestDataMaker;
import org.junit.Test;

import ij.ImagePlus;

public class StructureModelIndexTest {

	@Test
	public void testHildRuegRod() {
		final ImagePlus imp = TestDataMaker.rod(16384, 64);
		final double smi = StructureModelIndex.hildRueg(imp, 6, 0.5f);
		assertEquals(3.0, smi, 0.03);
	}

	@Test
	public void testHildRuegPlate() {
		final ImagePlus imp = TestDataMaker.brick(2048, 2048, 6);
		final double smi = StructureModelIndex.hildRueg(imp, 6, 0.5f);
		assertEquals(0.0, smi, 0.05);
	}

	@Test
	public void testHildRuegSphere() {
		final ImagePlus imp = TestDataMaker.sphere(256);
		final double smi = StructureModelIndex.hildRueg(imp, 6, 0.5f);
		assertEquals(4.0, smi, 0.01);
	}
}
