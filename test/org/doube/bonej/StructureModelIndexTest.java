package org.doube.bonej;

import static org.junit.Assert.*;
import ij.ImagePlus;

import org.doube.geometry.TestDataMaker;
import org.junit.Test;

public class StructureModelIndexTest {

	@Test
	public void testHildRuegRod() {
		ImagePlus imp = TestDataMaker.rod(16384, 64);
		double smi = StructureModelIndex.hildRueg(imp, 6, 0.5f);
		assertEquals(3.0, smi, 0.03);
	}
	
	@Test
	public void testHildRuegPlate() {
		ImagePlus imp = TestDataMaker.brick(2048, 2048, 6);
		double smi = StructureModelIndex.hildRueg(imp, 6, 0.5f);
		assertEquals(0.0, smi, 0.05);
	}
	
	@Test
	public void testHildRuegSphere() {
		ImagePlus imp = TestDataMaker.sphere(256);
		double smi = StructureModelIndex.hildRueg(imp, 6, 0.5f);
		assertEquals(4.0, smi, 0.01);
	}
}
