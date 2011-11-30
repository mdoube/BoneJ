package org.doube.bonej;

import static org.junit.Assert.*;
import ij.ImagePlus;

import org.doube.geometry.TestDataMaker;
import org.junit.Test;

public class StructureModelIndexTest {

	@Test
	public void testHildRuegRod() {
		int length = 16384;
		int diameter = 64;
		ImagePlus imp = TestDataMaker.rod(length, diameter);
		double smi = StructureModelIndex.hildRueg(imp, 6, 0.5f);
		assertEquals(3.0, smi, 0.03);
	}
	
	@Test
	public void testHildRuegPlate() {
		//TODO not yet implemented
	}
	
	@Test
	public void testHildRuegSphere() {
		ImagePlus imp = TestDataMaker.sphere(256);
		double smi = StructureModelIndex.hildRueg(imp, 6, 0.5f);
		assertEquals(4.0, smi, 0.01);
	}
}
