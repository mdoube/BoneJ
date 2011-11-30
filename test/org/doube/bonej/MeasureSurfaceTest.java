package org.doube.bonej;

import static org.junit.Assert.*;

import org.doube.geometry.TestDataMaker;
import org.junit.Test;

public class MeasureSurfaceTest {

	@Test
	public void testRun() {
		fail("Not yet implemented"); // TODO
	}

	@Test
	public void testGetSurfaceArea() {
		double area = MeasureSurface.getSurfaceArea(TestDataMaker.octahedron());
		assertEquals(4.4558146404, area, 1e-6);
	}

}