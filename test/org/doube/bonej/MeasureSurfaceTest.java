package org.doube.bonej;

import static org.junit.Assert.assertEquals;

import java.util.List;

import javax.vecmath.Point3f;

import org.doube.geometry.TestDataMaker;
import org.junit.Test;

import ij.ImagePlus;
import marchingcubes.MCTriangulator;

public class MeasureSurfaceTest {

	@Test
	public void testGetSurfaceAreaOctahedron() {
		final double area = MeasureSurface.getSurfaceArea(TestDataMaker.octahedron());
		assertEquals(4.4558146404, area, 1e-6);
	}

	@Test
	public void testGetSurfaceAreaSphere() {
		final int r = 64;
		final ImagePlus imp = TestDataMaker.sphere(r);
		final MCTriangulator mct = new MCTriangulator();
		final boolean[] channels = { true, false, false };
		@SuppressWarnings("unchecked")
		final List<Point3f> points = mct.getTriangles(imp, 128, channels, 4);
		final double area = MeasureSurface.getSurfaceArea(points);
		assertEquals(4 * Math.PI * r * r, area, area * 0.05);
	}

	@Test
	public void testGetSurfaceAreaBox() {
		final int d = 128;
		final ImagePlus imp = TestDataMaker.brick(d, d, d);
		final MCTriangulator mct = new MCTriangulator();
		final boolean[] channels = { true, false, false };
		@SuppressWarnings("unchecked")
		final List<Point3f> points = mct.getTriangles(imp, 128, channels, 4);
		final double area = MeasureSurface.getSurfaceArea(points);
		assertEquals(6 * d * d, area, area * 0.02);
	}
}