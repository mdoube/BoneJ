package org.doube.geometry;

import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;

import java.util.ArrayList;
import java.util.List;

import javax.vecmath.Point3f;

public class TestDataMaker {

	/**
	 * An octahedron with vertices defined by points using values of 2, 1.198039
	 * and 2.801961.
	 * 
	 * It has a calculated surface area of 4.4558146404 units^2
	 * 
	 * @return points defining the vertices of the octahedron's triangles
	 */
	static public List<Point3f> octahedron() {
		final float a = 2.0f;
		final float b = 1.198039f;
		final float c = 2.801961f;
		List<Point3f> points = new ArrayList<Point3f>();
		points.add(new Point3f(a, a, b));
		points.add(new Point3f(b, a, a));
		points.add(new Point3f(a, b, a));
		points.add(new Point3f(a, a, b));
		points.add(new Point3f(a, c, a));
		points.add(new Point3f(b, a, a));
		points.add(new Point3f(a, b, a));
		points.add(new Point3f(c, a, a));
		points.add(new Point3f(a, a, b));
		points.add(new Point3f(c, a, a));
		points.add(new Point3f(a, c, a));
		points.add(new Point3f(a, a, b));
		points.add(new Point3f(a, b, a));
		points.add(new Point3f(b, a, a));
		points.add(new Point3f(a, a, c));
		points.add(new Point3f(b, a, a));
		points.add(new Point3f(a, c, a));
		points.add(new Point3f(a, a, c));
		points.add(new Point3f(a, b, a));
		points.add(new Point3f(a, a, c));
		points.add(new Point3f(c, a, a));
		points.add(new Point3f(c, a, a));
		points.add(new Point3f(a, a, c));
		points.add(new Point3f(a, c, a));
		return points;
	}

	/**
	 * Generate a rod of circular cross-section. The rod is oriented with its
	 * axis in the z direction and in a stack of 2*diameter wide and high.
	 * 
	 * @param diameter
	 * @param length
	 * 
	 * @return
	 */
	public static ImagePlus rod(int length, int diameter) {
		ImageStack stack = new ImageStack(2 * diameter, 2 * diameter);
		for (int i = 0; i < length; i++) {
			ImageProcessor ip = new ByteProcessor(2 * diameter, 2 * diameter);
			ip.setColor(255);
			ip.fillOval(diameter / 2, diameter / 2, diameter, diameter);
			stack.addSlice("" + i, ip);
		}
		ImagePlus imp = new ImagePlus("rod", stack);
		return imp;
	}
}
