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
	 * @param length
	 * @param diameter
	 * 
	 * @return
	 */
	public static ImagePlus rod(int length, int diameter) {
		ImageStack stack = new ImageStack(2 * diameter, 2 * diameter);
		for (int i = 0; i < length; i++) {
			ImageProcessor ip = new ByteProcessor(2 * diameter, 2 * diameter);
			ip.setColor(255);
			ip.fillOval((int) Math.floor(diameter / 2),
					(int) Math.floor(diameter / 2), diameter, diameter);
			stack.addSlice("" + i, ip);
		}
		ImagePlus imp = new ImagePlus("rod", stack);
		return imp;
	}

	/**
	 * Draw a solid sphere in the foreground, padded with 1 voxel of background
	 * on all stack faces
	 * 
	 * @param radius
	 * @return stack containing solid binary sphere
	 */
	public static ImagePlus sphere(int radius) {
		final int side = 2 * radius + 2;
		ImageStack stack = new ImageStack(side, side);
		ImageProcessor ip = new ByteProcessor(side, side);
		stack.addSlice("", ip); // padding
		for (int zd = -radius; zd <= radius; zd++) {
			final int rc = (int) Math.round(Math
					.sqrt(radius * radius - zd * zd));
			ImageProcessor ipr = new ByteProcessor(side, side);
			ipr.setColor(255);
			ipr.fillOval(radius + 1 - rc, radius + 1 - rc, 2 * rc, 2 * rc);
			stack.addSlice("", ipr);
		}
		ImageProcessor ipe = new ByteProcessor(side, side);
		stack.addSlice("", ipe); // padding
		ImagePlus imp = new ImagePlus("sphere", stack);
		return imp;
	}

	/**
	 * Create a binary brick of arbitrary width, height and depth, padded with 1
	 * voxel of background on all faces.
	 * 
	 * @param width
	 * @param height
	 * @param depth
	 * @return image with brick in foreground
	 */
	public static ImagePlus brick(int width, int height, int depth) {
		ImageStack stack = new ImageStack(width + 2, height + 2);
		ImageProcessor ip = new ByteProcessor(width + 2, height + 2);
		stack.addSlice("", ip);
		for (int i = 0; i < depth; i++) {
			ImageProcessor ip2 = new ByteProcessor(width + 2, height + 2);
			ip2.setColor(255);
			ip2.setRoi(1, 1, width, height);
			ip2.fill();
			stack.addSlice("", ip2);
		}
		ImageProcessor ip3 = new ByteProcessor(width + 2, height + 2);
		stack.addSlice("", ip3);
		ImagePlus imp = new ImagePlus("brick", stack);
		return imp;
	}

}
