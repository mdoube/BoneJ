package org.doube.geometry;

import ij.ImagePlus;
import ij.ImageStack;
import org.doube.skeleton.Skeletonize3D;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;

import java.util.ArrayList;
import java.util.List;

import javax.vecmath.Point3f;

/**
 * Static methods to generate images for testing
 * 
 * @author Michael Doube
 */
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

	/**
	 * Draw a circle with vertical and horizontal crossing, then skeletonize it
	 * 
	 * @param size
	 *            width and height of the image, circle diameter is size/2
	 * @return image containing a white (255) circle on black (0) background
	 */
	public static ImagePlus crossedCircle(int size) {
		ImageProcessor ip = new ByteProcessor(size, size);
		ip.setColor(0);
		ip.fill();
		ip.setColor(255);
		ip.drawOval(size / 4, size / 4, size / 2, size / 2);
		ip.drawLine(size / 2, size / 4, size / 2, 3 * size / 4);
		ip.drawLine(size / 4, size / 2, 3 * size / 4, size / 2);
		ImagePlus imp = new ImagePlus("crossed-circle", ip);
		Skeletonize3D skel = new Skeletonize3D();
		return skel.getSkeleton(imp);
	}

	/**
	 * Draw the edges of a brick with 32 pixels of padding on all faces
	 * 
	 * @param width
	 *            Width of the box frame in pixels
	 * @param height
	 *            Height of the box frame in pixels
	 * @param depth
	 *            Depth of the box frame in pixels
	 * @return Image containing a 1-pixel wide outline of a 3D box
	 */
	public static ImagePlus boxFrame(int width, int height, int depth) {
		ImageStack stack = new ImageStack(width + 64, height + 64);
		for (int s = 1; s <= depth + 64; s++) {
			ImageProcessor ip = new ByteProcessor(width + 64, height + 64);
			ip.setColor(0);
			ip.fill();
			stack.addSlice(ip);
		}
		ImageProcessor ip = stack.getProcessor(32);
		ip.setColor(255);
		ip.drawRect(32, 32, width, height);
		ip = stack.getProcessor(32 + depth);
		ip.setColor(255);
		ip.drawRect(32, 32, width, height);
		for (int s = 33; s < 32 + depth; s++) {
			ip = stack.getProcessor(s);
			ip.setColor(255);
			ip.drawPixel(32, 32);
			ip.drawPixel(32, 31 + height);
			ip.drawPixel(31 + width, 32);
			ip.drawPixel(31 + width, 31 + height);
		}
		ImagePlus imp = new ImagePlus("box-frame", stack);
		return imp;
	}

	public static ImagePlus binaryNoise(int width, int height, int depth,
			double ratio) {
		final int npixels = width * height;
		ImageStack stack = new ImageStack(width, height);
		for (int i = 0; i < depth; i++) {
			ByteProcessor bp = new ByteProcessor(width, height);
			for (int index = 0; index < npixels; index ++){
				double random = Math.random();
				if (random > ratio)
					bp.set(index, 255);
			}
			stack.addSlice(bp);
		}
		ImagePlus imp = new ImagePlus("binary-noise", stack);
		return imp;
	}

	public static ImagePlus plates(int width, int height, int depth, int spacing) {
		ImageStack stack = new ImageStack(width, height);
		for (int i = 0; i < depth; i++)
			stack.addSlice(new ByteProcessor(width, height));
		
		for (int i = 1; i <= depth; i += spacing){
			ByteProcessor bp = (ByteProcessor) stack.getProcessor(i);
			bp.add(255);
		}
		ImagePlus imp = new ImagePlus("plates", stack);
		return imp; 
	}
}
