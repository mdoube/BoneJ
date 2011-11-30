package org.doube.geometry;

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
}
