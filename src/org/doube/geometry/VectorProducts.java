package org.doube.geometry;

import javax.vecmath.Point3f;

public class VectorProducts {
	/**
	 * Calculate the cross product of 3 Point3f's, which describe two vectors
	 * joined at the tails. Can be used to find the plane / surface normal of a
	 * triangle. Half of its magnitude is the area of the triangle.
	 * 
	 * @param point0
	 *            both vectors' tails
	 * @param point1
	 *            vector 1's head
	 * @param point2
	 *            vector 2's head
	 * @return
	 */
	public static Point3f crossProduct(Point3f point0, Point3f point1,
			Point3f point2) {
		final float x1 = point1.x - point0.x;
		final float y1 = point1.y - point0.y;
		final float z1 = point1.z - point0.z;

		final float x2 = point2.x - point0.x;
		final float y2 = point2.y - point0.y;
		final float z2 = point2.z - point0.z;

		Point3f crossVector = new Point3f();
		crossVector.x = y1 * z2 - z1 * y2;
		crossVector.y = z1 * x2 - x1 * z2;
		crossVector.z = x1 * y2 - y1 * x2;

		return crossVector;
	}

}
