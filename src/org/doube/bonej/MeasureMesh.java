package org.doube.bonej;

import java.util.List;

import javax.vecmath.Point3f;

import ij.IJ;
import ij.ImagePlus;
import ij.plugin.PlugIn;

import marchingcubes.MCTriangulator;

/**
 * Make a mesh from a binary image and get some measurements
 * from it.  
 * 
 * @author Michael Doube
 *
 */
public class MeasureMesh implements PlugIn {

	public void run(String arg) {
		ImagePlus imp = IJ.getImage();
		if (null == imp) {
			IJ.noImage();
			return;
		}

		int threshold = 128;
		boolean[] channels = { true, true, true };
		int resamplingF = 2;

		MCTriangulator mct = new MCTriangulator();
		List<Point3f> points = mct.getTriangles(imp, threshold, channels,
				resamplingF);

		double area = getSurfaceArea(points);
		IJ.log("Area of mesh is " + area + " "
				+ imp.getCalibration().getUnits());
	}

	private double getSurfaceArea(List<Point3f> points) {
		double sumArea = 0;
		for (int n = 0; n < points.size(); n += 3) {
			Point3f point0 = points.get(n);
			Point3f point1 = points.get(n + 1);
			Point3f point2 = points.get(n + 2);

			float x1 = point1.x - point0.x;
			float y1 = point1.y - point0.y;
			float z1 = point1.z - point0.z;

			float x2 = point2.x - point0.x;
			float y2 = point2.y - point0.y;
			float z2 = point2.z - point0.z;

			// area of triangle is half magnitude
			// of cross product of 2 edge vectors

			double x = y1 * z2 - z1 * y2;
			double y = z1 * x2 - x1 * z2;
			double z = x1 * y2 - y1 * x2;

			double deltaArea = 0.5 * Math.sqrt(x * x + y * y + z * z);

			sumArea += deltaArea;
		}
		return sumArea;

	}
}
