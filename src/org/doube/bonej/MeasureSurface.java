package org.doube.bonej;

import java.util.List;

import javax.vecmath.Color3f;
import javax.vecmath.Point3f;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
import ij3d.Content;
import ij3d.Image3DUniverse;

import marchingcubes.MCTriangulator;
import customnode.CustomTriangleMesh;

/**
 * Make a mesh from a binary image and get some measurements from it.
 * 
 * @author Michael Doube
 * 
 */
public class MeasureSurface implements PlugIn {

	@SuppressWarnings("unchecked")
	public void run(String arg) {
		ImagePlus imp = IJ.getImage();
		if (null == imp) {
			IJ.noImage();
			return;
		}
		ImageCheck ic = new ImageCheck();
		if (!ic.isBinary(imp)) {
			IJ.error("Mesh creation requires a binary image.");
			return;
		}

		GenericDialog gd = new GenericDialog("Options");
		gd.addNumericField("Resampling", 6, 0);
		gd.addCheckbox("Show surface", true);
		gd.showDialog();
		int resamplingF = (int) Math.floor(gd.getNextNumber());
		boolean doSurfaceRendering = gd.getNextBoolean();

		final int threshold = 128;
		final boolean[] channels = { true, false, false };

		MCTriangulator mct = new MCTriangulator();
		List<Point3f> points = mct.getTriangles(imp, threshold, channels,
				resamplingF);

		IJ.log("Isosurface contains " + (points.size() / 3) + " triangles");
		
		ResultInserter ri = new ResultInserter();
		double area = getSurfaceArea(points);
		ri.setResultInRow(imp,
				"BS (" + imp.getCalibration().getUnits() + "^2)", area);

		if (doSurfaceRendering) {
			renderSurface(points, "Surface of " + imp.getTitle());
		}
		points.clear();
		ri.updateTable();
		return;
	}

	/**
	 * Show the surface in the ImageJ 3D Viewer
	 * 
	 * @param points
	 * @param title
	 */
	private void renderSurface(List<Point3f> points, String title) {
		CustomTriangleMesh mesh = new CustomTriangleMesh(points);

		// Create a universe and show it
		Image3DUniverse univ = new Image3DUniverse();
		univ.show();

		// Add the mesh
		Content c = univ.addCustomMesh(mesh, title);
		Color3f green = new Color3f(0.0f, 0.5f, 0.0f);
		c.getColor();
		c.getColor();
		c.setColor(green);
		c.setTransparency((float) 0.33);
		c.setSelected(true);
	}

	/**
	 * Calculate surface area of the isosurface
	 * 
	 * @param points
	 * @return
	 */
	private double getSurfaceArea(List<Point3f> points) {
		double sumArea = 0;
		for (int n = 0; n < points.size(); n += 3) {
			Point3f point0 = points.get(n);
			Point3f point1 = points.get(n + 1);
			Point3f point2 = points.get(n + 2);

			// TODO reject triangle and continue if it is flush
			// with a cut face / image side

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
