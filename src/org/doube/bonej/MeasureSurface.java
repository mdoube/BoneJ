package org.doube.bonej;

import java.util.List;

import javax.vecmath.Color3f;
import javax.vecmath.Point3f;

import org.doube.geometry.Trig;
import org.doube.geometry.VectorProducts;
import org.doube.util.ImageCheck;
import org.doube.util.ResultInserter;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
import ij3d.Content;
import ij3d.Image3DUniverse;

import marchingcubes.MCTriangulator;
import customnode.CustomTriangleMesh;

/**
 * Make a mesh from a binary or 8-bit image and get surface area measurements
 * from it.
 * 
 * @author Michael Doube
 * 
 */
public class MeasureSurface implements PlugIn {

	@SuppressWarnings("unchecked")
	public void run(String arg) {
		if (!ImageCheck.checkIJVersion())
			return;
		ImagePlus imp = IJ.getImage();
		if (null == imp) {
			IJ.noImage();
			return;
		}
		int threshold = 128;
		ImageCheck ic = new ImageCheck();
		if (ic.isBinary(imp)) {
			threshold = 128;
		} else if (imp.getBitDepth() == 8) {
			ThresholdMinConn tmc = new ThresholdMinConn();
			threshold = imp.getProcessor().getAutoThreshold(
					tmc.getStackHistogram(imp));
		} else {
			IJ.error("Isosurface", "Image type not supported");
			return;
		}

		GenericDialog gd = new GenericDialog("Options");
		gd.addNumericField("Resampling", 6, 0);
		gd.addNumericField("Threshold", threshold, 0);
		gd.addCheckbox("Show surface", true);
		gd.showDialog();
		int resamplingF = (int) Math.floor(gd.getNextNumber());
		threshold = (int) Math.floor(gd.getNextNumber());
		boolean doSurfaceRendering = gd.getNextBoolean();

		final boolean[] channels = { true, false, false };

		MCTriangulator mct = new MCTriangulator();
		List<Point3f> points = mct.getTriangles(imp, threshold, channels,
				resamplingF);

		IJ.log("Isosurface contains " + (points.size() / 3) + " triangles");

		ResultInserter ri = ResultInserter.getInstance();
		double area = getSurfaceArea(points);
		ri.setResultInRow(imp, "BS (" + imp.getCalibration().getUnits() + "²)",
				area);
		ri.updateTable();

		if (points.size() == 0) {
			IJ.error("Isosurface contains no points");
			return;
		}

		if (doSurfaceRendering) {
			renderSurface(points, "Surface of " + imp.getTitle());
		}
		points.clear();
		return;
	}

	/**
	 * Show the surface in the ImageJ 3D Viewer
	 * 
	 * @param points
	 * @param title
	 */
	private void renderSurface(List<Point3f> points, String title) {
		IJ.showStatus("Rendering surface...");
		CustomTriangleMesh mesh = new CustomTriangleMesh(points);

		// Create a universe and show it
		Image3DUniverse univ = new Image3DUniverse();
		univ.show();

		// Add the mesh
		Content c = univ.addCustomMesh(mesh, title);
		Color3f green = new Color3f(0.0f, 0.5f, 0.0f);
		c.getColor();
		c.setColor(green);
		c.setTransparency((float) 0.33);
		c.setSelected(true);
	}

	/**
	 * Calculate surface area of the isosurface
	 * 
	 * @param points
	 *            in 3D triangle mesh
	 * @return surface area
	 */
	public static double getSurfaceArea(List<Point3f> points) {
		double sumArea = 0;
		final int nPoints = points.size();
		for (int n = 0; n < nPoints; n += 3) {
			IJ.showStatus("Calculating surface area...");
			final Point3f point0 = points.get(n);
			final Point3f point1 = points.get(n + 1);
			final Point3f point2 = points.get(n + 2);

			// TODO reject triangle and continue if it is flush
			// with a cut face / image side

			// area of triangle is half magnitude
			// of cross product of 2 edge vectors
			Point3f cp = VectorProducts.crossProduct(point0, point1, point2);

			final double deltaArea = 0.5 * Trig.distance3D(cp);

			sumArea += deltaArea;
		}
		return sumArea;
	}
}
