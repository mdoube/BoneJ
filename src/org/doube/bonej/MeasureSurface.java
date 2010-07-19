package org.doube.bonej;

/**
 * MeasureSurface plugin for ImageJ
 * Copyright 2009 2010 Michael Doube
 * 
 *This program is free software: you can redistribute it and/or modify
 *it under the terms of the GNU General Public License as published by
 *the Free Software Foundation, either version 3 of the License, or
 *(at your option) any later version.
 *
 *This program is distributed in the hope that it will be useful,
 *but WITHOUT ANY WARRANTY; without even the implied warranty of
 *MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *GNU General Public License for more details.
 *
 *You should have received a copy of the GNU General Public License
 *along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

import java.util.List;

import javax.vecmath.Color3f;
import javax.vecmath.Point3f;

import org.doube.geometry.Vectors;
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
		if (!ImageCheck.checkEnvironment())
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
		gd.addHelp("http://bonej.org/isosurface");
		gd.showDialog();
		int resamplingF = (int) Math.floor(gd.getNextNumber());
		threshold = (int) Math.floor(gd.getNextNumber());
		boolean doSurfaceRendering = gd.getNextBoolean();
		if (gd.wasCanceled())
			return;

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
		Point3f origin = new Point3f(0.0f, 0.0f, 0.0f);
		for (int n = 0; n < nPoints; n += 3) {
			IJ.showStatus("Calculating surface area...");
			final Point3f point0 = points.get(n);
			final Point3f point1 = points.get(n + 1);
			final Point3f point2 = points.get(n + 2);

			// TODO reject triangle and continue if it is flush
			// with a cut face / image side

			// area of triangle is half magnitude
			// of cross product of 2 edge vectors
			Point3f cp = Vectors.crossProduct(point0, point1, point2);

			final double deltaArea = 0.5 * cp.distance(origin);

			sumArea += deltaArea;
		}
		return sumArea;
	}
}
