package org.doube.bonej;

/**
 *  StructureModelIndex Copyright 2010 Michael Doube
 * 
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

import java.util.ArrayList;
import java.util.List;

import javax.media.j3d.Geometry;
import javax.media.j3d.GeometryArray;
import javax.vecmath.Color3f;
import javax.vecmath.Point3f;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
import ij3d.Image3DUniverse;
import isosurface.MeshEditor;

import marchingcubes.MCTriangulator;

import org.doube.bonej.Dilate;
import org.doube.util.ImageCheck;
import org.doube.util.ResultInserter;

import customnode.CustomTriangleMesh;

/**
 * Calculates the structure model index (SMI), a measure of how plate-like or
 * rod-like a structure is.
 * 
 * @author mdoube
 * @see <p>
 *      Hildebrand T, Rüegsegger P. Quantification of Bone Microarchitecture
 *      with the Structure Model Index. Comput Methods Biomech Biomed Engin
 *      1997;1(1):15-23. <a
 *      href="http://dx.doi.org/10.1080/01495739708936692">doi:
 *      10.1080/01495739708936692</a>
 *      </p>
 * 
 */
public class StructureModelIndex implements PlugIn {

	public void run(String arg) {
		ImagePlus imp = IJ.getImage();
		if (null == imp) {
			IJ.noImage();
			return;
		}
		ImageCheck ic = new ImageCheck();
		if (!ic.isBinary(imp)) {
			IJ.error("SMI needs a binary image.");
			return;
		}

		GenericDialog gd = new GenericDialog("Mesh Parameters");
		gd.addNumericField("Voxel resampling", 6, 0, 5, "voxels");
		gd.addNumericField("Mesh smoothing (0-1)", 0.5, 3, 5, "");
		gd.showDialog();
		if (gd.wasCanceled()) {
			return;
		}

		int voxelResampling = (int) Math.floor(gd.getNextNumber());
		float meshSmoothing = (float) gd.getNextNumber();

		// double smi = skyScan(imp, voxelResampling, meshSmoothing);
		double smi = hildRug(imp, voxelResampling, meshSmoothing);
		ResultInserter ri = ResultInserter.getInstance();
		ri.setResultInRow(imp, "SMI", smi);
		ri.updateTable();

		return;
	}

	/**
	 * <p>
	 * Calculate the SMI according to the SkyScan technical manual: make a
	 * surface mesh, get the surface area (s<sub>1</sub>) and volume (v), dilate
	 * the voxel model, get the new surface area (s<sub>2</sub>), calculate SMI
	 * as <br/>
	 * smi = 6 * (s<sub>2</sub>-s<sub>1</sub>)*v / s<sub>1</sub><sup>2</sup>
	 * </p>
	 * 
	 * @param imp
	 *            binary ImagePlus
	 * @param voxelResampling
	 *            resampling of voxels for smoother mesh creation. If set too
	 *            large will obliterate small features, if set too small, will
	 *            get too many surface triangles with jagging artefact.
	 * @param meshSmoothing
	 *            factor between 0 and 1 for smoothing of the surface mesh
	 * @return SMI. This is 0 for a plate, 3 for a rod and 4 for a sphere. SMI
	 *         can be negative if the structure contains many concave surfaces.
	 * @see <a href="http://www.skyscan.be/products/downloads.htm">The
	 *      description of measured parameters</a>
	 */
	@SuppressWarnings("unchecked")
	public static double skyScan(ImagePlus imp, int voxelResampling,
			float meshSmoothing) {
		int threshold = 128;
		final boolean[] channels = { true, false, false };
		MCTriangulator mct = new MCTriangulator();
		IJ.showStatus("Finding surface points...");
		List<Point3f> points = mct.getTriangles(imp, threshold, channels,
				voxelResampling);
		final Color3f colour = new Color3f(0.0f, 0.0f, 0.0f);
		IJ.showStatus("Creating surface mesh...");
		CustomTriangleMesh surface = new CustomTriangleMesh(points, colour,
				0.0f);
		IJ.showStatus("Smoothing surface mesh...");
		MeshEditor.smooth(surface, meshSmoothing);
		IJ.showStatus("Calculating volume...");
		double v = surface.getVolume();

		double s1 = MeasureSurface.getSurfaceArea(surface.getMesh());

		IJ.showStatus("Dilating voxel model...");
		Dilate d = new Dilate();
		ImagePlus imp2 = d.dilate(imp, 255);

		IJ.showStatus("Finding surface points...");
		points = mct.getTriangles(imp2, threshold, channels, voxelResampling);
		imp2.changes = false;
		imp2.close();
		IJ.showStatus("Creating surface mesh...");
		surface = new CustomTriangleMesh(points, colour, 0.0f);
		IJ.showStatus("Smoothing surface mesh...");
		MeshEditor.smooth(surface, meshSmoothing);
		double s2 = MeasureSurface.getSurfaceArea(surface.getMesh());
		double smi = 6 * ((s2 - s1) * v / (s1 * s1));
		IJ.showStatus("SMI calculated.");
		return smi;
	}

	@SuppressWarnings("unchecked")
	public static double hildRug(ImagePlus imp, int voxelResampling,
			float meshSmoothing) {
		int threshold = 128;
		final boolean[] channels = { true, false, false };
		MCTriangulator mct = new MCTriangulator();
		IJ.showStatus("Finding surface points...");
		List<Point3f> triangles = mct.getTriangles(imp, threshold, channels,
				voxelResampling);
		final Color3f colour = new Color3f(0.0f, 0.0f, 0.0f);
		IJ.showStatus("Creating surface mesh...");
		CustomTriangleMesh surface = new CustomTriangleMesh(triangles, colour,
				0.0f);
		IJ.showStatus("Smoothing surface mesh...");
		MeshEditor.smooth(surface, meshSmoothing);
		IJ.showStatus("Calculating volume...");
		double v = surface.getVolume();

		double s1 = MeasureSurface.getSurfaceArea(surface.getMesh());
		IJ.log("Original surface area is "+s1);
		GeometryArray ga = (GeometryArray) surface.getGeometry();
		float[] normals = new float[ga.getVertexCount()];
		ga.getNormals(0, normals);

		// get all the unique vertices
		ArrayList<Point3f> vertices = new ArrayList<Point3f>();
		final int nPoints = triangles.size();
		for (int p = 0; p < nPoints; p++) {
			Point3f testPoint = triangles.get(p);
			final int uniquePoints = vertices.size();
			boolean unique = true;
			for (int u = 0; u < uniquePoints; u++) {
				if (testPoint.equals(vertices.get(u))) {
					unique = false;
					break;
				}
			}
			if (unique) {
				vertices.add(testPoint);
			}
		}
		IJ.log("number of vertices = " + vertices.size());

		// associate each unique vertex with the triangles around it
		final int nVertices = vertices.size();
		ArrayList<ArrayList<Integer>> trianglesAtVertex = new ArrayList<ArrayList<Integer>>(nVertices);
		ArrayList<Integer> emptyArrayList = new ArrayList<Integer>();
		for (int i = 0; i < nVertices; i++){
			trianglesAtVertex.add(emptyArrayList);
		}
		final int nTriangles = nPoints / 3;
		IJ.log("nTriangles = " + nTriangles);
		for (int i = 0; i < nVertices; i++) {
			Point3f vertex = vertices.get(i);
			for (int p = 0; p < nPoints; p++) {
				Point3f testPoint = triangles.get(p);
				if (vertex.equals(testPoint)) {
					trianglesAtVertex.get(i).add(p);
				}
			}
		}

		// get the normals of the triangles around each vertex
		// and calculate the normal of the vertex as the mean triangle normal
		Point3f[] vertexNormals = new Point3f[nVertices];
		for (int i = 0; i < nVertices; i++) {
			final int vT = trianglesAtVertex.get(i).size();
			Point3f sumNormals = new Point3f();
			for (int j = 0; j < vT; j++) {
				final int pointIndex = trianglesAtVertex.get(i).get(j);
				final int corner = pointIndex % 3;
				Point3f point0 = new Point3f();
				Point3f point1 = new Point3f();
				Point3f point2 = new Point3f();
				switch (corner) {
				case 0:
					point0 = triangles.get(pointIndex);
					point1 = triangles.get(pointIndex + 1);
					point2 = triangles.get(pointIndex + 2);
					break;
				case 1:
					point0 = triangles.get(pointIndex - 1);
					point1 = triangles.get(pointIndex);
					point2 = triangles.get(pointIndex + 1);
					break;
				case 2:
					point0 = triangles.get(pointIndex - 2);
					point1 = triangles.get(pointIndex - 1);
					point2 = triangles.get(pointIndex);
					break;
				}
				Point3f surfaceNormal = crossProduct(point0, point1, point2);
				sumNormals.x += surfaceNormal.x;
				sumNormals.y += surfaceNormal.y;
				sumNormals.z += surfaceNormal.z;
			}
			vertexNormals[i] = new Point3f();
			vertexNormals[i].x = 100 * sumNormals.x / vT;
			vertexNormals[i].y = 100 * sumNormals.y / vT;
			vertexNormals[i].z = 100 * sumNormals.z / vT;
		}

		// move all the points by their normal
		ArrayList<Point3f> movedTriangles = new ArrayList<Point3f>();
		for (int t = 0; t < nPoints; t++) {
			Point3f point = triangles.get(t);
			Point3f newPoint = (Point3f) point.clone();
			for (int i = 0; i < nVertices; i++) {
				if (point.equals(vertices.get(i))) {
					newPoint.x += vertexNormals[i].x;
					newPoint.y += vertexNormals[i].y;
					newPoint.z += vertexNormals[i].z;
					movedTriangles.add(newPoint);
					break;
				}
			}
		}
		CustomTriangleMesh surface2 = new CustomTriangleMesh(triangles, colour,
				0.0f);
		IJ.showStatus("Smoothing surface mesh...");
		MeshEditor.smooth(surface, meshSmoothing);
		IJ.showStatus("Calculating volume...");

		double s2 = MeasureSurface.getSurfaceArea(surface2.getMesh());
		IJ.log("Dilated surface area is "+s2);

		// Add the mesh
		Image3DUniverse univ = new Image3DUniverse();
		univ.show();
		Color3f red = new Color3f(1.0f, 0.0f, 0.0f);
		Color3f green = new Color3f(0.0f, 1.0f, 0.0f);
		try {
			univ.addTriangleMesh(triangles, green, "Original").setLocked(true);
		} catch (NullPointerException npe) {
			IJ.log("3D Viewer was closed before rendering completed.");
		}
		try {
			univ.addTriangleMesh(movedTriangles, red, "Dilated").setLocked(
					true);
		} catch (NullPointerException npe) {
			IJ.log("3D Viewer was closed before rendering completed.");
		}

		return 0;
	}

	private static Point3f crossProduct(Point3f point0, Point3f point1,
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
