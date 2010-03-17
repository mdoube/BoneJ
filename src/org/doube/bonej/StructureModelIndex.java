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
import java.util.Hashtable;
import java.util.List;

import javax.vecmath.Color3f;
import javax.vecmath.Point3f;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
import isosurface.MeshEditor;

import marchingcubes.MCTriangulator;

import org.doube.bonej.Dilate;
import org.doube.geometry.Trig;
import org.doube.geometry.VectorProduct;
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
		if (!ImageCheck.checkIJVersion()) {
			return;
		}
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
		String[] smiMethods = { "Hildebrand & Rüegsegger", "SkyScan" };
		gd.addChoice("SMI Method", smiMethods, "Hildebrand & Rüegsegger");
		gd.addNumericField("Voxel resampling", 6, 0, 5, "voxels");
		gd.addNumericField("Mesh smoothing (0-1)", 0.5, 3, 5, "");
		gd.showDialog();
		if (gd.wasCanceled()) {
			return;
		}

		String smiMethod = gd.getNextChoice();
		int voxelResampling = (int) Math.floor(gd.getNextNumber());
		float meshSmoothing = (float) gd.getNextNumber();

		double smi = 0;
		if (smiMethod.equals(smiMethods[0])) {
			smi = hildRueg(imp, voxelResampling, meshSmoothing);
		} else {
			smi = skyScan(imp, voxelResampling, meshSmoothing);
		}
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

	/**
	 * <p>
	 * Calculate the structure model index according to the description by
	 * Hildebrand and Rugsegger. Creates a surface model, dilates it by a small
	 * increment and compares the areas before and after dilation.
	 * </p>
	 * 
	 * <p>
	 * This method is preferred to voxel-dilating methods, like skyScan() in
	 * this class
	 * </p>
	 * 
	 * @see <p>
	 *      Hildebrand T, Rüegsegger P. Quantification of Bone Microarchitecture
	 *      with the Structure Model Index. Comput Methods Biomech Biomed Engin
	 *      1997;1(1):15-23. <a
	 *      href="http://dx.doi.org/10.1080/01495739708936692">doi:
	 *      10.1080/01495739708936692</a>
	 *      </p>
	 * 
	 * @param imp
	 *            binary 3D image
	 * @param voxelResampling
	 *            how much to resample voxels while creatign surface mesh
	 * @param meshSmoothing
	 *            how much smoothign to apply to the mesh
	 * @return SMI
	 */
	@SuppressWarnings("unchecked")
	public static double hildRueg(ImagePlus imp, int voxelResampling,
			float meshSmoothing) {
		int threshold = 128;
		final boolean[] channels = { true, false, false };
		final double r = imp.getCalibration().pixelWidth / 100;
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

		// get all the unique vertices
		ArrayList<Point3f> vertices = new ArrayList<Point3f>();
		final int nPoints = triangles.size();
		for (int p = 0; p < nPoints; p++) {
			IJ.showStatus("Finding unique vertices...");
			IJ.showProgress(p, nPoints);
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

		// associate each unique vertex with the triangles around it
		final int nVertices = vertices.size();
		final int nTriangles = nPoints / 3;
		ArrayList<ArrayList<Integer>> trianglesAtVertex = new ArrayList<ArrayList<Integer>>(
				nVertices);
		for (int i = 0; i < nVertices; i++) {
			IJ.showStatus("Associating vertices with triangles...");
			IJ.showProgress(i, nVertices);
			Point3f vertex = vertices.get(i);
			ArrayList<Integer> tr = new ArrayList<Integer>();
			for (int p = 0; p < nPoints; p++) {
				Point3f testPoint = triangles.get(p);
				if (vertex.equals(testPoint)) {
					tr.add(p);
				}
			}
			trianglesAtVertex.add(tr);
		}

		// get the normals of the triangles around each vertex
		// and calculate the normal of the vertex as the mean triangle normal
		IJ.showStatus("Calculating vertex normals...");
		Point3f[] triangleCentroids = new Point3f[nTriangles];
		Point3f[] triangleNormals = new Point3f[nTriangles];
		Hashtable hash = new Hashtable();
		for (int i = 0; i < nTriangles; i++) {
			triangleCentroids[i] = new Point3f();
			triangleNormals[i] = new Point3f();
		}
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
				Point3f surfaceNormal = VectorProduct.crossProduct(point0,
						point1, point2);
				sumNormals.x += surfaceNormal.x;
				sumNormals.y += surfaceNormal.y;
				sumNormals.z += surfaceNormal.z;
			}

			Point3f vertex = vertices.get(i);

			Point3f normal = new Point3f();
			normal.x = sumNormals.x / vT;
			normal.y = sumNormals.y / vT;
			normal.z = sumNormals.z / vT;

			// Turn normal into a unit vector
			final double length = Trig.distance3D(normal);
			normal.x /= length;
			normal.y /= length;
			normal.z /= length;

			hash.put(vertex, normal);
		}

		// move all the points by the unit normal * small increment r
		ArrayList<Point3f> movedTriangles = new ArrayList<Point3f>();
		for (int t = 0; t < nPoints; t++) {
			IJ.showStatus("Dilating surface mesh...");
			IJ.showProgress(t, nPoints);
			Point3f point = triangles.get(t);
			Point3f newPoint = (Point3f) point.clone();
			Point3f normal = (Point3f) hash.get(point);
			newPoint.x += normal.x * r;
			newPoint.y += normal.y * r;
			newPoint.z += normal.z * r;
			movedTriangles.add(newPoint);
		}

		CustomTriangleMesh surface2 = new CustomTriangleMesh(movedTriangles,
				colour, 0.0f);
		IJ.showStatus("Smoothing surface mesh...");
		MeshEditor.smooth(surface, meshSmoothing);

		double s2 = MeasureSurface.getSurfaceArea(surface2.getMesh());
		double sR = (s2 - s1) / r;
		double smi = 6 * sR * v / (s1 * s1);
		IJ.showStatus("SMI calculated.");
		return smi;
	}
}
