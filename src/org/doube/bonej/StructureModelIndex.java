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

		double smi = skyScan(imp, voxelResampling, meshSmoothing);
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
}
