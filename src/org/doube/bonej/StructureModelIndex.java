package org.doube.bonej;

import java.util.List;

import javax.vecmath.Color3f;
import javax.vecmath.Point3f;

import ij.IJ;
import ij.ImagePlus;
import ij.plugin.PlugIn;

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
		double smi = skyScan(imp);
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
	 * @return
	 */
	public static double skyScan(ImagePlus imp) {
		int threshold = 128;
		final boolean[] channels = { true, false, false };
		int resamplingF = 2;
		MCTriangulator mct = new MCTriangulator();
		List<Point3f> points = mct.getTriangles(imp, threshold, channels,
				resamplingF);
		final Color3f colour = new Color3f(0.0f, 0.0f, 0.0f);
		CustomTriangleMesh surface = new CustomTriangleMesh(points, colour,
				0.0f);
		double v = surface.getVolume();

		double s1 = MeasureSurface.getSurfaceArea(points);

		Dilate d = new Dilate();
		ImagePlus imp2 = d.dilate(imp, 255);

		points = mct.getTriangles(imp2, threshold, channels, resamplingF);
		double s2 = MeasureSurface.getSurfaceArea(points);
		double smi = 6 * ((s2 - s1) * v / (s1 * s1));
		return smi;
	}
}
