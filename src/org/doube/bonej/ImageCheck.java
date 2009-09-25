package org.doube.bonej;

import ij.IJ;
import ij.ImagePlus;
import ij.measure.Calibration;
import ij.process.ImageStatistics;

/**
 * Check if an image conforms to the type defined by each method. Methods return
 * true if image satisfies conditions and false if it doesn't.
 * 
 * @author mdoube
 * 
 */
public class ImageCheck {

	/**
	 * Check if image is binary
	 * 
	 * @param imp
	 * @return true if image is binary
	 */
	public boolean isBinary(ImagePlus imp) {
		if (imp == null) {
			IJ.noImage();
			return false;
		}
		if (imp.getType() != ImagePlus.GRAY8)
			return false;

		ImageStatistics stats = imp.getStatistics();
		if (stats.histogram[0] + stats.histogram[255] != stats.pixelCount)
			return false;
		return true;
	}

	/**
	 * Check if an image is a multi-slice image stack
	 * 
	 * @param imp
	 * @return true if the image has >= 2 slices
	 */
	public boolean isMultiSlice(ImagePlus imp) {
		if (imp == null) {
			IJ.noImage();
			return false;
		}

		if (imp.getStackSize() < 2)
			return false;
		return true;
	}

	/**
	 * Check if the image's voxels are isotropic in all 3 dimensions (i.e. are
	 * placed on a cubic grid)
	 * 
	 * @param imp
	 * @return true if voxel width == height == depth
	 */
	public boolean isVoxelIsotropic(ImagePlus imp) {
		if (imp == null){
			IJ.noImage();
			return false;
		}
		Calibration cal = imp.getCalibration();
		double vW = cal.pixelWidth;
		double vH = cal.pixelHeight;
		double vD = cal.pixelDepth;

		if (vW != vH)
			return false;
		if (vW != vD)
			return false;
		if (vH != vD)
			return false;

		return true;
	}
}
