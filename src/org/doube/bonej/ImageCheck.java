package org.doube.bonej;

import ij.ImagePlus;
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
		if (imp == null)
			return false;
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
		if (imp == null)
			return false;
		if (imp.getStackSize() < 2)
			return false;
		return true;
	}
}
