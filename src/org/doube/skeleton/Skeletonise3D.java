package org.doube.skeleton;

/**
 * Skeletonise3D plugin for ImageJ(C).
 * Copyright (C) 2008 Ignacio Arganda-Carreras 
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation (http://www.gnu.org/licenses/gpl.txt )
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 * 
 */

import org.doube.util.ImageCheck;
import org.doube.util.UsageReporter;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.plugin.PlugIn;
import ij.plugin.filter.Skeletonize3D;

/**
 * Main class. This class is a plugin for the ImageJ interface for 2D and 3D
 * thinning (skeletonization) of binary images (2D/3D).
 * 
 * <p>
 * This work is an implementation by Ignacio Arganda-Carreras of the 3D thinning
 * algorithm from Lee et al., which is included in the core of ImageJ since v1.46j.
 * </p>
 * 
 * @see <p>
 *      Lee, Kashyap and Chu (1994) Building skeleton models via 3-D medial
 *      surface/axis thinning algorithms. Computer Vision, Graphics, and Image
 *      Processing, 56(6):462–478 <a
 *      href="http://dx.doi.org/10.1006/cgip.1994.1042"
 *      >doi:10.1006/cgip.1994.1042</a>.
 *      </p>
 *      <p>
 *      Based on the ITK version from Hanno Homann <a
 *      href="http://hdl.handle.net/1926/1292">
 *      http://hdl.handle.net/1926/1292</a>
 *      </p>
 *      <p>
 *      <a href=
 *      "http://imagejdocu.tudor.lu/doku.php?id=plugin:morphology:skeletonize3d:start"
 *      >Skeletonise3D homepage</a>
 *      </p>
 * 
 * @version 1.0 11/19/2008
 * @author Ignacio Arganda-Carreras <ignacio.arganda@uam.es>
 * 
 */
public class Skeletonise3D implements PlugIn {

	public void run(String run) {
		if (!ImageCheck.checkEnvironment())
			return;
		ImagePlus imp = IJ.getImage();
		if (null == imp) {
			IJ.noImage();
			return;
		}
		ImageCheck ic = new ImageCheck();
		if (!ic.isBinary(imp)) {
			IJ.error("Skeletonise 3D requires a binary image");
			return;
		}

		ImagePlus skeleton = getSkeleton(imp);

		skeleton.show();
		if (imp.isInvertedLut() && !skeleton.isInvertedLut())
			IJ.run("Invert LUT");
		UsageReporter.reportEvent(this).send();
		return;
	}

	/**
	 * Gets a medial axis skeleton from a binary imp using a topology-preserving
	 * iterative algorithm
	 * 
	 * @param imp
	 *            input image
	 * @return skeletonised image
	 */
	public ImagePlus getSkeleton(ImagePlus imp) {
		ImageStack inputImage = imp.getStack();

		// Prepare data
		ImageStack outputImage = prepareData(inputImage);
		ImagePlus imp2 = new ImagePlus("Skeleton of " + imp.getTitle(),
				outputImage);

		Skeletonize3D sk = new Skeletonize3D();
		sk.setup("", imp2);
		sk.run(imp2.getProcessor());

		imp2.setCalibration(imp.getCalibration());

		return imp2;
	}

	/* ----------------------------------------------------------------------- */
	/**
	 * Prepare data for computation. Copy the input image to the output image,
	 * changing from the input type to the output type.
	 * 
	 * @param outputImage
	 *            output image stack
	 */
	private ImageStack prepareData(ImageStack inputImage) {
		final int width = inputImage.getWidth();
		final int height = inputImage.getHeight();
		final int depth = inputImage.getSize();
		// IJ.write("Prepare Data: Copy input to output");
		IJ.showStatus("Prepare Data: Copy input to output ...");

		ImageStack outputImage = new ImageStack(width, height, depth);
		// Copy the input to the output, changing all foreground pixels to
		// have value 1 in the process.
		for (int z = 0; z < depth; z++) {
			byte[] pixels = (byte[]) inputImage.getPixels(z + 1);
			outputImage.setPixels(pixels.clone(), z + 1);
			outputImage.setSliceLabel(inputImage.getSliceLabel(z + 1), z + 1);
			for (int x = 0; x < width; x++)
				for (int y = 0; y < height; y++)
					if (((byte[]) inputImage.getPixels(z + 1))[x + y * width] != 0) {
						((byte[]) outputImage.getPixels(z + 1))[x + y * width] = 1;
					} else {
						((byte[]) outputImage.getPixels(z + 1))[x + y * width] = 0;
					}
		}
		// IJ.write("Prepare Data End");
		IJ.showStatus("Prepare Data End.");
		return outputImage;
	} /* end prepareData */

} /* end Skeletonize3D_ */
