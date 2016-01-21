package org.doube.bonej;

import java.awt.image.ColorModel;

import org.doube.util.ImageCheck;
import org.doube.util.UsageReporter;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
import ij.process.ByteProcessor;

/**
 * This class implements the dilation filter. The kernel size is fixed with a
 * diameter of 3 pixels. This makes sense the operation gets rapidly
 * computationally more expensive with increasing diameter. Computational
 * complexity is related to the third power of the diameter, so 2-fold diameter
 * means 8-fold computation time. For complexity reasons, the implementation
 * uses a 6-neighbour- hood and not a 27-neighborhood.
 *
 * Imported from Fiji's VIB_.jar on 2009-09-21
 *
 * @author Benjamin Schimd
 *
 */
public class Dilate implements PlugIn {

	private int w, h, d;
	private byte[][] pixels_in;
	private byte[][] pixels_out;

	public void run(final String arg) {
		if (!ImageCheck.checkEnvironment())
			return;
		final ImagePlus imp = IJ.getImage();
		if (null == imp) {
			IJ.noImage();
			return;
		}
		final GenericDialog gd = new GenericDialog("Dilate");
		gd.addNumericField("Iso value", 255, 0);
		gd.addHelp("http://pacific.mpi-cbg.de/wiki/index.php/3D_Binary_Filters");
		gd.showDialog();
		if (gd.wasCanceled())
			return;
		final ImagePlus imp2 = dilate(imp, (int) gd.getNextNumber());
		imp.setStack(null, imp2.getImageStack());
		UsageReporter.reportEvent(this).send();
		return;
	}

	public ImagePlus dilate(final ImagePlus image, final int threshold) {

		// Determine dimensions of the image
		w = image.getWidth();
		h = image.getHeight();
		d = image.getStackSize();

		this.pixels_in = new byte[d][];
		this.pixels_out = new byte[d][];
		for (int z = 0; z < d; z++) {
			this.pixels_in[z] = (byte[]) image.getStack().getPixels(z + 1);
			this.pixels_out[z] = new byte[w * h];
		}

		// iterate
		for (int z = 0; z < d; z++) {
			IJ.showProgress(z, d - 1);
			for (int y = 0; y < h; y++) {
				for (int x = 0; x < w; x++) {
					if (get(x, y, z) == threshold || get(x - 1, y, z) == threshold || get(x + 1, y, z) == threshold
							|| get(x, y - 1, z) == threshold || get(x, y + 1, z) == threshold
							|| get(x, y, z - 1) == threshold || get(x, y, z + 1) == threshold)

						set(x, y, z, threshold);
					else
						set(x, y, z, get(x, y, z));
				}
			}
		}

		final ColorModel cm = image.getStack().getColorModel();

		// create output image
		final ImageStack stack = new ImageStack(w, h);
		for (int z = 0; z < d; z++) {
			stack.addSlice(image.getImageStack().getSliceLabel(z + 1), new ByteProcessor(w, h, this.pixels_out[z], cm));
		}
		final ImagePlus imp = new ImagePlus();
		imp.setCalibration(image.getCalibration());
		imp.setStack(null, stack);
		return imp;
	}

	public int get(int x, int y, int z) {
		x = x < 0 ? 0 : x;
		x = x >= w ? w - 1 : x;
		y = y < 0 ? 0 : y;
		y = y >= h ? h - 1 : y;
		z = z < 0 ? 0 : z;
		z = z >= d ? d - 1 : z;
		return this.pixels_in[z][y * w + x] & 0xff;
	}

	public void set(final int x, final int y, final int z, final int v) {
		this.pixels_out[z][y * w + x] = (byte) v;
		return;
	}
}