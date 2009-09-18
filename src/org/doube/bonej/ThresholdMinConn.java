package org.doube.bonej;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;

public class ThresholdMinConn implements PlugInFilter {

	private ImagePlus imp;

	public int setup(String arg, ImagePlus imp) {
		this.imp = imp;
		if (imp == null || imp.getNSlices() < 2) {
			IJ.showMessage("A stack must be open");
			return DONE;
		}
		
		return DOES_8G + DOES_16 + STACK_REQUIRED;
	}

	public void run(ImageProcessor ip) {
		int[] histogram = getStackHistogram(this.imp);
		
		
		ip = this.imp.getProcessor();
		ip.getAutoThreshold(histogram);
	}

	/**
	 * Get a 256 bin histogram of stack's pixel values
	 * 
	 * @param imp2
	 * @return
	 */
	private int[] getStackHistogram(ImagePlus imp2) {
		int[] histogram = new int[256];
		int w = imp2.getWidth();
		int h = imp2.getHeight();
		int d = imp2.getStackSize();
		ImageStack stack = imp2.getStack();
		if (imp2.getBitDepth() == 8){
			for (int z = 1; z <= d; z++){
				ImageProcessor sliceIP = stack.getProcessor(z);
				int[] sliceHistogram = sliceIP.getHistogram();
				for (int i = 0; i < 256; i++){
					histogram[i] += sliceHistogram[i];
				}
			}
		}

		else if (imp2.getBitDepth() == 16){
			long[] stackHistogram = new long[(int)Float.MAX_VALUE];
		}
		return histogram;
	}

}
