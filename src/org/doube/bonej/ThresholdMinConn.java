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
		int startThreshold = ip.getAutoThreshold(histogram);

		// get a range of thresholds to test
		int nTests = 11;
		double testRange = 0.1; // as a fraction of the autothreshold
		double testStep = 2 * testRange * startThreshold / (nTests - 1);
		int[] testThreshold = new int[nTests];
		for (int i = 0; i < nTests; i++) {
			testThreshold[i] = (int) Math.round(startThreshold
					* (1 - testRange) + i * testStep);
		}

		// for each test threshold, Purify, get connectivity
		double[] conns = getConns(this.imp, testThreshold);

	}

	/**
	 * 
	 * @param imp2
	 * @param testThreshold
	 * @return
	 */
	private double[] getConns(ImagePlus imp2, int[] testThreshold) {
		int nTests = testThreshold.length;
		double[] conns = new double[nTests];
		
		ImageStack stack2 = imp2.getImageStack();
		ImageStack stack3 = new ImageStack(stack2.getWidth(), stack2.getHeight());
		ImagePlus imp3 = new ImagePlus();
		// for each test
		for (int i = 0; i < nTests; i++) {
			int thresh = testThreshold[i];
			for (int z = 1; z < imp2.getNSlices(); z++){
				short[] pixels = (short[])stack2.getPixels(z);
				byte[] tpixels = new byte[pixels.length];
				for (int j = 0; j < pixels.length; j++){
					if (pixels[j] > thresh){
						tpixels[j] = (byte) 255;
					} else {
						tpixels[j] = (byte) 0;
					}
				}
				if (z > stack3.getSize()) stack3.addSlice(""+z, tpixels);
				else stack3.setPixels(tpixels, z);
			}
			// purify
			imp3.setStack("Threshold "+(i+1)+"/"+nTests, stack3);
			imp3.show();
			if (!imp3.isInvertedLut()) IJ.run("Invert LUT");
			IJ.run("Purify", "chunk=4");

			// get the connectivity			
			Connectivity con = new Connectivity();
			double sumEuler = con.getSumEuler(imp3);
			double deltaChi = con.getDeltaChi(imp3, sumEuler);
			double connectivity = con.getConnectivity(deltaChi);
			// add connectivity to the array
			conns[i] = connectivity;
			IJ.log("Connectivity is "+conns[i]+" for threshold "+testThreshold[i]);
		}

		return conns;
	}

	/**
	 * Get a histogram of stack's pixel values
	 * 
	 * @param imp2
	 * @return
	 */
	private int[] getStackHistogram(ImagePlus imp2) {

		int w = imp2.getWidth();
		int h = imp2.getHeight();
		int d = imp2.getStackSize();
		ImageStack stack = imp2.getStack();
		if (imp2.getBitDepth() == 8) {
			int[] histogram = new int[256];
			for (int z = 1; z <= d; z++) {
				IJ.showStatus("Getting stack histogram...");
				IJ.showProgress(z, d);
				ImageProcessor sliceIP = stack.getProcessor(z);
				int[] sliceHistogram = sliceIP.getHistogram();
				for (int i = 0; i < 256; i++) {
					histogram[i] += sliceHistogram[i];
				}
			}
			return histogram;
		}

		else if (imp2.getBitDepth() == 16) {
			short adjustment = 0;
			if (imp2.getCalibration().isSigned16Bit())
				adjustment = Short.MIN_VALUE;
			int[] histogram = new int[Short.MAX_VALUE - Short.MIN_VALUE + 1];
			for (int z = 1; z <= d; z++) {
				IJ.showStatus("Getting stack histogram...");
				IJ.showProgress(z, d);
				ImageProcessor sliceIP = stack.getProcessor(z);
				for (int y = 0; y < h; y++) {
					for (int x = 0; x < w; x++) {
						int i = sliceIP.get(x, y) - adjustment;
						histogram[i]++;
					}
				}
			}
			return histogram;
		} else
			return null;
	}

}
