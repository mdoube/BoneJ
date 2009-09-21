package org.doube.bonej;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.gui.Plot;
import ij.measure.CurveFitter;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;

public class ThresholdMinConn implements PlugInFilter {

	private ImagePlus imp;
	private boolean doPlot = false, applyThreshold = false;

	public int setup(String arg, ImagePlus imp) {
		this.imp = imp;
		if (imp == null || imp.getNSlices() < 2) {
			IJ.showMessage("A stack must be open");
			return DONE;
		}

		return DOES_8G + DOES_16 + STACK_REQUIRED;
	}

	public void run(ImageProcessor ip) {
		if (!showDialog()){
			return;
		}

		int[] testThreshold = getTestThreshold(this.imp);
		double[] conns = getConns(this.imp, testThreshold);
		double minimum = getMinimum(testThreshold, conns);

		if (doPlot)
			showPlot(testThreshold, conns);

		if (applyThreshold) {
			ip.setThreshold(minimum, getMaximum(this.imp), ImageProcessor.BLACK_AND_WHITE_LUT);
			IJ.run("Convert to Mask", " ");
		}
	}

	private double getMaximum(ImagePlus imp2) {
		double max = Double.MIN_VALUE;
		ImageStack stack2 = imp2.getImageStack();
		for (int z = 1; z <= imp2.getStackSize(); z++){
			short[] pixels = (short[])stack2.getPixels(z);
			for (int i = 0; i < pixels.length; i++){
				max = Math.max(max, pixels[i]);
			}
		}
		return max;
	}

	private int[] getTestThreshold(ImagePlus imp2) {
		int[] histogram = getStackHistogram(imp2);
		ImageProcessor ip = imp2.getProcessor();
		int startThreshold = ip.getAutoThreshold(histogram);

		// get a range of thresholds to test
		int nTests = 21;
		double testRange = 0.2; // as a fraction of the autothreshold
		double testStep = 2 * testRange * startThreshold / (nTests - 1);
		int[] testThreshold = new int[nTests];
		for (int i = 0; i < nTests; i++) {
			testThreshold[i] = (int) Math.round(startThreshold
					* (1 - testRange) + i * testStep);
		}

		return testThreshold;
	}

	/**
	 * Fit a parabola to the threshold and connectivity data and return its
	 * minimum
	 * 
	 * @param testThreshold
	 * @param conns
	 * @return
	 */
	private double getMinimum(int[] testThreshold, double[] conns) {
		// convert xData to double
		int nPoints = testThreshold.length;
		double[] xData = new double[nPoints];
		for (int i = 0; i < nPoints; i++) {
			xData[i] = (double) testThreshold[i];
		}
		CurveFitter cf = new CurveFitter(xData, conns);
		cf.doFit(CurveFitter.POLY2);
		String parabola = cf.getResultString();
		IJ.log(parabola);
		double[] params = cf.getParams();
		double a = params[0], b = params[1], c = params[2];
		double xmin = -b / (2 * c);
		double ymin = a + b * xmin + c * xmin * xmin;
		IJ.log("minimum connectivity is at (" + xmin + ", " + ymin + ")");
		return xmin;
	}

	/**
	 * Display a graph showing connectivity vs. threshold
	 * 
	 * @param testThreshold
	 * @param conns
	 */
	private void showPlot(int[] testThreshold, double[] conns) {
		// convert arrays to floats
		int nPoints = testThreshold.length;
		float[] xData = new float[nPoints];
		float[] yData = new float[nPoints];
		double xMin = Double.MAX_VALUE;
		double xMax = Double.MIN_VALUE;
		double yMax = Double.MIN_VALUE;
		for (int i = 0; i < nPoints; i++) {
			xData[i] = (float) testThreshold[i];
			yData[i] = (float) conns[i];
			xMin = Math.min(xMin, xData[i]);
			xMax = Math.max(xMax, xData[i]);
			yMax = Math.max(yMax, yData[i]);
		}
		Plot plot = new Plot("Connectivity vs. Threshold", "Threshold",
				"Connectivity", xData, yData);
		plot.addPoints(xData, yData, Plot.CIRCLE);
		plot.setLimits(xMin, xMax, 0, yMax);
		plot.draw();

		ImageProcessor plotIp = plot.getProcessor();
		ImagePlus plotImage = new ImagePlus();
		plotImage.setProcessor("Connectivity vs. Threshold", plotIp);
		plotImage.show();
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
		ImageStack stack3 = new ImageStack(stack2.getWidth(), stack2
				.getHeight());
		ImagePlus imp3 = new ImagePlus();
		// for each test
		for (int i = 0; i < nTests; i++) {
			int thresh = testThreshold[i];
			for (int z = 1; z < imp2.getNSlices(); z++) {
				short[] pixels = (short[]) stack2.getPixels(z);
				byte[] tpixels = new byte[pixels.length];
				for (int j = 0; j < pixels.length; j++) {
					if (pixels[j] > thresh) {
						tpixels[j] = (byte) 255;
					} else {
						tpixels[j] = (byte) 0;
					}
				}
				if (z > stack3.getSize())
					stack3.addSlice("" + z, tpixels);
				else
					stack3.setPixels(tpixels, z);
			}
			// purify
			imp3.setStack("Threshold " + (i + 1) + "/" + nTests, stack3);
			imp3.show();
			if (!imp3.isInvertedLut())
				IJ.run("Invert LUT");
			IJ.run("Purify", "chunk=4");
			IJ.run("Erode (3D)", "iso=255");
			IJ.run("Purify", "chunk=4");
			IJ.run("Dilate (3D)", "iso=255");

			// get the connectivity
			Connectivity con = new Connectivity();
			double sumEuler = con.getSumEuler(imp3);
			double deltaChi = con.getDeltaChi(imp3, sumEuler);
			double connectivity = con.getConnectivity(deltaChi);
			// add connectivity to the array
			conns[i] = connectivity;
			IJ.log("Connectivity is " + conns[i] + " for threshold "
					+ testThreshold[i]);
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

	private boolean showDialog() {
		GenericDialog gd = new GenericDialog("Options");
		gd.addCheckbox("Show Plot", true);
		gd.addCheckbox("Apply Threshold", false);
		gd.showDialog();
		if (gd.wasCanceled()) {
			return false;
		} else {
			doPlot = gd.getNextBoolean();
			applyThreshold = gd.getNextBoolean();
			return true;
		}
	}
}
