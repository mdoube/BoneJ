package org.doube.bonej;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.gui.Plot;
import ij.measure.CurveFitter;
import ij.plugin.PlugIn;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;

public class ThresholdMinConn implements PlugIn {

	// private ImagePlus imp;
	private int testCount = 11, subVolume = 256;
	private double testRange = 0.2;

	/** Show a plot of connectivity vs. threshold */
	private boolean doPlot = false;

	/** Apply the threshold to the stack once it is found */
	private boolean applyThreshold = false;

	/**
	 * Return the autothreshold for the stack histogram without doing
	 * connectivity analysis
	 */
	private boolean thresholdOnly = false;

	public void run(String arg) {
		if (!ImageCheck.checkIJVersion())
			return;
		ImagePlus imp = IJ.getImage();
		ImageProcessor ip = imp.getProcessor();
		ImageCheck ic = new ImageCheck();
		if (ic.isBinary(imp)) {
			IJ.error("Can't threshold a binary image");
			return;
		}
		if (!showDialog()) {
			return;
		}

		if (ic.dicomVoxelDepth(imp) != imp.getCalibration().pixelDepth) {
			IJ.run("Properties...");
		}

		int[] histogram = getStackHistogram(imp);
		double threshold = (double) ip.getAutoThreshold(histogram);

		if (!thresholdOnly) {
			int[] testThreshold = getTestThreshold(imp, histogram);
			double[] conns = getConns(imp, testThreshold, subVolume);
			double minimum = getMinimum(testThreshold, conns);
			threshold = checkMinimum(imp, minimum, histogram);
			if (doPlot)
				showPlot(testThreshold, conns);
		}
		IJ.log(imp.getTitle() + " threshold  = " + IJ.d2s(threshold, 1));

		if (applyThreshold) {
			final int w = imp.getWidth();
			final int h = imp.getHeight();
			final int d = imp.getStackSize();
			// final int nPixels = w * h;
			ImageStack stack = imp.getImageStack();
			ImageStack stack2 = new ImageStack(w, h);
			for (int z = 1; z <= d; z++) {
				// byte[] slice = new byte[nPixels];
				ip = stack.getProcessor(z);
				ByteProcessor bp = new ByteProcessor(w, h);
				for (int y = 0; y < h; y++) {
					for (int x = 0; x < w; x++) {
						final double pixel = (double) ip.get(x, y);
						if (pixel > threshold) {
							bp.set(x, y, 255);
						} else {
							bp.set(x, y, 0);
						}
					}
				}
				// short[] pixels = (short[]) ip.getPixels();
				// for (int i = 0; i < nPixels; i++) {
				// int value = pixels[i] & 0xffff;
				// if (value > threshold) {
				// slice[i] = (byte) 255;
				// } else {
				// slice[i] = (byte) 0;
				// }
				// }
				stack2.addSlice(stack.getSliceLabel(z), bp);
			}
			imp.setStack(imp.getTitle(), stack2);
			IJ.selectWindow(imp.getTitle());
			if (!imp.isInvertedLut())
				IJ.run("Invert LUT");
		}
		return;
	}

	/**
	 * Check the calculated optimal (parabolic minimum) value for sanity. If the
	 * test passes, the original minimum value is returned, otherwise the
	 * autoThreshold of the histogram is returned.
	 * 
	 * @param minimum
	 * @param histogram
	 * @return
	 */
	private double checkMinimum(ImagePlus imp, double minimum, int[] histogram) {
		double threshold = minimum;
		ImageProcessor ip = imp.getProcessor();

		// threshold cannot be greater or less than min and max of histogram
		int i = 0;
		while (histogram[i] == 0)
			i++;
		int histogramMin = i;

		i = histogram.length - 1;
		while (histogram[i] == 0)
			i--;
		int histogramMax = i;

		if (minimum < histogramMin || minimum > histogramMax) {
			threshold = ip.getAutoThreshold(histogram);
			IJ.log("Calculated threshold is outside bounds of pixel values. "
					+ "Using histogram-based auto threshold.");
		}

		return threshold;
	}

	/**
	 * 
	 * @param imp2
	 * @param histogram
	 * @return
	 */
	private int[] getTestThreshold(ImagePlus imp2, int[] histogram) {
		ImageProcessor ip = imp2.getProcessor();
		int startThreshold = ip.getAutoThreshold(histogram);

		// get a range of thresholds to test
		int nTests = testCount;
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
		// String parabola = cf.getResultString();
		// IJ.log(parabola);
		double[] params = cf.getParams();
		double b = params[1], c = params[2];
		double xmin = -b / (2 * c);
		// double ymin = a + b * xmin + c * xmin * xmin;
		// IJ.log("minimum connectivity is at (" + xmin + ", " + ymin + ")");
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
	 * Calculate connectivity after threshold-purify-erode-purify-dilate for
	 * several threshold values. 
	 * 
	 * @param imp2
	 * @param testThreshold array of test threshold values (from getTestThreshold)
	 * @return array containing connectivity resulting from each test threshold
	 */
	private double[] getConns(ImagePlus imp2, int[] testThreshold, int subVolume) {
		int nTests = testThreshold.length;
		double[] conns = new double[nTests];

		// make a stack out of imp2 that is no greater than
		// subvolume pixels in any dimension
		ImageStack stack = imp2.getImageStack();

		final int width = (int) Math.min(imp2.getWidth(), subVolume);
		final int height = (int) Math.min(imp2.getHeight(), subVolume);
		final int depth = (int) Math.min(imp2.getStackSize(), subVolume);

		ImageStack stack2 = new ImageStack(width, height);
		for (int z = 1; z <= depth; z++) {
			byte[] nullPixels = new byte[width * height];
			stack2.addSlice(stack.getSliceLabel(z), nullPixels);
			ImageProcessor ip = stack.getProcessor(z);
			ImageProcessor ip2 = stack2.getProcessor(z);
			for (int y = 0; y < height; y++) {
				for (int x = 0; x < width; x++) {
					ip2.set(x, y, ip.get(x, y));
				}
			}
		}
		
		ImagePlus imp3 = new ImagePlus();
		for (int i = 0; i < nTests; i++) {
			//apply threshold
			final int thresh = testThreshold[i];
			ImageStack stack3 = new ImageStack(width, height);
			for (int z = 1; z <= depth; z++) {
				ImageProcessor ip2 = stack2.getProcessor(z);
				ByteProcessor bp = new ByteProcessor(width, height);
				for (int y = 0; y < height; y++) {
					for (int x = 0; x < width; x++) {
						if (ip2.get(x, y) > thresh) {
							bp.set(x, y, 255);
						} else {
							bp.set(x, y, 0);
						}
					}
				}
				if (z > stack3.getSize())
					stack3.addSlice(stack2.getSliceLabel(z), bp);
				else
					stack3.setPixels(bp, z);
			}
			// purify
			imp3.setStack("Threshold " + (i + 1) + "/" + nTests, stack3);
			imp3.setCalibration(imp2.getCalibration());
			imp3.show();
			if (!imp3.isInvertedLut())
				IJ.run("Invert LUT");
			Purify p = new Purify();
			Erode e = new Erode();
			Dilate d = new Dilate();
			Object[] result = p.purify(imp3, 4, false);
			replaceImage(imp3, (ImagePlus) result[1]);
			e.erode(imp3, 255).show();
			result = p.purify(imp3, 4, false);
			replaceImage(imp3, (ImagePlus) result[1]);
			d.dilate(imp3, 255).show();

			// get the connectivity
			Connectivity con = new Connectivity();
			double sumEuler = con.getSumEuler(imp3);
			double deltaChi = con.getDeltaChi(imp3, sumEuler);
			double connectivity = con.getConnectivity(deltaChi);
			// add connectivity to the array
			conns[i] = connectivity;
		}
		imp3.close();
		return conns;
	}

	/**
	 * Replace the image in imp with imp2
	 * 
	 * @param imp
	 * @param imp2
	 */
	private void replaceImage(ImagePlus imp, ImagePlus imp2) {
		ImageStack stack2 = imp2.getStack();
		imp.setStack(null, stack2);
		imp.show();
		if (!imp.isInvertedLut())
			IJ.run("Invert LUT");
	}

	/**
	 * Get a histogram of stack's pixel values
	 * 
	 * @param imp2
	 * @return
	 */
	public int[] getStackHistogram(ImagePlus imp2) {
		final int d = imp2.getStackSize();
		ImageStack stack = imp2.getStack();

		if (stack.getSize() == 1) {
			return stack.getProcessor(1).getHistogram();
		}

		else if (imp2.getBitDepth() == 8) {
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
			int[] histogram = new int[65536];
			for (int z = 1; z <= d; z++) {
				IJ.showStatus("Getting stack histogram...");
				IJ.showProgress(z, d);
				ImageProcessor sliceIP = stack.getProcessor(z);
				int[] sliceHistogram = sliceIP.getHistogram();
				for (int i = 0; i < 65536; i++) {
					histogram[i] += sliceHistogram[i];
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
		gd.addCheckbox("Threshold Only", false);
		gd.addNumericField("Tests", testCount, 0);
		gd.addNumericField("Range", testRange, 2);
		gd.addNumericField("Subvolume Size", subVolume, 0);
		gd.showDialog();
		if (gd.wasCanceled()) {
			return false;
		} else {
			doPlot = gd.getNextBoolean();
			applyThreshold = gd.getNextBoolean();
			thresholdOnly = gd.getNextBoolean();
			testCount = (int) Math.floor(gd.getNextNumber());
			testRange = gd.getNextNumber();
			subVolume = (int) Math.floor(gd.getNextNumber());
			return true;
		}
	}
}
