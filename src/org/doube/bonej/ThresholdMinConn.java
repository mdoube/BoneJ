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

	public int setup(String arg, ImagePlus imp) {
		this.imp = imp;
		if (imp == null || imp.getNSlices() < 2) {
			IJ.showMessage("A stack must be open");
			return DONE;
		}

		return DOES_8G + DOES_16 + STACK_REQUIRED;
	}

	public void run(ImageProcessor ip) {
		if (!showDialog()) {
			return;
		}
		ImageCheck ic = new ImageCheck();
		if (ic.isBinary(imp)) {
			IJ.error("");
			return;
		}

		int[] histogram = getStackHistogram(imp);
		double threshold = (double) ip.getAutoThreshold(histogram);

		if (!thresholdOnly) {
			int[] testThreshold = getTestThreshold(imp, histogram);
			double[] conns = getConns(imp, testThreshold, subVolume);
			double minimum = getMinimum(testThreshold, conns);
			threshold = checkMinimum(minimum, histogram);
			if (doPlot)
				showPlot(testThreshold, conns);
		}
		IJ.log(imp.getTitle() + " threshold  = " + IJ.d2s(threshold, 1));

		if (applyThreshold) {
			int w = imp.getWidth();
			int h = imp.getHeight();
			ImageStack stack = imp.getImageStack();
			ImageStack stack2 = new ImageStack(w, h);
			int nPixels = w * h;
			for (int z = 1; z <= imp.getStackSize(); z++) {
				byte[] slice = new byte[nPixels];
				ip = stack.getProcessor(z);
				short[] pixels = (short[]) ip.getPixels();
				for (int i = 0; i < nPixels; i++) {
					int value = pixels[i] & 0xffff;
					if (value > threshold) {
						slice[i] = (byte) 255;
					} else {
						slice[i] = (byte) 0;
					}
				}
				stack2.addSlice("", slice);
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
	 * test passes, the original minumum value is returned, otherwise the
	 * autoThreshold of the histogram is returned.
	 * 
	 * @param minimum
	 * @param histogram
	 * @return
	 */
	private double checkMinimum(double minimum, int[] histogram) {
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
	 * 
	 * @param imp2
	 * @param testThreshold
	 * @return
	 */
	private double[] getConns(ImagePlus imp2, int[] testThreshold, int subVolume) {
		int nTests = testThreshold.length;
		double[] conns = new double[nTests];

		// make a stack out of imp2 that is no greater than 256 pixels in
		// any dimension
		ImageStack stack = imp2.getImageStack();

		int width = (int) Math.min(imp2.getWidth(), subVolume);
		int height = (int) Math.min(imp2.getHeight(), subVolume);
		int depth = (int) Math.min(imp2.getStackSize(), subVolume);
		int oldWidth = imp2.getWidth();

		ImageStack stack2 = new ImageStack(width, height);
		for (int z = 1; z <= depth; z++) {
			short[] pixels = (short[]) stack.getPixels(z);
			short[] newPixels = new short[width * height];
			for (int y = 0; y < height; y++) {
				int offset = y * width;
				int oldOffset = y * oldWidth;
				for (int x = 0; x < width; x++) {
					newPixels[offset + x] = pixels[oldOffset + x];
				}
			}
			stack2.addSlice(stack.getSliceLabel(z), newPixels);
		}

		ImageStack stack3 = new ImageStack(stack2.getWidth(), stack2
				.getHeight());
		ImagePlus imp3 = new ImagePlus();
		// for each test
		for (int i = 0; i < nTests; i++) {
			int thresh = testThreshold[i];
			for (int z = 1; z <= stack2.getSize(); z++) {
				short[] pixels = (short[]) stack2.getPixels(z);
				byte[] tpixels = new byte[pixels.length];
				for (int j = 0; j < pixels.length; j++) {
					int value = pixels[j] & 0xffff;
					if (value > thresh) {
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
			imp3.setCalibration(imp2.getCalibration());
			imp3.show();
			if (!imp3.isInvertedLut())
				IJ.run("Invert LUT");
			Purify p = new Purify();
			Erode e = new Erode();
			Dilate d = new Dilate();
			Object[] result = p.purify(imp3, 4, false, false);
			replaceImage(imp3, (ImagePlus) result[1]);
			e.erode(imp3, 255, false).show();
			result = p.purify(imp3, 4, false, false);
			replaceImage(imp3, (ImagePlus) result[1]);
			d.dilate(imp3, 255, false).show();

			// get the connectivity
			Connectivity con = new Connectivity();
			double sumEuler = con.getSumEuler(imp3);
			double deltaChi = con.getDeltaChi(imp3, sumEuler);
			double connectivity = con.getConnectivity(deltaChi);
			// add connectivity to the array
			conns[i] = connectivity;
			// IJ.log("Connectivity is " + conns[i] + " for threshold "
			// + testThreshold[i]);
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
			// scale the histogram so that getAutoThreshold()
			//returns a sensible value
			int mode = 0;
			for (int i = 0; i < histogram.length; i++) {
				if (histogram[i] > mode) {
					mode = histogram[i];
				}
			}
			if (mode > 1e5) {
				double factor = 1e5 / mode;
				for (int i = 0; i < histogram.length; i++) {
					histogram[i] = (int) Math.round(histogram[i] * factor);
				}
				mode = 0;
				for (int i = 0; i < histogram.length; i++) {
					if (histogram[i] > mode) {
						mode = histogram[i];
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
