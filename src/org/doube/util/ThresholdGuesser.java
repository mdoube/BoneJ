package org.doube.util;

import ij.ImagePlus;
import ij.measure.Calibration;

public class ThresholdGuesser {
	private static final double airHU = -1000;

	/**
	 * Set up default thresholds and report them to the user as HU values if the
	 * image has HU calibration or plain values if not. Used as a first guess
	 * for dialogs that have to handle both HU and uncalibrated images.
	 *
	 * @return double[2] containing minimum and maximum thresholds
	 */
	public static double[] setDefaultThreshold(final ImagePlus imp) {
		final Calibration cal = imp.getCalibration();
		double min = 0;
		double max = 0;
		if (!ImageCheck.huCalibrated(imp)) {
			// set some sensible thresholding defaults
			final int[] histogram = StackStats.getStackHistogram(imp);
			final int histoLength = histogram.length;
			int histoMax = histoLength - 1;
			for (int i = histoLength - 1; i >= 0; i--) {
				if (histogram[i] > 0) {
					histoMax = i;
					break;
				}
			}
			min = imp.getProcessor().getAutoThreshold(histogram);
			max = histoMax;
			if (cal.isSigned16Bit() && cal.getCValue(0) == 0) {
				min += Short.MIN_VALUE;
				max += Short.MIN_VALUE;
			}
		} else {
			// default bone thresholds are 0 and 4000 HU
			min = airHU + 1000;
			max = airHU + 5000;
		}
		final double[] thresholds = { min, max };
		return thresholds;
	}
}
