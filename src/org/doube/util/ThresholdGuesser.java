package org.doube.util;

import ij.ImagePlus;
import ij.measure.Calibration;

import org.doube.bonej.ThresholdMinConn;

public class ThresholdGuesser {

	/**
	 * Set up default thresholds and report them to the user, using auto
	 * thresholding of a histogram generated from all pixels in a stack
	 * 
	 * @return double[2] containing minimum and maximum thresholds
	 */
	public static double[] setDefaultThreshold(ImagePlus imp) {
		Calibration cal = imp.getCalibration();
		double min = 0;
		double max = 0;
		// set some sensible thresholding defaults
		ThresholdMinConn tmc = new ThresholdMinConn();
		int[] histogram = tmc.getStackHistogram(imp);
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
		double[] thresholds = { min, max };
		return thresholds;
	}
}
