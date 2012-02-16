package org.doube.util;

import ij.ImagePlus;
import ij.ImageStack;

public class StackStats {
	/**
	 * Work out some summary stats
	 * 
	 * @param imp
	 *            32-bit thickness image
	 * @return double[] containing mean, standard deviation and maximum as its
	 *         0th and 1st and 2nd elements respectively
	 * 
	 */
	public static double[] meanStdDev(ImagePlus imp) {
		final int w = imp.getWidth();
		final int h = imp.getHeight();
		final int d = imp.getStackSize();
		final int wh = w * h;
		final ImageStack stack = imp.getStack();
		long pixCount = 0;
		double sumThick = 0;
		double maxThick = 0;

		for (int s = 1; s <= d; s++) {
			final float[] slicePixels = (float[]) stack.getPixels(s);
			for (int p = 0; p < wh; p++) {
				final double pixVal = slicePixels[p];
				if (pixVal > 0) {
					sumThick += pixVal;
					maxThick = Math.max(maxThick, pixVal);
					pixCount++;
				}
			}
		}
		final double meanThick = sumThick / pixCount;

		double sumSquares = 0;
		for (int s = 1; s <= d; s++) {
			final float[] slicePixels = (float[]) stack.getPixels(s);
			for (int p = 0; p < wh; p++) {
				final double pixVal = slicePixels[p];
				if (pixVal > 0) {
					final double residual = meanThick - pixVal;
					sumSquares += residual * residual;
				}
			}
		}
		final double stDev = Math.sqrt(sumSquares / pixCount);
		double[] stats = { meanThick, stDev, maxThick };
		return stats;
	}

}
