package org.doube.util;

import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Roi;
import ij.measure.Calibration;
import ij.plugin.frame.RoiManager;
import ij.process.ImageProcessor;

import java.awt.List;
import java.awt.Rectangle;
import java.util.ArrayList;

/**
 * Do useful things with ImageJ's ROI Manager
 * 
 * @author Michael Doube
 * */
public class RoiMan {
	/**
	 * Get the calibrated 3D coordinates of point ROIs from the ROI manager
	 * 
	 * @param imp
	 * @param roiMan
	 * @return double[n][3] containing n (x, y, z) coordinates
	 */
	public static double[][] getRoiManPoints(ImagePlus imp, RoiManager roiMan) {
		Calibration cal = imp.getCalibration();
		double vW = cal.pixelWidth;
		double vH = cal.pixelHeight;
		double vD = cal.pixelDepth;
		int nPoints = 0;
		List listRoi = roiMan.getList();
		Roi[] roiList = roiMan.getRoisAsArray();
		for (int i = 0; i < roiMan.getCount(); i++) {
			Roi roi = roiList[i];
			if (roi.getType() == 10) {
				nPoints++;
			}
		}
		double[][] dataPoints = new double[nPoints][3];
		int j = 0;
		for (int i = 0; i < roiMan.getCount(); i++) {
			Roi roi = roiList[i];
			if (roi.getType() == 10) {
				String label = listRoi.getItem(i);
				Rectangle xy = roi.getBounds();
				dataPoints[j][0] = xy.getX() * vW;
				dataPoints[j][1] = xy.getY() * vH;
				dataPoints[j][2] = roiMan.getSliceNumber(label) * vD;
				j++;
			}
		}
		return dataPoints;
	}

	/**
	 * Return a list of ROIs that are active in the given slice, s. ROIs without
	 * a slice number are assumed to be active in all slices.
	 * 
	 * @param roiMan
	 * @param s
	 * @return
	 */
	public static ArrayList<Roi> getSliceRoi(RoiManager roiMan, int s) {
		ArrayList<Roi> roiList = new ArrayList<Roi>();
		Roi[] rois = roiMan.getRoisAsArray();
		for (Roi roi : rois) {
			int sliceNumber = roiMan.getSliceNumber(roi.getName());
			if (sliceNumber == -1)
				sliceNumber = roi.getPosition();
			if (sliceNumber == s || sliceNumber == 0)
				roiList.add(roi);
		}
		return roiList;
	}

	/**
	 * Find the x, y and z limits of the ROIs in the ROI Manager
	 * 
	 * @param roiMan
	 * @return int[] containing x min, x max, y min, y max, z min and z max, or
	 *         null if there is no ROI Manager or if the ROI Manager is empty.
	 *         If any of the ROIs contains no slice information, z min is set to
	 *         1 and z max is set to Integer.MAX_VALUE
	 */
	public static int[] getLimits(RoiManager roiMan) {
		if (roiMan == null || roiMan.getCount() == 0)
			return null;
		int xmin = Integer.MAX_VALUE;
		int xmax = 0;
		int ymin = Integer.MAX_VALUE;
		int ymax = 0;
		int zmin = Integer.MAX_VALUE;
		int zmax = 1;
		boolean noZroi = false;
		Roi[] rois = roiMan.getRoisAsArray();
		for (Roi roi : rois) {
			Rectangle r = roi.getBounds();
			xmin = Math.min(r.x, xmin);
			xmax = Math.max(r.x + r.width, xmax);
			ymin = Math.min(r.y, ymin);
			ymax = Math.max(r.y + r.height, ymax);
			int slice = roiMan.getSliceNumber(roi.getName());
			if (slice > 0) {
				zmin = Math.min(slice, zmin);
				zmax = Math.max(slice, zmax);
			} else
				noZroi = true; // found a ROI with no Z info
		}
		if (noZroi) {
			int[] limits = { xmin, xmax, ymin, ymax, 1, Integer.MAX_VALUE };
			return limits;
		} else {
			int[] limits = { xmin, xmax, ymin, ymax, zmin, zmax };
			return limits;
		}
	}

	/**
	 * Crop a stack to the limits of the ROIs in the ROI Manager and optionally
	 * fill the background with a single pixel value.
	 * 
	 * @param roiMan
	 *            ROI Manager containing ROIs
	 * @param stack
	 *            input stack
	 * @param fillBackground
	 *            if true, background will be set to value
	 * @param fillValue
	 *            value to set background to
	 * @param padding
	 *            empty pixels to pad faces of cropped stack with
	 * @return cropped copy of input stack
	 */
	public static ImageStack cropStack(RoiManager roiMan, ImageStack stack,
			boolean fillBackground, int fillValue, int padding) {
		int[] limits = getLimits(roiMan);
		final int xmin = limits[0];
		final int xmax = limits[1];
		final int ymin = limits[2];
		final int ymax = limits[3];
		final int zmin = limits[4];
		final int zmax = (limits[5] == Integer.MAX_VALUE) ? stack.getSize()
				: limits[5] + 1 + padding;
		final int w = xmax - xmin;
		final int h = ymax - ymin;

		ImageStack out = new ImageStack(w + 2 * padding, h + 2 * padding);
		for (int z = zmin - 1 - padding; z <= zmax; z++) {
			ImageProcessor ip = stack.getProcessor(z);
			ImageProcessor ipOut = ip.duplicate();
			final int length = ipOut.getPixelCount();
			if (fillBackground)
				for (int i = 0; i < length; i++)
					ipOut.set(i, fillValue);
			ArrayList<Roi> rois = getSliceRoi(roiMan, z);
			for (Roi roi : rois) {
				ip.setRoi(roi);
				Rectangle r = roi.getBounds();
				ImageProcessor mask = ip.getMask();
				final int rh = r.y + r.height;
				final int rw = r.x + r.width;
				for (int y = r.y; y < rh; y++) {
					for (int x = r.x; x < rw; x++) {
						if (mask.get(x - r.x, y - r.y) > 0)
							ipOut.set(x, y, ip.get(x, y));
					}
				}
			}
			ipOut.setRoi(xmin - padding, ymin - padding, w + 2 * padding, h + 2
					* padding);
			ipOut = ipOut.crop();
			out.addSlice(stack.getSliceLabel(z), ipOut);
		}
		return out;
	}

	/**
	 * Remove all ROIs from the ROI manager
	 */
	public static void deleteAll(RoiManager roiMan) {
		Roi[] rois = roiMan.getRoisAsArray();
		for (int i = 0; i < rois.length; i++) {
			if (roiMan.getCount() == 0)
				break;
			roiMan.select(i);
			roiMan.runCommand("delete");
		}
	}
}
