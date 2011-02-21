package org.doube.util;

import ij.ImagePlus;
import ij.gui.Roi;
import ij.measure.Calibration;
import ij.plugin.frame.RoiManager;

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
	 * Return a list of ROIs that are active in the given slice, s.
	 * ROIs without a slice number are assumed to be active in all slices.
	 * 
	 * @param roiMan
	 * @param s
	 * @return
	 */
	public static ArrayList<Roi> getSliceRoi(RoiManager roiMan, int s) {
		ArrayList<Roi> roiList = new ArrayList<Roi>();
		Roi[] rois = roiMan.getRoisAsArray();
		for (Roi roi : rois){
			int sliceNumber = roiMan.getSliceNumber(roi.getName());
			if (sliceNumber == s || sliceNumber == -1)
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
	 *         If any of the ROIs contains no slice information, z min is set to 1
	 *         and z max is set to Integer.MAX_VALUE
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
				noZroi = true; //found a ROI with no Z info
		}
		if (noZroi) {
			int[] limits = { xmin, xmax, ymin, ymax, 1, Integer.MAX_VALUE };
			return limits;
		} else {
			int[] limits = { xmin, xmax, ymin, ymax, zmin, zmax };
			return limits;
		}
	}
}
