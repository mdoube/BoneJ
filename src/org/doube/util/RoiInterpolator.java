package org.doube.util;

import java.awt.Rectangle;
import java.util.ArrayList;

import org.doube.geometry.BinaryInterpolator;

import ij.plugin.PlugIn;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;
import ij.gui.Roi;
import ij.plugin.filter.ThresholdToSelection;
import ij.plugin.frame.RoiManager;

public class RoiInterpolator implements PlugIn {
	int[][] idt;
	int w, h;

	public void run(String arg) {
		RoiManager roiman = RoiManager.getInstance();
		if (roiman == null){
			IJ.error("Please populate the ROI Manager with ROIs");
			return;
		}
		Roi[] rois = roiman.getRoisAsArray();
		int xmax = 0;
		int xmin = Integer.MAX_VALUE;
		int ymax = 0;
		int ymin = Integer.MAX_VALUE;
		int zmax = 1;
		int zmin = Integer.MAX_VALUE;
		ArrayList<Integer> templateSlices = new ArrayList<Integer>();
		for (Roi roi : rois){
			final int slice = roiman.getSliceNumber(roi.getName());
			templateSlices.add(new Integer(slice));
			if (slice == 0) //ignore non-slice associated ROIs
				continue;
			zmin = Math.min(slice, zmin);
			zmax = Math.max(slice, zmax);
			Rectangle bounds = roi.getBounds();
			xmin = Math.min(xmin, bounds.x);
			ymin = Math.min(ymin, bounds.y);
			xmax = Math.max(xmax, bounds.x + bounds.width);
			ymax = Math.max(ymax, bounds.y + bounds.height);
		}
		//create the binary stack
		final int nSlices = zmax - zmin + 1;
		ImageStack stack = new ImageStack(xmax, ymax);
		for (int s = 0; s < nSlices; s++){
			ByteProcessor bp = new ByteProcessor(xmax, ymax);
			bp.setColor(255);
			for (Roi roi : rois){
				final int slice = roiman.getSliceNumber(roi.getName());
				if (slice == zmin + s){
					bp.setRoi(roi);
					if (roi.getType() == Roi.RECTANGLE)
						bp.fill();
					else 
						bp.fill(roi);
				}
			}
			stack.addSlice(""+s, bp);
		}
		//do the binary interpolation
		BinaryInterpolator bi = new BinaryInterpolator();
		bi.run(stack);
		ImagePlus binary = new ImagePlus("interpolated", stack);

		//get the ROIs
		ThresholdToSelection ts = new ThresholdToSelection();
		ts.setup("", binary);
		for(int s = 0; s < nSlices; s++) {
			if (templateSlices.contains(new Integer(s + zmin)))
				continue;
			ImageProcessor bp = stack.getProcessor(s+1);
			int threshold = 255;
			bp.setThreshold(threshold, threshold, ImageProcessor.NO_LUT_UPDATE);
			Roi roi = ts.convert(bp);
			roi.setPosition(s + zmin);
			roiman.addRoi(roi);
		}
		IJ.showStatus("ROIs interpolated");
	}
}
