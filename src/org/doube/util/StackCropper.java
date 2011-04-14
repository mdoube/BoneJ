package org.doube.util;

import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
import ij.plugin.frame.RoiManager;

/**
 * Implements the stack cropping method from RoiMan as a plugin
 * 
 * @author Michael Doube
 * 
 */
public class StackCropper implements PlugIn {

	public void run(String arg) {
		ImagePlus imp = WindowManager.getCurrentImage();
		RoiManager roiMan = RoiManager.getInstance();
		if (imp == null || roiMan == null)
			return;

		GenericDialog gd = new GenericDialog("Crop Stack by ROI");
		gd.addCheckbox("Replace Original", false);
		gd.addCheckbox("Fill outside", false);
		gd.addNumericField("Fill_value", 0, 0, 6, "");
		gd.addNumericField("Padding", 0, 0);
		gd.showDialog();
		if (gd.wasCanceled())
			return;
		boolean doReplace = gd.getNextBoolean();
		boolean doFill = gd.getNextBoolean();
		int fillValue = (int) gd.getNextNumber();
		int padding = (int) gd.getNextNumber();
		ImageStack stack = RoiMan.cropStack(roiMan, imp.getImageStack(),
				doFill, fillValue, padding);
		if (doReplace) {
			imp.setStack(stack);
			imp.show();
		} else {
			ImagePlus out = new ImagePlus(imp.getTitle() + "-crop");
			out.setStack(stack);
			out.show();
		}
	}
}
