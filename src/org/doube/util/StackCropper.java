package org.doube.util;

import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
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

		ImageStack stack = RoiMan.cropStack(roiMan, imp.getImageStack(), false,
				0, 0);
		ImagePlus out = new ImagePlus(imp.getTitle()+"-crop");
		out.setStack(stack);
		out.show();
	}
}
