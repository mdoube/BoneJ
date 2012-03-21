package org.doube.util;

import ij.IJ;
import ij.ImagePlus;
import ij.plugin.PlugIn;

public class VoxelDepthChecker implements PlugIn {

	public void run(String arg) {
		ImagePlus imp = IJ.getImage();
		if (null == imp)
			return;
			
		ImageCheck ic = new ImageCheck();
		ic.dicomVoxelDepth(imp);
		UsageReporter.reportEvent(this).send();
		return;
	}

}
