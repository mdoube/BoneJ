package org.doube.util;

import ij.IJ;
import ij.ImagePlus;
import ij.plugin.PlugIn;

public class VoxelDepthChecker implements PlugIn {

	@Override
	public void run(final String arg) {
		final ImagePlus imp = IJ.getImage();
		if (null == imp)
			return;

		final ImageCheck ic = new ImageCheck();
		ImageCheck.dicomVoxelDepth(imp);
		UsageReporter.reportEvent(this).send();
		return;
	}

}
