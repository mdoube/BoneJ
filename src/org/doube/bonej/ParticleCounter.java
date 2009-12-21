package org.doube.bonej;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
import org.doube.bonej.Purify;

public class ParticleCounter implements PlugIn{

	
	public void run(String arg) {
		ImagePlus imp = IJ.getImage();
		if (null == imp){
			IJ.noImage();
			return;
		}
		ImageCheck ic = new ImageCheck();
		if (!ic.isBinary(imp)){
			IJ.error("Binary image required");
			return;	
		}
		GenericDialog gd = new GenericDialog("Setup");
		gd.addNumericField("Slices per chunk", 2, 0);
		if (gd.wasCanceled()){
			return;
		}
		final int slicesPerChunk = (int) Math.floor(gd.getNextNumber());
		
		Purify purify = new Purify();
		Object[] result = purify.getParticles(imp, slicesPerChunk, Purify.FORE);
		byte[][] workArray = (byte[][]) result[0];
		int[][] particleLabels = (int[][]) result[1];
		long[] particleSizes = (long[]) result[2];
		final int w = imp.getWidth();
		final int h = imp.getHeight();
	}
}
