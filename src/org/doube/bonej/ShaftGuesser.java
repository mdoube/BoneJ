package org.doube.bonej;

import org.doube.util.ImageCheck;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;

/**
 * <p>
 * Make an informed guess at which slices of an image stack contain the
 * proximal and distal "ends" of a femoral shaft.
 * </p>
 * 
 * @author Nick Powell
 *
 */

public class ShaftGuesser implements PlugIn {
	
	/** List of slice numbers */
	private double[] slices;
	
	private int al, startSlice, ss, iss;
	
	public void run(String arg) {
		
		/* ImageJ version checking */
		if (!ImageCheck.checkEnvironment())
			return;
		ImagePlus imp = IJ.getImage();
		if (null == imp) {
			IJ.noImage();
			return;
		}
		
		this.al = imp.getStackSize() + 1;
		this.startSlice = 1;
		this.iss = imp.getImageStackSize();
		
		this.slices = new double[this.al];
		for (int s = this.startSlice; s <= this.iss; s++) {
			slices[s] = (double) s;
		}
		
		
		GenericDialog gd = new GenericDialog("Shaft Guesser Options");
		gd.addCheckbox("Proximal end of bone has lower slice number", true);
		gd.addMessage("Choose how to estimate the shaft parameters");
		gd.addNumericField("Smooth over # slices (+/-): ", Math.round(this.al / 50), 0);
		gd.addNumericField("Calculate gradient over # slices: ", Math.round(this.al / 50), 0);
		
		gd.addCheckbox("Graph output", true);
	}

}
