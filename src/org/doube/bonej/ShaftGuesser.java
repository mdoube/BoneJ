package org.doube.bonej;

import org.doube.util.DeleteSliceRange;
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
	
	private boolean proximalLow, doGraph;
	
	private double smoothOver;
	private double gradientOver;
	
	/** List of slice numbers */
	private double[] slices;
	
	private int al, startSlice, endSlice, ss, iss;
	
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
		this.iss = imp.getImageStackSize();
		
		this.slices = new double[this.al];
		for (int s = this.startSlice; s <= this.iss; s++) {
			slices[s] = (double) s;
		}
		
		GenericDialog gd = new GenericDialog("Shaft Guesser Options");
		gd.addCheckbox("Proximal end of bone has lower slice number", true);
		gd.addNumericField("Bone start slice", 1, 2);
		gd.addNumericField("Bone end slice", this.iss, 2);
		gd.addMessage("Choose how to estimate the shaft parameters: ");
		gd.addNumericField("Smooth over # slices (+/-): ", Math.round(this.al / 50), 0);
		gd.addNumericField("Calculate gradient over # slices (+/-): ", Math.round(this.al / 50), 0);
		gd.addCheckbox("Graph output", true);
		gd.showDialog();
		
		this.proximalLow = gd.getNextBoolean();
		this.startSlice = (int) gd.getNextNumber();
		this.endSlice = (int) gd.getNextNumber();
		this.smoothOver = gd.getNextNumber();
		this.gradientOver = gd.getNextNumber();
		this.doGraph = gd.getNextBoolean();
		
		/* Initial setup */
		if(!proximalLow) {
			IJ.run("Flip Z");
		}
		
		if(this.startSlice > 1) {
			DeleteSliceRange dsr1 = new DeleteSliceRange();
			dsr1.deleteSliceRange(imp.getImageStack(), 1, this.startSlice - 1);
			this.startSlice = 1;
		}
		if(this.endSlice < this.iss) {
			DeleteSliceRange dsr2 = new DeleteSliceRange();
			dsr2.deleteSliceRange(imp.getImageStack(), this.endSlice + 1, this.iss);
			this.endSlice = this.iss;
		}
		
	}

}
