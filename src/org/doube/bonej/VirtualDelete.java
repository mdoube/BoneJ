package org.doube.bonej;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.VirtualStack;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;

public class VirtualDelete implements PlugIn {
	public void run(String arg) {

		ImagePlus imp = IJ.getImage();
		if (null == imp || !imp.getStack().isVirtual()) {
			IJ.error("Plugin requires a virtual stack");
			return;
		}
		VirtualStack stack = (VirtualStack) imp.getStack();
		GenericDialog gd = new GenericDialog("Delete slices");
		gd.addMessage("Inclusive range of slices to delete");
		gd.addNumericField("First", imp.getCurrentSlice(), 0);
		gd.addNumericField("Last", imp.getCurrentSlice(), 0);
		gd.showDialog();
		if (gd.wasCanceled()) {
			return;
		}
		int first = (int)Math.floor(gd.getNextNumber());
		int last = (int)Math.floor(gd.getNextNumber());
		
		deleteVirtualSlices(stack, first, last);
		imp.setStack(null, stack);
		imp.show();
	}

	/**
	 * Delete a range of slices from a virtual stack
	 * @param first
	 * @param last
	 */
	public void deleteVirtualSlices(VirtualStack stack, int first, int last) {
		for (int s = first; s <= last; s++){
			stack.deleteSlice(first);
		}
	}
}
