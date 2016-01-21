package org.doube.util;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.VirtualStack;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;

public class DeleteSliceRange implements PlugIn {
	public void run(final String arg) {
		if (!ImageCheck.checkEnvironment())
			return;
		final ImagePlus imp = IJ.getImage();
		if (null == imp) {
			return;
		}
		final GenericDialog gd = new GenericDialog("Delete slices");
		gd.addMessage("Inclusive range of slices to delete");
		gd.addNumericField("First", imp.getCurrentSlice(), 0);
		gd.addNumericField("Last", imp.getCurrentSlice(), 0);
		gd.showDialog();
		if (gd.wasCanceled()) {
			return;
		}
		final int first = (int) Math.floor(gd.getNextNumber());
		final int last = (int) Math.floor(gd.getNextNumber());

		// check sanity of first and last values
		if (first < 1) {
			IJ.error("First slice cannot be less than 1.");
			return;
		} else if (last > imp.getStackSize()) {
			IJ.error("Last slice cannot be greater than the number of slices.");
			return;
		} else if (first > last) {
			IJ.error("First slice cannot be after last slice");
			return;
		}

		if (imp.getStack().isVirtual()) {
			final VirtualStack stack = (VirtualStack) imp.getStack();
			deleteSliceRange(stack, first, last);
			imp.setStack(null, stack);
		} else {
			final ImageStack stack = imp.getStack();
			deleteSliceRange(stack, first, last);
			imp.setStack(null, stack);
		}
		imp.show();
		UsageReporter.reportEvent(this).send();
	}

	/**
	 * Delete a range of slices from a stack
	 *
	 * @param stack
	 * @param first
	 *            the first slice to remove
	 * @param last
	 *            the last slice to remove
	 */
	public void deleteSliceRange(final ImageStack stack, final int first, final int last) {
		for (int s = first; s <= last; s++) {
			stack.deleteSlice(first);
		}
		return;
	}

	/**
	 * Delete a range of slices from a virtual stack
	 *
	 * @param stack
	 * @param first
	 *            the first slice to remove
	 * @param last
	 *            the last slice to remove
	 */
	public void deleteSliceRange(final VirtualStack stack, final int first, final int last) {
		for (int s = first; s <= last; s++) {
			stack.deleteSlice(first);
		}
		return;
	}
}
