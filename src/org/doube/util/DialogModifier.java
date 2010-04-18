package org.doube.util;

import ij.gui.GenericDialog;

import java.awt.Label;
import java.awt.Panel;

public class DialogModifier {
	/**
	 * Replace the unit string to the right of all numeric textboxes in a
	 * GenericDialog
	 * 
	 * @param gd
	 *            the dialog window
	 * @param oldUnits
	 *            original unit string
	 * @param newUnits
	 *            new unit string
	 */
	public static void replaceUnitString(GenericDialog gd, String oldUnits,
			String newUnits) {
		for (int n = 0; n < gd.getComponentCount(); n++) {
			if (gd.getComponent(n) instanceof Panel) {
				Panel panel = (Panel) gd.getComponent(n);
				if (panel.getComponent(1) instanceof Label) {
					Label u = (Label) panel.getComponent(1);
					String unitString = u.getText();
					u.setText(unitString.replace(oldUnits, newUnits));
				}
			}
		}
	}

}
