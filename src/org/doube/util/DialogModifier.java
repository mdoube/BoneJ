package org.doube.util;

import ij.IJ;
import ij.gui.GenericDialog;

import java.awt.Label;
import java.awt.Panel;
import java.util.Vector;

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

	/**
	 * Count the number of each input type in a GenericDialog and call
	 * the appropriate get method the correct number of times, registering
	 * the value of the field in the macro recorder.
	 * 
	 * @param gd
	 */
	public static void registerMacroValues(GenericDialog gd) {
		// count the different input types
		Vector<?> numbers = gd.getNumericFields();
		if (numbers != null) {
			IJ.log("Number of numeric fields = " + numbers.size());
			int nNumeric = numbers.size();
			for (int i = 0; i < nNumeric; i++) {
				gd.getNumericFields();
			}
		}
		Vector<?> checkboxes = gd.getCheckboxes();
		if (checkboxes != null) {
			IJ.log("Number of checkboxes = " + checkboxes.size());
			int nCheckboxes = checkboxes.size();
			for (int i = 0; i < nCheckboxes; i++) {
				gd.getNextBoolean();
			}
		}
		Vector<?> choices = gd.getChoices();
		if (choices != null) {
			IJ.log("Number of choices = " + choices.size());
			int nChoices = choices.size();
			for (int i = 0; i < nChoices; i++) {
				gd.getNextChoice();
			}
		}
		Vector<?> strings = gd.getStringFields();
		if (strings != null) {
			IJ.log("Number of strings = " + strings.size());
			int nStrings = strings.size();
			for (int i = 0; i < nStrings; i++) {
				gd.getNextString();
			}
		}
		return;
	}
}
