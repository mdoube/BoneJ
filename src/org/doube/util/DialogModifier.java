package org.doube.util;

import ij.IJ;
import ij.gui.GenericDialog;

import java.awt.Checkbox;
import java.awt.Choice;
import java.awt.Component;
import java.awt.Label;
import java.awt.Panel;
import java.awt.TextField;
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
	 * Go through all values in a GenericDialog's Components and call the
	 * appropriate get method. Recursively enter Panel Components. Will throw an
	 * ArrayIndexOutOfBounds exception if gd.getNext... is called elsewhere in
	 * dialogItemChanged().
	 * 
	 * @param gd
	 * @param comps
	 */
	public static void registerMacroValues(GenericDialog gd, Component[] comps) {
		try {
			for (Component c : comps) {
				if (c instanceof Checkbox)
					gd.getNextBoolean();
				else if (c instanceof Choice)
					gd.getNextChoice();
				else if (c instanceof TextField) {
					String text = ((TextField) c).getText();
					try {
						Double.parseDouble(text);
						gd.getNextNumber();
					} catch (NumberFormatException e) {
						gd.getNextString();
					}
				} else if (c instanceof Panel)
					registerMacroValues(gd, ((Panel) c).getComponents());
				else
					continue;
			}
		} catch (Exception e) {
			IJ.log("This plugin causes an exception\n" + e.toString());
		}
		return;
	}

	/**
	 * Check all the numeric text fields in a dialog and return false if any of
	 * them cannot be parsed into a number. Accepts any decimal number,
	 * "Infinity" and "NaN". Rejects strings of 0 length or that contain any
	 * non-decimal characters.
	 * 
	 * 
	 * @param textFields
	 *            e.g. result of GenericDialog.getNumericFields();
	 * @return true if all numeric text fields contain a valid
	 *         number
	 */
	public static boolean allNumbersValid(Vector<?> textFields) {
		for (Object text : textFields) {
			String string = ((TextField) text).getText();
			if (string.length() == 0)
				return false;
			try {
				Double.parseDouble(string);
			} catch (NumberFormatException e) {
				return false;
			}
		}
		return true;
	}
}
