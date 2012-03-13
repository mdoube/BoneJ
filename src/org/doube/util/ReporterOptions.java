package org.doube.util;

import ij.Prefs;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;

public class ReporterOptions implements PlugIn {
	
	public static final String OPTOUTSET = "bonej.report.option.set";
	public static final String OPTOUTKEY = "bonej.allow.reporter";
	
	public void run(String arg) {

		GenericDialog dialog = new GenericDialog("BoneJ");
		dialog.addMessage("Allow usage data collection?");
		dialog.addMessage("BoneJ would like to collect data on \n"
				+ "which plugins are being used, to direct development\n"
				+ "and promote BoneJ to funders.");
		dialog.addMessage("If you agree to participate please hit OK\n"
				+ "otherwise, cancel. For more information click Help.");
		dialog.addHelp("http://bonej.org/stats");
		dialog.showDialog();
		if (dialog.wasCanceled())
			Prefs.set(OPTOUTKEY, false);
		else
			Prefs.set(OPTOUTKEY, true);

		Prefs.set(OPTOUTSET , true);
		Prefs.savePreferences();
		return;
	}

}
