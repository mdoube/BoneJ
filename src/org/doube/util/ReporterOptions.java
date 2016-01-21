package org.doube.util;

import java.util.Random;

import ij.Prefs;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;

public class ReporterOptions implements PlugIn {

	public static final String OPTOUTSET = "bonej.report.option.set";
	/** Set to false if reporting is not allowed */
	public static final String OPTOUTKEY = "bonej.allow.reporter";
	public static final String COOKIE = "bonej.report.cookie";
	public static final String COOKIE2 = "bonej.report.cookie2";
	public static final String FIRSTTIMEKEY = "bonej.report.firstvisit";
	public static final String SESSIONKEY = "bonej.report.bonejsession";
	public static final String IJSESSIONKEY = "bonej.report.ijsession";

	public void run(final String arg) {

		final GenericDialog dialog = new GenericDialog("BoneJ");
		dialog.addMessage("Allow usage data collection?");
		dialog.addMessage("BoneJ would like to collect data on \n"
				+ "which plugins are being used, to direct development\n" + "and promote BoneJ to funders.");
		dialog.addMessage(
				"If you agree to participate please hit OK\n" + "otherwise, cancel. For more information click Help.");
		dialog.addHelp("http://bonej.org/stats");
		dialog.showDialog();
		if (dialog.wasCanceled()) {
			Prefs.set(OPTOUTKEY, false);
			Prefs.set(ReporterOptions.COOKIE, "");
			Prefs.set(ReporterOptions.COOKIE2, "");
			Prefs.set(ReporterOptions.FIRSTTIMEKEY, "");
			Prefs.set(ReporterOptions.SESSIONKEY, "");
			Prefs.set(ReporterOptions.IJSESSIONKEY, "");
		} else {
			Prefs.set(OPTOUTKEY, true);
			Prefs.set(ReporterOptions.COOKIE, new Random().nextInt(Integer.MAX_VALUE));
			Prefs.set(ReporterOptions.COOKIE2, new Random().nextInt(Integer.MAX_VALUE));
			final long time = System.currentTimeMillis() / 1000;
			Prefs.set(ReporterOptions.FIRSTTIMEKEY, Long.toString(time));
			Prefs.set(SESSIONKEY, 1);
		}

		Prefs.set(OPTOUTSET, true);
		Prefs.savePreferences();
		UsageReporter.reportEvent(this).send();
		return;
	}
}
