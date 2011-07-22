package org.doube.util;

import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.measure.Calibration;
import ij.plugin.PlugIn;
import ij.util.DicomTools;

public class DensityCalibrator implements PlugIn {

	public void run(String arg) {
		ImagePlus imp = WindowManager.getCurrentImage();
		if (arg.equals("scanco"))
			try {
				scanco(imp);
			} catch (Exception e) {
				IJ.error(e.getMessage());
				return;
			}
	}

	private void scanco(ImagePlus imp) throws IllegalArgumentException {
		String manufacturer = DicomTools.getTag(imp, "0008,0070");
		IJ.log("manufacturer = "+manufacturer);
		if (!manufacturer.contains("SCANCO")) {
			throw new IllegalArgumentException(
					"File not a SCANCO Medical DICOM");
		}
		double slope = Double.parseDouble(DicomTools.getTag(imp, "0029,1004"));
		double intercept = Double.parseDouble(DicomTools.getTag(imp,
				"0029,1005"));
		double scaling = Double
				.parseDouble(DicomTools.getTag(imp, "0029,1000"));
		double c = intercept - 32768 * slope / scaling;
		double m = slope / scaling;
		double[] coef = {c, m};
		Calibration cal = imp.getCalibration();
		cal.setFunction(Calibration.STRAIGHT_LINE, coef, "mg HA/ccm");
		imp.setCalibration(cal);
		imp.updateAndDraw();
	}

}
