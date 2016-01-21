package org.doube.util;

import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.measure.Calibration;
import ij.plugin.PlugIn;
import ij.util.DicomTools;

public class DensityCalibrator implements PlugIn {

	public void run(final String arg) {
		final ImagePlus imp = WindowManager.getCurrentImage();
		if (imp == null) {
			IJ.noImage();
			return;
		}
		if (arg.equals("scanco"))
			try {
				scanco(imp);
			} catch (final NumberFormatException e) {
				IJ.error("Calibration data missing from DICOM header");
				return;
			} catch (final NullPointerException e) {
				IJ.error("Calibration data missing from DICOM header");
				return;
			} catch (final IllegalArgumentException e) {
				IJ.error(e.getMessage());
				return;
			} catch (final Exception e) {
				IJ.error("Can't calibrate image\n" + e.getMessage());
				return;
			}
		UsageReporter.reportEvent(this).send();
	}

	private void scanco(final ImagePlus imp) throws IllegalArgumentException {
		final String manufacturer = DicomTools.getTag(imp, "0008,0070");
		if (manufacturer == null || !manufacturer.contains("SCANCO")) {
			throw new IllegalArgumentException("File is not a SCANCO Medical DICOM");
		}
		final double slope = Double.parseDouble(DicomTools.getTag(imp, "0029,1004"));
		final double intercept = Double.parseDouble(DicomTools.getTag(imp, "0029,1005"));
		final double scaling = Double.parseDouble(DicomTools.getTag(imp, "0029,1000"));
		final double c = intercept - 32768 * slope / scaling;
		final double m = slope / scaling;
		final double[] coef = { c, m };
		final Calibration cal = imp.getCalibration();
		cal.setFunction(Calibration.STRAIGHT_LINE, coef, "mg HA/ccm");
		imp.setCalibration(cal);
		imp.updateAndDraw();
	}

}
