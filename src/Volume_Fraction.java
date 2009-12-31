import java.awt.Rectangle;

import org.doube.bonej.ResultInserter;
import org.doube.util.ImageCheck;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.measure.Calibration;
import ij.process.ImageProcessor;
import ij.plugin.PlugIn;
import ij.gui.*;

public class Volume_Fraction implements PlugIn {

	public void run(String arg) {
		if (!ImageCheck.checkIJVersion())
			return;
		final ImagePlus imp = IJ.getImage();
		if (null == imp) {
			IJ.noImage();
			return;
		}
		ImageProcessor ip = imp.getProcessor();
		final ImageStack stack = imp.getImageStack();
		final ImageProcessor mask = ip.getMask();
		final boolean hasMask = (mask != null);
		Rectangle r = ip.getRoi();
		final int rLeft = r.x;
		final int rTop = r.y;
		final int rRight = rLeft + r.width;
		final int rBottom = rTop + r.height;
		final int nSlices = imp.getStackSize();

		long volTotal = 0;
		long volBone = 0;
		final double[] thresholds = setThreshold(imp);
		final double minT = thresholds[0];
		final double maxT = thresholds[1];
		for (int s = 1; s <= nSlices; s++) {
			ImageProcessor ipSlice = stack.getProcessor(s);
			for (int v = rTop; v < rBottom; v++) {
				final int vrTop = v - rTop;
				for (int u = rLeft; u < rRight; u++) {
					if (!hasMask || mask.get(u - rLeft, vrTop) > 0) {
						volTotal++;
						final double pixel = ipSlice.get(u, v);
						if (pixel >= minT && pixel <= maxT) {
							volBone++;
						}
					}
				}
			}
		}
		double p = (double) volBone / (double) volTotal;
		Calibration cal = imp.getCalibration();
		double voxelVol = cal.pixelWidth * cal.pixelHeight * cal.pixelDepth;
		ResultInserter ri = ResultInserter.getInstance();
		ri.setResultInRow(imp, "BV (" + cal.getUnits() + "^3)", volBone
				* voxelVol);
		ri.setResultInRow(imp, "TV (" + cal.getUnits() + "^3)", volTotal
				* voxelVol);
		ri.setResultInRow(imp, "BV/TV", p);
		ri.updateTable();
		return;
	}

	private double[] setThreshold(ImagePlus imp) {
		double[] thresholds = new double[2];
		ImageCheck ic = new ImageCheck();
		if (ic.isBinary(imp)) {
			thresholds[0] = 128;
			thresholds[1] = 255;
		} else {
			IJ.run("Threshold...");
			new WaitForUserDialog("Set the threshold, then click OK.").show();
			thresholds[0] = imp.getProcessor().getMinThreshold();
			thresholds[1] = imp.getProcessor().getMaxThreshold();
		}
		return thresholds;
	}
}
