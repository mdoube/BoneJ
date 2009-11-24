/**
 * Fit_Sphere plugin for ImageJ
 * Copyright 2008 2009 Michael Doube
 * 
 *This program is free software: you can redistribute it and/or modify
 *it under the terms of the GNU General Public License as published by
 *the Free Software Foundation, either version 3 of the License, or
 *(at your option) any later version.
 *
 *This program is distributed in the hope that it will be useful,
 *but WITHOUT ANY WARRANTY; without even the implied warranty of
 *MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *GNU General Public License for more details.
 *
 *You should have received a copy of the GNU General Public License
 *along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ImageProcessor;
import ij.plugin.filter.PlugInFilter;
import ij.gui.*;
import ij.plugin.frame.*;
import ij.measure.Calibration;

import org.doube.bonej.ImageCheck;
import org.doube.bonej.ResultInserter;
import org.doube.bonej.FitSphere;

/**
 *<p>
 * Takes point selections from ROI manager and returns the centroid and radius
 * of a best fit sphere Ported from Angelo Tardugno's C++
 * </p>
 * 
 * 
 *@author Michael Doube and Angelo Tardugno
 *@version 0.1
 */
public class Fit_Sphere implements PlugInFilter {
	ImagePlus imp;
	ImageProcessor ip;
	RoiManager roiMan = RoiManager.getInstance();
	protected ImageStack sourceStack;
	public boolean doCopy, doInnerCube, doOuterCube;
	public int padding;
	public double cropFactor;

	public int setup(String arg, ImagePlus imp) {
		this.imp = imp;
		if (imp == null || imp.getNSlices() < 2) {
			IJ.showMessage("A stack must be open");
			return DONE;
		}
		if (roiMan == null && imp != null) {
			IJ.run("ROI Manager...");
			IJ.error("Please populate ROI Manager with point ROIs");
			return DONE;
		}
		sourceStack = imp.getStack();
		return DOES_ALL + STACK_REQUIRED + NO_CHANGES;
	}

	public void run(ImageProcessor ip) {
		if (!ImageCheck.checkIJVersion())
			return;
		double[] voxDim = getVoxDim(imp);
		if (voxDim[2] > voxDim[0] * 20) {
			if (!IJ
					.showMessageWithCancel(
							"Voxel depth problem",
							"The voxel depth (slice thickness)"
									+ "is unusually large.\nClick OK if you know it is correct.")) {
				imp.unlock();
				IJ.run("Properties...");
			}
		}
		if (!imp.lock())
			imp.lock(); // if we have unlocked the image to reset properties,
		// relock it.
		if (!showDialog()) {
			return;
		}
		FitSphere fs = new FitSphere();
		double[][] points = fs.getRoiManPoints(imp, roiMan);
		double[] sphereDim = fs.fitSphere(points);

		String units = imp.getCalibration().getUnits();
		ResultInserter ri = ResultInserter.getInstance();
		ri.setResultInRow(imp, "X centroid (" + units + ")", sphereDim[0]);
		ri.setResultInRow(imp, "Y centroid (" + units + ")", sphereDim[1]);
		ri.setResultInRow(imp, "Z centroid (" + units + ")", sphereDim[2]);
		ri.setResultInRow(imp, "Radius (" + units + ")", sphereDim[3]);
		ri.updateTable();

		if (doCopy)
			copySphere(imp, ip, padding, cropFactor, sphereDim);
		if (doInnerCube)
			copyInnerCube(imp, ip, cropFactor, sphereDim);
		if (doOuterCube)
			copyOuterCube(imp, ip, cropFactor, sphereDim);
	}

	private boolean showDialog() {
		GenericDialog gd = new GenericDialog("Setup");
		gd.addMessage("");
		gd.addCheckbox("Copy Sphere", true);
		gd.addCheckbox("Inner Cube", true);
		gd.addCheckbox("Outer Cube", true);
		gd.addNumericField("Padding", 2, 0, 2, "voxels");
		gd.addNumericField("Crop Factor", 1.0, 2, 4, "");
		gd.showDialog();
		if (gd.wasCanceled()) {
			return false;
		}
		doCopy = gd.getNextBoolean();
		doInnerCube = gd.getNextBoolean();
		doOuterCube = gd.getNextBoolean();
		padding = (int) gd.getNextNumber();
		cropFactor = gd.getNextNumber();
		return true;
	}

	private static double[] getVoxDim(ImagePlus imp) {
		Calibration cal = imp.getCalibration();
		double[] voxDim = { cal.pixelWidth, cal.pixelHeight, cal.pixelDepth };
		return voxDim;
	}

	public void copySphere(ImagePlus imp, ImageProcessor ip, int padding,
			double cropFactor, double[] sphereDim) {
		double[] voxDim = getVoxDim(imp);
		final double vW = voxDim[0];
		final double vH = voxDim[1];
		final double vD = voxDim[2];
		final double xC = sphereDim[0];
		final double yC = sphereDim[2];
		final double zC = sphereDim[2];
		final double r = sphereDim[3];
		int startX = (int) Math.round((xC - r * cropFactor) / vW) - padding;
		int startY = (int) Math.round((yC - r * cropFactor) / vH) - padding;
		int startZ = (int) Math.round((zC - r * cropFactor) / vD) - padding;
		int roiWidth = (int) Math.round(2 * r * cropFactor / vW) + 2 * padding;
		int roiHeight = (int) Math.round(2 * r * cropFactor / vH) + 2 * padding;
		int roiDepth = (int) Math.round(2 * r * cropFactor / vD) + 2 * padding;
		// bounds checking
		if (startX < 0)
			startX = 0;
		if (startY < 0)
			startY = 0;
		if (startZ < 1)
			startZ = 1;
		if (startX + roiWidth > imp.getWidth())
			roiWidth = imp.getWidth() - startX;
		if (startY + roiHeight > imp.getHeight())
			roiHeight = imp.getHeight() - startY;
		if (startZ + roiDepth > imp.getStackSize())
			roiDepth = imp.getStackSize() - startZ;
		ImageStack sourceStack = imp.getImageStack();
		ImageStack targetStack = new ImageStack(roiWidth, roiHeight);
		final int endZ = startZ + roiDepth;
		final int roiArea = roiWidth * roiHeight;
		final int endY = startY + roiHeight;
		final int endX = startX + roiWidth;
		final double maxRad = r * cropFactor;
		for (int z = startZ; z <= endZ; z++) {
			IJ.showProgress(z - startZ, roiDepth);
			IJ.showStatus("Copying sphere to new stack");
			short[] targetSlice = new short[roiArea];
			int nRows = 0;
			ip = sourceStack.getProcessor(z);
			for (int y = startY; y < endY; y++) {
				int index = nRows * roiWidth;
				int nCols = 0;
				for (int x = startX; x < endX; x++) {
					double distance = Math.sqrt((x * vW - xC) * (x * vW - xC)
							+ (y * vH - yC) * (y * vH - yC) + (z * vD - zC)
							* (z * vD - zC));
					if (distance < maxRad) {
						targetSlice[index + nCols] = (short) ip.get(x, y);
					} else {
						targetSlice[index + nCols] = 0;
					}
					nCols++;
				}
				nRows++;
			}
			targetStack.addSlice(sourceStack.getSliceLabel(z), targetSlice);
		}
		ImagePlus target = new ImagePlus("Sphere", targetStack);
		target.setCalibration(imp.getCalibration());
		target.setDisplayRange(imp.getDisplayRangeMin(), imp
				.getDisplayRangeMax());
		target.show();
		return;
	}

	// TODO make this go faster by getting slice pixels and iterating through
	// its array rather than using setSlice and getPixel
	public void copyInnerCube(ImagePlus imp, ImageProcessor ip,
			double cropFactor, double[] sphereDim) {
		Calibration cal = imp.getCalibration();
		double[] voxDim = getVoxDim(imp);
		final double vW = voxDim[0];
		final double vH = voxDim[1];
		final double vD = voxDim[2];
		final double xC = sphereDim[0];
		final double yC = sphereDim[2];
		final double zC = sphereDim[2];
		final double r = sphereDim[3];
		final double h = r * cropFactor / Math.sqrt(3);
		int startX = (int) Math.round((xC - h) / vW);
		int startY = (int) Math.round((yC - h) / vH);
		int startZ = (int) Math.round((zC - h) / vD);
		int roiWidth = (int) Math.round(2 * h / vW);
		int roiHeight = (int) Math.round(2 * h / vH);
		int roiDepth = (int) Math.round(2 * h / vD);
		// bounds checking
		if (startX < 0)
			startX = 0;
		if (startY < 0)
			startY = 0;
		if (startZ < 1)
			startZ = 1;
		if (startX + roiWidth > imp.getWidth())
			roiWidth = imp.getWidth() - startX;
		if (startY + roiHeight > imp.getHeight())
			roiHeight = imp.getHeight() - startY;
		if (startZ + roiDepth > imp.getStackSize())
			roiDepth = imp.getStackSize() - startZ;
		ImageStack sourceStack = imp.getStack();
		ImageStack targetStack = new ImageStack(roiWidth, roiHeight);
		final int endZ = startZ + roiDepth;
		final int roiArea = roiWidth * roiHeight;
		final int endY = startY + roiHeight;
		final int endX = startX + roiWidth;
		for (int z = startZ; z <= endZ; z++) {
			IJ.showProgress(z - startZ, roiDepth);
			IJ.showStatus("Copying largest enclosed cube");
			short[] targetSlice = new short[roiArea];
			int nRows = 0;
			ip = sourceStack.getProcessor(z);
			for (int y = startY; y < endY; y++) {
				int index = nRows * roiWidth;
				int nCols = 0;
				for (int x = startX; x < endX; x++) {
					targetSlice[index + nCols] = (short) ip.get(x, y);
					nCols++;
				}
				nRows++;
			}
			targetStack.addSlice(sourceStack.getSliceLabel(z), targetSlice);
		}
		ImagePlus target = new ImagePlus("Inner Cube", targetStack);
		target.setCalibration(cal);
		target.setDisplayRange(imp.getDisplayRangeMin(), imp
				.getDisplayRangeMax());
		target.show();
		return;
	}

	public void copyOuterCube(ImagePlus imp, ImageProcessor ip,
			double cropFactor, double[] sphereDim) {
		Calibration cal = imp.getCalibration();
		double[] voxDim = getVoxDim(imp);
		final double vW = voxDim[0];
		final double vH = voxDim[1];
		final double vD = voxDim[2];
		final double xC = sphereDim[0];
		final double yC = sphereDim[2];
		final double zC = sphereDim[2];
		final double r = sphereDim[3];
		final double h = r * cropFactor;
		int startX = (int) Math.round((xC - h) / vW);
		int startY = (int) Math.round((yC - h) / vH);
		int startZ = (int) Math.round((zC - h) / vD);
		int roiWidth = (int) Math.round(2 * h / vW);
		int roiHeight = (int) Math.round(2 * h / vH);
		int roiDepth = (int) Math.round(2 * h / vD);
		// bounds checking
		if (startX < 0)
			startX = 0;
		if (startY < 0)
			startY = 0;
		if (startZ < 1)
			startZ = 1;
		if (startX + roiWidth > imp.getWidth())
			roiWidth = imp.getWidth() - startX;
		if (startY + roiHeight > imp.getHeight())
			roiHeight = imp.getHeight() - startY;
		if (startZ + roiDepth > imp.getStackSize())
			roiDepth = imp.getStackSize() - startZ;
		ImageStack sourceStack = imp.getImageStack();
		ImageStack targetStack = new ImageStack(roiWidth, roiHeight);
		final int endZ = startZ + roiDepth;
		final int roiArea = roiWidth * roiHeight;
		final int endY = startY + roiHeight;
		final int endX = startX + roiWidth;
		for (int z = startZ; z <= endZ; z++) {
			IJ.showProgress(z - startZ, roiDepth);
			IJ.showStatus("Copying smallest enclosing cube");
			short[] targetSlice = new short[roiArea];
			int nRows = 0;
			ip = sourceStack.getProcessor(z);
			for (int y = startY; y < endY; y++) {
				int index = nRows * roiWidth;
				int nCols = 0;
				for (int x = startX; x < endX; x++) {
					targetSlice[index + nCols] = (short) ip.get(x, y);
					nCols++;
				}
				nRows++;
			}
			targetStack.addSlice(sourceStack.getSliceLabel(z), targetSlice);
		}
		ImagePlus target = new ImagePlus("Outer Cube", targetStack);
		target.setCalibration(cal);
		target.setDisplayRange(imp.getDisplayRangeMin(), imp
				.getDisplayRangeMax());
		target.show();
		return;
	}
}