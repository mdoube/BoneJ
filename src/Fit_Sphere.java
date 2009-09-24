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
		if (IJ.versionLessThan("1.42h")) {
			IJ.error("Your version of ImageJ is too old/n"
					+ "Please update (Help->Update ImageJ...)");
			return DONE;
		}
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
		return DOES_ALL + STACK_REQUIRED;
	}

	public void run(ImageProcessor ip) {
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
		if (doCopy)
			copySphere(imp, ip, padding, cropFactor, sphereDim);
		if (doInnerCube)
			copyInnerCube(imp, ip, cropFactor, sphereDim);
		if (doOuterCube)
			copyOuterCube(imp, ip, cropFactor, sphereDim);

		String units = imp.getCalibration().getUnits();
		ResultInserter ri = new ResultInserter();
		ri.setResultInRow(imp, "X centroid (" + units + ")", sphereDim[0]);
		ri.setResultInRow(imp, "Y centroid (" + units + ")", sphereDim[1]);
		ri.setResultInRow(imp, "Z centroid (" + units + ")", sphereDim[2]);
		ri.setResultInRow(imp, "Radius (" + units + ")", sphereDim[3]);
		ri.updateTable();
	}

	public boolean showDialog() {
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

	public static double[] getVoxDim(ImagePlus imp) {
		Calibration cal = imp.getCalibration();
		double[] voxDim = new double[3];
		voxDim[0] = cal.pixelWidth;
		voxDim[1] = cal.pixelHeight;
		voxDim[2] = cal.pixelDepth;
		return voxDim;
	}

	// TODO make this go faster by getting slice pixels and iterating through
	// it's array rather than using setSlice and getPixel
	public void copySphere(ImagePlus imp, ImageProcessor ip, int padding,
			double cropFactor, double[] sphereDim) {
		double[] voxDim = getVoxDim(imp);
		int startX = (int) Math
				.round((sphereDim[0] - sphereDim[3] * cropFactor) / voxDim[0])
				- padding;
		int startY = (int) Math
				.round((sphereDim[1] - sphereDim[3] * cropFactor) / voxDim[1])
				- padding;
		int startZ = (int) Math
				.round((sphereDim[2] - sphereDim[3] * cropFactor) / voxDim[2])
				- padding;
		int roiWidth = (int) Math.round(2 * sphereDim[3] * cropFactor
				/ voxDim[0])
				+ 2 * padding;
		int roiHeight = (int) Math.round(2 * sphereDim[3] * cropFactor
				/ voxDim[1])
				+ 2 * padding;
		int roiDepth = (int) Math.round(2 * sphereDim[3] * cropFactor
				/ voxDim[2])
				+ 2 * padding;
		ImageStack sourceStack = imp.getImageStack();
		ImageStack targetStack = new ImageStack(roiWidth, roiHeight);
		for (int z = startZ; z <= startZ + roiDepth; z++) {
			IJ.showProgress(z - startZ, roiDepth);
			IJ.showStatus("Copying sphere to new stack");
			short[] targetSlice = new short[roiWidth * roiHeight];
			int nRows = 0;
			for (int y = startY; y < startY + roiHeight; y++) {
				int index = nRows * roiWidth;
				int nCols = 0;
				ip = sourceStack.getProcessor(z);
				for (int x = startX; x < startX + roiWidth; x++) {
					double distance = Math.sqrt((x * voxDim[0] - sphereDim[0])
							* (x * voxDim[0] - sphereDim[0])
							+ (y * voxDim[1] - sphereDim[1])
							* (y * voxDim[1] - sphereDim[1])
							+ (z * voxDim[2] - sphereDim[2])
							* (z * voxDim[2] - sphereDim[2]));
					if (distance < sphereDim[3] * cropFactor) {
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
		double h = sphereDim[3] * cropFactor / Math.sqrt(3);
		int startX = (int) Math.round((sphereDim[0] - h) / voxDim[0]);
		int startY = (int) Math.round((sphereDim[1] - h) / voxDim[1]);
		int startZ = (int) Math.round((sphereDim[2] - h) / voxDim[2]);
		int roiWidth = (int) Math.round(2 * h / voxDim[0]);
		int roiHeight = (int) Math.round(2 * h / voxDim[1]);
		int roiDepth = (int) Math.round(2 * h / voxDim[2]);
		ImageStack sourceStack = imp.getStack();
		ImageStack targetStack = new ImageStack(roiWidth, roiHeight);
		for (int z = startZ; z <= startZ + roiDepth; z++) {
			IJ.showProgress(z - startZ, roiDepth);
			IJ.showStatus("Copying largest enclosed cube");
			short[] targetSlice = new short[roiWidth * roiHeight];
			int nRows = 0;
			for (int y = startY; y < startY + roiHeight; y++) {
				int index = nRows * roiWidth;
				int nCols = 0;
				ip = sourceStack.getProcessor(z);
				for (int x = startX; x < startX + roiWidth; x++) {
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
		double h = sphereDim[3] * cropFactor;
		int startX = (int) Math.round((sphereDim[0] - h) / voxDim[0]);
		int startY = (int) Math.round((sphereDim[1] - h) / voxDim[1]);
		int startZ = (int) Math.round((sphereDim[2] - h) / voxDim[2]);
		int roiWidth = (int) Math.round(2 * h / voxDim[0]);
		int roiHeight = (int) Math.round(2 * h / voxDim[1]);
		int roiDepth = (int) Math.round(2 * h / voxDim[2]);
		ImageStack sourceStack = imp.getImageStack();
		ImageStack targetStack = new ImageStack(roiWidth, roiHeight);
		for (int z = startZ; z <= startZ + roiDepth; z++) {
			IJ.showProgress(z - startZ, roiDepth);
			IJ.showStatus("Copying smallest enclosing cube");
			short[] targetSlice = new short[roiWidth * roiHeight];
			int nRows = 0;
			for (int y = startY; y < startY + roiHeight; y++) {
				int index = nRows * roiWidth;
				int nCols = 0;
				ip = sourceStack.getProcessor(z);
				for (int x = startX; x < startX + roiWidth; x++) {
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