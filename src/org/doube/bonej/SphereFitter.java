package org.doube.bonej;

/**
 * SphereFitter plugin for ImageJ
 * Copyright 2008 2009 2010 Michael Doube
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

import java.awt.AWTEvent;
import java.awt.Checkbox;
import java.awt.TextField;
import java.util.Vector;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
//import ij.WindowManager;
import ij.process.ImageProcessor;
import ij.plugin.PlugIn;
import ij.gui.*;
import ij.plugin.frame.*;
import ij.measure.Calibration;

import org.doube.geometry.FitSphere;
import org.doube.util.DialogModifier;
import org.doube.util.ImageCheck;
import org.doube.util.ResultInserter;
import org.doube.util.RoiMan;
import org.doube.util.UsageReporter;

/**
 * <p>
 * Takes point selections from ROI manager and returns the centroid and radius
 * of a best fit sphere Ported from Angelo Tardugno's C++
 * </p>
 * 
 * 
 * @author Michael Doube and Angelo Tardugno
 */
public class SphereFitter implements PlugIn, DialogListener {

	public void run(String arg) {
		if (!ImageCheck.checkEnvironment())
			return;
		final ImagePlus imp = IJ.getImage();
		if (null == imp) {
			IJ.noImage();
			return;
		}
		ImageCheck ic = new ImageCheck();
		if (!ic.isMultiSlice(imp)) {
			IJ.error("Stack required");
			return;
		}
		RoiManager roiMan = RoiManager.getInstance();
		if (roiMan == null && imp != null) {
			IJ.error("Please populate ROI Manager with point ROIs");
			IJ.run("ROI Manager...");
			return;
		}

		GenericDialog gd = new GenericDialog("Setup");
		gd.addCheckbox("Copy Sphere", true);
		gd.addNumericField("Padding", 2, 0, 2, "voxels");
		gd.addCheckbox("Inner Cube", true);
		gd.addCheckbox("Outer Cube", true);
		gd.addNumericField("Crop Factor", 1.0, 2, 4, "");
		gd.addCheckbox("Add_to_ROI_Manager", false);
		gd.addCheckbox("Clear_ROI_Manager", false);
		gd.addHelp("http://bonej.org/sphere");
		gd.addDialogListener(this);
		gd.showDialog();
		if (gd.wasCanceled()) {
			return;
		}
		final boolean doCopy = gd.getNextBoolean();
		final int padding = (int) gd.getNextNumber();
		final boolean doInnerCube = gd.getNextBoolean();
		final boolean doOuterCube = gd.getNextBoolean();
		final double cropFactor = gd.getNextNumber();
		final boolean doRoiMan = gd.getNextBoolean();
		final boolean clearRois = gd.getNextBoolean();

		final double[][] points = RoiMan.getRoiManPoints(imp, roiMan);
		double i = points[0][2];
		double j = i;
		for (double[] p : points)
			if (p[2] != i)
				j = p[2];
		if (j == i){
			IJ.showMessage("Can't fit sphere to points.\n"
					+ "All points are on the same slice.\n" +
					"Check the ROI Manager option More >> Options... >\n" +
					"Associate ROis with stack positions.\n");
			return;
		}
			
		
		double[] sphereDim = new double[4];
		try {
			sphereDim = FitSphere.fitSphere(points);
		} catch (IllegalArgumentException ia) {
			IJ.showMessage(ia.getMessage());
			return;
		} catch (RuntimeException re) {
			IJ.showMessage("Can't fit sphere to points.\n"
					+ "Add more point ROI's to the ROI Manager and try again.\n");
			return;
		}

		String units = imp.getCalibration().getUnits();
		ResultInserter ri = ResultInserter.getInstance();
		ri.setResultInRow(imp, "X centroid (" + units + ")", sphereDim[0]);
		ri.setResultInRow(imp, "Y centroid (" + units + ")", sphereDim[1]);
		ri.setResultInRow(imp, "Z centroid (" + units + ")", sphereDim[2]);
		ri.setResultInRow(imp, "Radius (" + units + ")", sphereDim[3]);
		ri.updateTable();

		if (doCopy)
			copySphere(imp, padding, cropFactor, sphereDim).show();
		if (doInnerCube)
			copyInnerCube(imp, cropFactor, sphereDim).show();
		if (doOuterCube)
			copyOuterCube(imp, cropFactor, sphereDim).show();
		if (doRoiMan) {
			addToRoiManager(imp, roiMan, sphereDim, clearRois);
		}
		UsageReporter.reportEvent(this).send();
		return;
	}

	/**
	 * 
	 * @param imp
	 * @param padding
	 * @param cropFactor
	 * @param sphereDim
	 * @return
	 */
	public ImagePlus copySphere(ImagePlus imp, int padding, double cropFactor,
			double[] sphereDim) {
		Calibration cal = imp.getCalibration();
		final double vW = cal.pixelWidth;
		final double vH = cal.pixelHeight;
		final double vD = cal.pixelDepth;
		final double xC = sphereDim[0];
		final double yC = sphereDim[1];
		final double zC = sphereDim[2];
		final double r = sphereDim[3];
		final double maxRad = r * cropFactor;
		int startX = (int) Math.round((xC - maxRad) / vW) - padding;
		int startY = (int) Math.round((yC - maxRad) / vH) - padding;
		int startZ = (int) Math.round((zC - maxRad) / vD) - padding;
		int roiWidth = (int) Math.round(2 * maxRad / vW) + 2 * padding;
		int roiHeight = (int) Math.round(2 * maxRad / vH) + 2 * padding;
		int roiDepth = (int) Math.round(2 * maxRad / vD) + 2 * padding;
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
		final ImageStack sourceStack = imp.getImageStack();
		ImageStack targetStack = new ImageStack(roiWidth, roiHeight, roiDepth);
		final int endZ = startZ + roiDepth;
		final int endY = startY + roiHeight;
		final int endX = startX + roiWidth;
		for (int z = startZ; z < endZ; z++) {
			IJ.showProgress(z - startZ, roiDepth);
			IJ.showStatus("Copying sphere to new stack");
			final int tZ = z - startZ + 1;
			final double dZ = z * vD - zC;
			final double dZ2 = dZ * dZ;
			targetStack.setPixels(
					Moments.getEmptyPixels(roiWidth, roiHeight,
							imp.getBitDepth()), tZ);
			targetStack.setSliceLabel(sourceStack.getShortSliceLabel(z), tZ);
			final ImageProcessor ip = sourceStack.getProcessor(z);
			ImageProcessor targetIP = targetStack.getProcessor(tZ);
			for (int y = startY; y < endY; y++) {
				final int tY = y - startY;
				final double dY = y * vH - yC;
				final double dY2 = dY * dY;
				for (int x = startX; x < endX; x++) {
					final double dX = x * vW - xC;
					final double distance = Math.sqrt(dX * dX + dY2 + dZ2);
					if (distance > maxRad)
						continue;
					else
						targetIP.set(x - startX, tY, ip.get(x, y));
				}
			}
		}
		ImagePlus target = new ImagePlus("Sphere", targetStack);
		target.setCalibration(cal);
		target.setDisplayRange(imp.getDisplayRangeMin(),
				imp.getDisplayRangeMax());
		return target;
	}

	/**
	 * 
	 * @param imp
	 * @param cropFactor
	 * @param sphereDim
	 * @return
	 */
	public ImagePlus copyInnerCube(ImagePlus imp, double cropFactor,
			double[] sphereDim) {
		Calibration cal = imp.getCalibration();
		final double vW = cal.pixelWidth;
		final double vH = cal.pixelHeight;
		final double vD = cal.pixelDepth;
		final double xC = sphereDim[0];
		final double yC = sphereDim[1];
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
		final ImageStack sourceStack = imp.getStack();
		ImageStack targetStack = new ImageStack(roiWidth, roiHeight, roiDepth);
		final int endZ = startZ + roiDepth;
		final int endY = startY + roiHeight;
		final int endX = startX + roiWidth;
		for (int z = startZ; z < endZ; z++) {
			IJ.showProgress(z - startZ, roiDepth);
			IJ.showStatus("Copying largest enclosed cube");
			final int tZ = z - startZ + 1;
			targetStack.setPixels(
					Moments.getEmptyPixels(roiWidth, roiHeight,
							imp.getBitDepth()), tZ);
			targetStack.setSliceLabel(sourceStack.getShortSliceLabel(z), tZ);
			final ImageProcessor ip = sourceStack.getProcessor(z);
			ImageProcessor targetIP = targetStack.getProcessor(tZ);
			for (int y = startY; y < endY; y++) {
				final int tY = y - startY;
				for (int x = startX; x < endX; x++) {
					targetIP.set(x - startX, tY, ip.get(x, y));
				}
			}
		}
		ImagePlus target = new ImagePlus("Inner Cube", targetStack);
		target.setCalibration(cal);
		target.setDisplayRange(imp.getDisplayRangeMin(),
				imp.getDisplayRangeMax());
		return target;
	}

	/**
	 * 
	 * @param imp
	 * @param cropFactor
	 * @param sphereDim
	 * @return
	 */
	public ImagePlus copyOuterCube(ImagePlus imp, double cropFactor,
			double[] sphereDim) {
		Calibration cal = imp.getCalibration();
		final double vW = cal.pixelWidth;
		final double vH = cal.pixelHeight;
		final double vD = cal.pixelDepth;
		final double xC = sphereDim[0];
		final double yC = sphereDim[1];
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
		final ImageStack sourceStack = imp.getImageStack();
		ImageStack targetStack = new ImageStack(roiWidth, roiHeight, roiDepth);
		final int endZ = startZ + roiDepth;
		final int endY = startY + roiHeight;
		final int endX = startX + roiWidth;
		for (int z = startZ; z < endZ; z++) {
			final int tZ = z - startZ + 1;
			IJ.showProgress(z - startZ, roiDepth);
			IJ.showStatus("Copying smallest enclosing cube");
			targetStack.setPixels(
					Moments.getEmptyPixels(roiWidth, roiHeight,
							imp.getBitDepth()), tZ);
			targetStack.setSliceLabel(sourceStack.getShortSliceLabel(z), tZ);
			final ImageProcessor ip = sourceStack.getProcessor(z);
			ImageProcessor targetIP = targetStack.getProcessor(tZ);
			for (int y = startY; y < endY; y++) {
				final int tY = y - startY;
				for (int x = startX; x < endX; x++) {
					targetIP.set(x - startX, tY, ip.get(x, y));
				}
			}
		}
		ImagePlus target = new ImagePlus("Outer Cube", targetStack);
		target.setCalibration(cal);
		target.setDisplayRange(imp.getDisplayRangeMin(),
				imp.getDisplayRangeMax());
		return target;
	}

	/**
	 * Add series of circular ROIs to the ROI Manager based on the centre and
	 * radius of a sphere
	 * 
	 * @param imp
	 *            Needed for decalibration of calibrated (x,y,z) r values
	 * @param roiMan
	 *            Instance of the ROI Manager i.e. RoiManager.getInstance().
	 * @param sphereDim
	 *            calibrated centroid (x, y, z) and radius
	 * @param clearRois
	 * @throws IllegalArgumentException
	 *             if roiMan is null, rather than instantiating RoiManager.
	 */
	public static void addToRoiManager(ImagePlus imp, RoiManager roiMan,
			double[] sphereDim, boolean clearRois)
			throws IllegalArgumentException {
		if (roiMan == null)
			throw new IllegalArgumentException(
					"ROI Manager has not been instantiated");
		if (clearRois) {
			RoiMan.deleteAll(roiMan);
		}
		Calibration cal = imp.getCalibration();
		final double xs = sphereDim[0];
		final int xi = (int) (xs / cal.pixelWidth);
		final double ys = sphereDim[1];
		final int yi = (int) (ys / cal.pixelHeight);
		final double r = sphereDim[3];
		final int zc = (int) Math.round(sphereDim[2] / cal.pixelDepth);
		final int rz = (int) Math.round(r / cal.pixelDepth);
		final int zStart = Math.max(zc - rz, 1);
		final int zEnd = Math.min(zc + rz, imp.getImageStackSize());
		for (int z = zStart; z <= zEnd; z++) { // iterate through z slices
			final double zd = (zc - z) * cal.pixelDepth;
			final double rc = Math.sqrt(r * r - zd * zd);
			final int wi = (int) (rc / cal.pixelWidth);
			final int hi = (int) (rc / cal.pixelHeight);
			OvalRoi ellipse = new OvalRoi(xi - wi, yi - hi, wi * 2, hi * 2);
			ellipse.setPosition(z);
			roiMan.addRoi(ellipse);
		}
	}

	public boolean dialogItemChanged(GenericDialog gd, AWTEvent e) {
		if (!DialogModifier.allNumbersValid(gd.getNumericFields()))
			return false;
		Vector<?> checkboxes = gd.getCheckboxes();
		Vector<?> numbers = gd.getNumericFields();
		Checkbox box = (Checkbox) checkboxes.get(0);
		TextField num = (TextField) numbers.get(0);
		if (box.getState()) {
			num.setEnabled(true);
		} else {
			num.setEnabled(false);
		}
		DialogModifier.registerMacroValues(gd, gd.getComponents());
		return true;
	}
}