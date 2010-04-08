package org.doube.bonej;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.Plot;
import ij.macro.Interpreter;
import ij.measure.CurveFitter;
import ij.plugin.PlugIn;
import ij.util.Tools;

import java.util.ArrayList;

import org.doube.util.ImageCheck;
import org.doube.util.ResultInserter;

/**
 * <p>
 * FractalCount_.java,v 1.42 2005/05/30 07:52:59
 * </p>
 * <p>
 * Estimates the fractal dimension of 2D and 3D binary images. <br />
 * Free Software in the Public domain.
 * </p>
 * 
 * @author Per Christian Henden
 * @author Jens Bache-Wiig
 * @author Michael Doube (refactoring)
 * @see <p>
 *      <a
 *      href="http://www.pvv.org/~perchrh/imagej/fractal.html">http://www.pvv.
 *      org/~perchrh/imagej/fractal.html</a>
 *      </p>
 * 
 */

public class FractalBoxCounter implements PlugIn {

	boolean noGo = false;
	final int autoDiv = 4;

	// User-changeable defaults :
	boolean plotGraph = true;
	boolean verboseOutput = true;
	int threshold = 128;
	int maxBox = 24;
	int minBox = 6;
	double divBox = 1.2;
	int numOffsets = 1;
	boolean autoParam = true;

	// TODO split run method into more sensible methods
	public void run(String arg) {
		if (!ImageCheck.checkIJVersion())
			return;
		ImagePlus imp = IJ.getImage();
		if (null == imp) {
			IJ.noImage();
			return;
		}
		ImageCheck ic = new ImageCheck();
		if (!ic.isBinary(imp)) {
			IJ.showMessage("Fractal Count requires a binary image.");
			return;
		}
		if (noGo)
			return;
		try {
			// Fetch data
			final int width = imp.getWidth();
			final int height = imp.getHeight();
			final int depth = imp.getStackSize();

			if (autoParam) {
				maxBox = Math.max(width, Math.max(height, depth)) / autoDiv;
				if (verboseOutput) {
					IJ.log("Automatic max box size " + maxBox + " selected");
				}
			}

			// Create variables we need and set them
			long bestCount; // keep track of best count so far
			long count = 0; // current count
			ArrayList<Double> xList = new ArrayList<Double>();
			ArrayList<Double> yList = new ArrayList<Double>();
			int xPos, yPos, zPos, yPart;
			int xGrid, yGrid, zGrid;
			int xStart, yStart, zStart;
			int xEnd, yEnd, zEnd;

			// Start timer
			long startTime = System.currentTimeMillis();

			for (int boxSize = maxBox; boxSize >= minBox; boxSize /= divBox) {
				if (verboseOutput) {
					IJ.showStatus("Estimating dimension, box size: " + boxSize);
				}

				bestCount = Long.MAX_VALUE; // init count for this boxSize

				final int increment = Math.max(1, boxSize / numOffsets);

				for (int gridOffsetX = 0; (gridOffsetX < boxSize)
						&& (gridOffsetX < width); gridOffsetX += increment) {

					for (int gridOffsetY = 0; (gridOffsetY < boxSize)
							&& (gridOffsetY < height); gridOffsetY += increment) {

						for (int gridOffsetZ = 0; (gridOffsetZ < boxSize)
								&& (gridOffsetZ < depth); gridOffsetZ += increment) {

							count = 0;

							final int iMax = width + gridOffsetX;
							final int jMax = height + gridOffsetY;
							final int kMax = depth + gridOffsetZ;

							// Iterate over box-grid
							for (int i = 0; i <= iMax; i += boxSize) {
								xGrid = -gridOffsetX + i;
								for (int j = 0; j <= jMax; j += boxSize) {
									yGrid = -gridOffsetY + j;
									for (int k = 0; k <= kMax; k += boxSize) {
										zGrid = -gridOffsetZ + k;

										xStart = 0;
										if (xGrid < 0) {
											xStart = -xGrid;
										}
										if ((boxSize + xGrid) >= width) {
											xEnd = Math.min(width,
													(width - xGrid));
										} else {
											xEnd = boxSize;
										}

										yStart = 0;
										if (yGrid < 0) {
											yStart = -yGrid;
										}
										if ((boxSize + yGrid) >= height) {
											yEnd = Math.min(height,
													(height - yGrid));
										} else {
											yEnd = boxSize;
										}

										zStart = 0;
										if (zGrid < 0) {
											zStart = -zGrid;
										}
										if ((boxSize + zGrid) >= depth) {
											zEnd = Math.min(depth,
													(depth - zGrid));
										} else {
											zEnd = boxSize;
										}

										for (int x = xStart; x < xEnd; x++) {
											xPos = x + xGrid;

											for (int y = yStart; y < yEnd; y++) {
												yPos = (y + yGrid);
												yPart = yPos * width;

												for (int z = zStart; z < zEnd; z++) {
													zPos = z + zGrid;

													// If pixel inside region,
													// count it

													if ((0xff & ((byte[]) imp
															.getStack()
															.getPixels(zPos + 1))[xPos
															+ yPart]) >= threshold) {
														count++;
														z = x = y = boxSize; // stops
														// things
													}
												}
											}
										}
									}
								}
							}
							if (count < bestCount) {
								bestCount = count;
							}
						}
					}
				}
				xList.add(new Double(boxSize));
				yList.add(new Double(bestCount));
			}

			if (verboseOutput) {
				IJ.log("\nDuration: "
						+ (System.currentTimeMillis() - startTime) / 1000.0
						+ " seconds\n");
			}

			double[] boxSizes = new double[xList.size()];
			double[] boxCountSums = new double[yList.size()];
			for (int i = 0; i < boxSizes.length; i++) {
				boxSizes[i] = -Math.log((xList.get(i)).doubleValue());
				boxCountSums[i] = Math.log((yList.get(i)).doubleValue());
			}

			if (verboseOutput) {
				IJ.log("Used " + boxSizes.length
						+ " different box sizes, from " + maxBox + " to "
						+ minBox);
				IJ.log("with a reduction rate of " + divBox + " and "
						+ numOffsets + " translations of each box.");
			}

			if (boxSizes.length == 0) {
				IJ.error("\nError: No boxes!\n"
						+ "Make sure that starting and ending box size and"
						+ "reduction rate allow for at least one"
						+ " box size to exist!");
				return;
			}

			CurveFitter cf = new CurveFitter(boxSizes, boxCountSums);
			cf.doFit(CurveFitter.STRAIGHT_LINE);
			double[] p = cf.getParams();
			double RSq = cf.getRSquared();

			if (verboseOutput) {
				IJ.log(imp.getTitle() + ": Dimension estimate: "
						+ IJ.d2s(p[1], 4) + ": R²: " + RSq + ": Settings: "
						+ maxBox + ":" + minBox + ":" + divBox + ":"
						+ numOffsets);
			}

			ResultInserter ri = ResultInserter.getInstance();
			ri.setResultInRow(imp, "Fractal Dimension", p[1]);
			ri.setResultInRow(imp, "R²", RSq);
			ri.updateTable();

			if (plotGraph && !Interpreter.isBatchMode()) {
				drawGraph(p, boxSizes, boxCountSums);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	private void drawGraph(double[] params, double[] boxSizes,
			double[] boxCountSums) {

		final int samples = 100;
		float[] px = new float[samples];
		float[] py = new float[samples];
		double[] a = Tools.getMinMax(boxSizes);
		double xmin = a[0], xmax = a[1];

		a = Tools.getMinMax(boxCountSums);
		double ymin = a[0], ymax = a[1];
		final double inc = (xmax - xmin) / ((double) samples - 1);
		double tmp = xmin;

		for (int i = 0; i < samples; i++) {
			px[i] = (float) tmp;
			tmp += inc;
		}
		for (int i = 0; i < samples; i++) {
			py[i] = (float) CurveFitter.f(CurveFitter.STRAIGHT_LINE, params,
					px[i]);
		}
		a = Tools.getMinMax(py);
		ymin = Math.min(ymin, a[0]);
		ymax = Math.max(ymax, a[1]);
		Plot plot = new Plot("Plot", "-log(box size)", "log(box count)", px, py);
		ImagePlus plotImage = new ImagePlus("Plot", plot.getProcessor());
		plotImage.show();
		plot.setLimits(xmin, xmax * 0.9, ymin, ymax * 1.1);
		plot.addPoints(boxSizes, boxCountSums, Plot.CIRCLE);
		plot.addPoints(px, py, Plot.LINE);
		plot.addLabel(0.25, 0.25, "Slope: " + IJ.d2s(params[1], 4));
		return;
	}

	public void showAbout() {
		IJ.showMessage("About FractalCount..",
				"This plugin calculates the boxing dimension"
						+ " (an estimate of the fractal dimension) \n"
						+ "of 2D and 3D images.\n");
	}

}