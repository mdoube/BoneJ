package org.doube.skeleton;

/**
 * Skeletonize3D plugin for ImageJ(C).
 * Copyright (C) 2008 Ignacio Arganda-Carreras
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation (http://www.gnu.org/licenses/gpl.txt )
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 *
 */

import java.util.ArrayList;

import org.doube.util.ImageCheck;
import org.doube.util.UsageReporter;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.plugin.PlugIn;

/**
 * Main class. This class is a plugin for the ImageJ interface for 2D and 3D
 * thinning (skeletonization) of binary images (2D/3D).
 *
 * <p>
 * This work is an implementation by Ignacio Arganda-Carreras of the 3D thinning
 * algorithm from Lee et al.
 * </p>
 *
 * @see
 * 		<p>
 *      Lee, Kashyap and Chu (1994) Building skeleton models via 3-D medial
 *      surface/axis thinning algorithms. Computer Vision, Graphics, and Image
 *      Processing, 56(6):462–478
 *      <a href="http://dx.doi.org/10.1006/cgip.1994.1042" >doi:10.1006/cgip.
 *      1994.1042</a>.
 *      </p>
 *      <p>
 *      Based on the ITK version from Hanno Homann
 *      <a href="http://hdl.handle.net/1926/1292"> http://hdl.handle.net/1926/
 *      1292</a>
 *      </p>
 *      <p>
 *      <a href=
 *      "http://imagejdocu.tudor.lu/doku.php?id=plugin:morphology:skeletonize3d:start"
 *      >Skeletonize3D homepage</a>
 *      </p>
 *
 * @version 1.0 11/19/2008
 * @author Ignacio Arganda-Carreras <ignacio.arganda@uam.es>
 *
 */
public class Skeletonize3D implements PlugIn {

	public void run(final String run) {
		if (!ImageCheck.checkEnvironment())
			return;
		final ImagePlus imp = IJ.getImage();
		if (null == imp) {
			IJ.noImage();
			return;
		}
		final ImageCheck ic = new ImageCheck();
		if (!ic.isBinary(imp)) {
			IJ.error("Skeletonise 3D requires a binary image");
			return;
		}

		final ImagePlus skeleton = getSkeleton(imp);

		skeleton.show();
		if (imp.isInvertedLut() && !skeleton.isInvertedLut())
			IJ.run("Invert LUT");
		UsageReporter.reportEvent(this).send();
		return;
	}

	/**
	 * Gets a medial axis skeleton from a binary imp using a topology-preserving
	 * iterative algorithm
	 *
	 * @param imp
	 *            input image
	 * @return skeletonised image
	 */
	public ImagePlus getSkeleton(final ImagePlus imp) {
		final ImageStack inputImage = imp.getStack();

		// Prepare data
		final ImageStack outputImage = prepareData(inputImage);

		// Compute Thinning
		final ImagePlus imp2 = computeThinImage(outputImage);

		final ImageStack stack2 = imp2.getImageStack();

		// Convert image to binary 0-255
		for (int i = 1; i <= stack2.getSize(); i++)
			stack2.getProcessor(i).multiply(255);

		imp2.setCalibration(imp.getCalibration());
		imp2.setTitle("Skeleton of " + imp.getTitle());
		return imp2;
	}

	/*
	 * -----------------------------------------------------------------------
	 */
	/**
	 * Prepare data for computation. Copy the input image to the output image,
	 * changing from the input type to the output type.
	 *
	 * @param outputImage
	 *            output image stack
	 */
	private ImageStack prepareData(final ImageStack inputImage) {
		final int width = inputImage.getWidth();
		final int height = inputImage.getHeight();
		final int depth = inputImage.getSize();
		// IJ.write("Prepare Data: Copy input to output");
		IJ.showStatus("Prepare Data: Copy input to output ...");

		final ImageStack outputImage = new ImageStack(width, height, depth);
		// Copy the input to the output, changing all foreground pixels to
		// have value 1 in the process.
		for (int z = 0; z < depth; z++) {
			final byte[] pixels = (byte[]) inputImage.getPixels(z + 1);
			outputImage.setPixels(pixels.clone(), z + 1);
			outputImage.setSliceLabel(inputImage.getSliceLabel(z + 1), z + 1);
			for (int x = 0; x < width; x++)
				for (int y = 0; y < height; y++)
					if (((byte[]) inputImage.getPixels(z + 1))[x + y * width] != 0) {
						((byte[]) outputImage.getPixels(z + 1))[x + y * width] = 1;
					} else {
						((byte[]) outputImage.getPixels(z + 1))[x + y * width] = 0;
					}
		}
		// IJ.write("Prepare Data End");
		IJ.showStatus("Prepare Data End.");
		return outputImage;
	} /* end prepareData */

	/*
	 * -----------------------------------------------------------------------
	 */
	/**
	 * Post processing for computing thinning.
	 *
	 * @param outputImage
	 *            output image stack
	 */
	public ImagePlus computeThinImage(final ImageStack outputImage) {
		// IJ.write("Compute Thin Image Start");
		IJ.showStatus("Computing thin image ...");
		final int w = outputImage.getWidth();
		final int h = outputImage.getHeight();
		final int d = outputImage.getSize();

		final ArrayList<int[]> simpleBorderPoints = new ArrayList<int[]>();

		// Prepare Euler LUT [Lee94]
		final int eulerLUT[] = new int[256];
		fillEulerLUT(eulerLUT);

		int iter = 1;

		// Loop through the image several times until there is no change.
		int unchangedBorders = 0;
		while (unchangedBorders < 6) // loop until no change for all the six
										// border types
		{
			unchangedBorders = 0;
			for (int currentBorder = 1; currentBorder <= 6; currentBorder++) {
				IJ.showStatus("Thinning iteration " + iter + " (" + currentBorder + "/6 borders) ...");

				// Loop through the image.
				for (int z = 0; z < d; z++) {
					for (int y = 0; y < h; y++) {
						for (int x = 0; x < w; x++) {
							// check if point is foreground
							if (getPixel(outputImage, x, y, z, w, h, d) != 1) {
								continue; // current point is already background
							}
							// check 6-neighbors if point is a border point of
							// type currentBorder
							boolean isBorderPoint = false;
							// North
							if (currentBorder == 1 && N(outputImage, x, y, z, w, h, d) <= 0)
								isBorderPoint = true;
							// South
							if (currentBorder == 2 && S(outputImage, x, y, z, w, h, d) <= 0)
								isBorderPoint = true;
							// East
							if (currentBorder == 3 && E(outputImage, x, y, z, w, h, d) <= 0)
								isBorderPoint = true;
							// West
							if (currentBorder == 4 && W(outputImage, x, y, z, w, h, d) <= 0)
								isBorderPoint = true;
							// Up
							if (currentBorder == 5 && U(outputImage, x, y, z, w, h, d) <= 0)
								isBorderPoint = true;
							// Bottom
							if (currentBorder == 6 && B(outputImage, x, y, z, w, h, d) <= 0)
								isBorderPoint = true;

							if (!isBorderPoint) {
								continue; // current point is not deletable
							}

							// check if point is the end of an arc
							int numberOfNeighbors = -1; // -1 and not 0 because
														// the center pixel will
														// be counted as well
							final byte[] neighbor = getNeighborhood(outputImage, x, y, z, w, h, d);
							for (int i = 0; i < 27; i++) // i = 0..26
							{
								if (neighbor[i] == 1)
									numberOfNeighbors++;
							}

							if (numberOfNeighbors == 1) {
								continue; // current point is not deletable
							}

							// Check if point is Euler invariant
							if (!isEulerInvariant(getNeighborhood(outputImage, x, y, z, w, h, d), eulerLUT)) {
								continue; // current point is not deletable
							}
							// Check if point is simple (deletion does not
							// change connectivity in the 3x3x3 neighborhood)
							if (!isSimplePoint(getNeighborhood(outputImage, x, y, z, w, h, d))) {
								continue; // current point is not deletable
							}
							// add all simple border points to a list for
							// sequential re-checking
							final int[] index = new int[3];
							index[0] = x;
							index[1] = y;
							index[2] = z;
							simpleBorderPoints.add(index);
						}
					}
					IJ.showProgress(z, d);
				}

				// sequential re-checking to preserve connectivity when
				// deleting in a parallel way
				boolean noChange = true;
				int[] index = null;
				for (int i = 0; i < simpleBorderPoints.size(); i++) {
					index = simpleBorderPoints.get(i);
					// 1. Set simple border point to 0
					setPixel(outputImage, index[0], index[1], index[2], w, h, d, (byte) 0);

					// 2. Check if neighborhood is still connected
					if (!isSimplePoint(getNeighborhood(outputImage, index[0], index[1], index[2], w, h, d))) {
						// we cannot delete current point, so reset
						setPixel(outputImage, index[0], index[1], index[2], w, h, d, (byte) 1);
					} else {
						noChange = false;
					}
				}
				if (noChange)
					unchangedBorders++;

				// IJ.write("# simple border points = " +
				// simpleBorderPoints.size());
				/*
				 * if(index != null) IJ.write("# last point = [" + index[0] +
				 * ", " + index[1] + "," + index[2]+ "]"); else IJ.write(
				 * "# last point = [0, 0, 0]");
				 */
				simpleBorderPoints.clear();
				// IJ.write("# simple border points = " +
				// simpleBorderPoints.size() + "\n");

			} // end currentBorder for loop

			// Progress bar iterations
			iter++;
		}

		// IJ.write("Compute Thin Image End");
		IJ.showStatus("Computed thin image.");
		final ImagePlus impOut = new ImagePlus();
		impOut.setStack(outputImage);
		return impOut;
	} /* end computeThinImage */

	/*
	 * -----------------------------------------------------------------------
	 */
	/**
	 * Get neighborhood of a pixel in a 3D image (0 border conditions)
	 *
	 * @param image
	 *            3D image (ImageStack)
	 * @param x
	 *            x- coordinate
	 * @param y
	 *            y- coordinate
	 * @param z
	 *            z- coordinate (in image stacks the indexes start at 1)
	 * @return corresponding 27-pixels neighborhood (0 if out of image)
	 */
	private byte[] getNeighborhood(final ImageStack image, final int x, final int y, final int z, final int w,
			final int h, final int d) {
		final byte[] neighborhood = new byte[27];

		neighborhood[0] = getPixel(image, x - 1, y - 1, z - 1, w, h, d);
		neighborhood[1] = getPixel(image, x, y - 1, z - 1, w, h, d);
		neighborhood[2] = getPixel(image, x + 1, y - 1, z - 1, w, h, d);

		neighborhood[3] = getPixel(image, x - 1, y, z - 1, w, h, d);
		neighborhood[4] = getPixel(image, x, y, z - 1, w, h, d);
		neighborhood[5] = getPixel(image, x + 1, y, z - 1, w, h, d);

		neighborhood[6] = getPixel(image, x - 1, y + 1, z - 1, w, h, d);
		neighborhood[7] = getPixel(image, x, y + 1, z - 1, w, h, d);
		neighborhood[8] = getPixel(image, x + 1, y + 1, z - 1, w, h, d);

		neighborhood[9] = getPixel(image, x - 1, y - 1, z, w, h, d);
		neighborhood[10] = getPixel(image, x, y - 1, z, w, h, d);
		neighborhood[11] = getPixel(image, x + 1, y - 1, z, w, h, d);

		neighborhood[12] = getPixel(image, x - 1, y, z, w, h, d);
		neighborhood[13] = getPixel(image, x, y, z, w, h, d);
		neighborhood[14] = getPixel(image, x + 1, y, z, w, h, d);

		neighborhood[15] = getPixel(image, x - 1, y + 1, z, w, h, d);
		neighborhood[16] = getPixel(image, x, y + 1, z, w, h, d);
		neighborhood[17] = getPixel(image, x + 1, y + 1, z, w, h, d);

		neighborhood[18] = getPixel(image, x - 1, y - 1, z + 1, w, h, d);
		neighborhood[19] = getPixel(image, x, y - 1, z + 1, w, h, d);
		neighborhood[20] = getPixel(image, x + 1, y - 1, z + 1, w, h, d);

		neighborhood[21] = getPixel(image, x - 1, y, z + 1, w, h, d);
		neighborhood[22] = getPixel(image, x, y, z + 1, w, h, d);
		neighborhood[23] = getPixel(image, x + 1, y, z + 1, w, h, d);

		neighborhood[24] = getPixel(image, x - 1, y + 1, z + 1, w, h, d);
		neighborhood[25] = getPixel(image, x, y + 1, z + 1, w, h, d);
		neighborhood[26] = getPixel(image, x + 1, y + 1, z + 1, w, h, d);

		return neighborhood;
	} /* end getNeighborhood */

	/*
	 * -----------------------------------------------------------------------
	 */
	/**
	 * Get pixel in 3D image (0 border conditions)
	 *
	 * @param image
	 *            3D image
	 * @param x
	 *            x- coordinate
	 * @param y
	 *            y- coordinate
	 * @param z
	 *            z- coordinate (in image stacks the indexes start at 1)
	 * @return corresponding pixel (0 if out of image)
	 */
	private byte getPixel(final ImageStack image, final int x, final int y, final int z, final int w, final int h,
			final int d) {
		if (x >= 0 && x < w && y >= 0 && y < h && z >= 0 && z < d)
			return ((byte[]) image.getPixels(z + 1))[x + y * w];
		else
			return 0;
	} /* end getPixel */

	/*
	 * -----------------------------------------------------------------------
	 */
	/**
	 * Set pixel in 3D image
	 *
	 * @param image
	 *            3D image
	 * @param x
	 *            x- coordinate
	 * @param y
	 *            y- coordinate
	 * @param z
	 *            z- coordinate (in image stacks the indexes start at 1)
	 * @param value
	 *            pixel value
	 */
	private void setPixel(final ImageStack image, final int x, final int y, final int z, final int w, final int h,
			final int d, final byte value) {
		if (x >= 0 && x < w && y >= 0 && y < h && z >= 0 && z < d)
			((byte[]) image.getPixels(z + 1))[x + y * w] = value;
	} /* end getPixel */

	/*
	 * -----------------------------------------------------------------------
	 */
	/**
	 * North neighborhood (0 border conditions)
	 *
	 * @param image
	 *            3D image
	 * @param x
	 *            x- coordinate
	 * @param y
	 *            y- coordinate
	 * @param z
	 *            z- coordinate (in image stacks the indexes start at 1)
	 * @return corresponding north pixel
	 */
	private byte N(final ImageStack image, final int x, final int y, final int z, final int w, final int h,
			final int d) {
		return getPixel(image, x, y - 1, z, w, h, d);
	} /* end N */

	/*
	 * -----------------------------------------------------------------------
	 */
	/**
	 * South neighborhood (0 border conditions)
	 *
	 * @param image
	 *            3D image
	 * @param x
	 *            x- coordinate
	 * @param y
	 *            y- coordinate
	 * @param z
	 *            z- coordinate (in image stacks the indexes start at 1)
	 * @return corresponding south pixel
	 */
	private byte S(final ImageStack image, final int x, final int y, final int z, final int w, final int h,
			final int d) {
		return getPixel(image, x, y + 1, z, w, h, d);
	} /* end S */

	/*
	 * -----------------------------------------------------------------------
	 */
	/**
	 * East neighborhood (0 border conditions)
	 *
	 * @param image
	 *            3D image
	 * @param x
	 *            x- coordinate
	 * @param y
	 *            y- coordinate
	 * @param z
	 *            z- coordinate (in image stacks the indexes start at 1)
	 * @return corresponding east pixel
	 */
	private byte E(final ImageStack image, final int x, final int y, final int z, final int w, final int h,
			final int d) {
		return getPixel(image, x + 1, y, z, w, h, d);
	} /* end E */

	/*
	 * -----------------------------------------------------------------------
	 */
	/**
	 * West neighborhood (0 border conditions)
	 *
	 * @param image
	 *            3D image
	 * @param x
	 *            x- coordinate
	 * @param y
	 *            y- coordinate
	 * @param z
	 *            z- coordinate (in image stacks the indexes start at 1)
	 * @return corresponding west pixel
	 */
	private byte W(final ImageStack image, final int x, final int y, final int z, final int w, final int h,
			final int d) {
		return getPixel(image, x - 1, y, z, w, h, d);
	} /* end W */

	/*
	 * -----------------------------------------------------------------------
	 */
	/**
	 * Up neighborhood (0 border conditions)
	 *
	 * @param image
	 *            3D image
	 * @param x
	 *            x- coordinate
	 * @param y
	 *            y- coordinate
	 * @param z
	 *            z- coordinate (in image stacks the indexes start at 1)
	 * @return corresponding up pixel
	 */
	private byte U(final ImageStack image, final int x, final int y, final int z, final int w, final int h,
			final int d) {
		return getPixel(image, x, y, z + 1, w, h, d);
	} /* end U */

	/*
	 * -----------------------------------------------------------------------
	 */
	/**
	 * Bottom neighborhood (0 border conditions)
	 *
	 * @param image
	 *            3D image
	 * @param x
	 *            x- coordinate
	 * @param y
	 *            y- coordinate
	 * @param z
	 *            z- coordinate (in image stacks the indexes start at 1)
	 * @return corresponding bottom pixel
	 */
	private byte B(final ImageStack image, final int x, final int y, final int z, final int w, final int h,
			final int d) {
		return getPixel(image, x, y, z - 1, w, h, d);
	} /* end N */

	/*
	 * -----------------------------------------------------------------------
	 */
	/**
	 * Fill Euler LUT
	 *
	 * @param LUT
	 *            Euler LUT
	 */
	private void fillEulerLUT(final int[] LUT) {
		LUT[1] = 1;
		LUT[3] = -1;
		LUT[5] = -1;
		LUT[7] = 1;
		LUT[9] = -3;
		LUT[11] = -1;
		LUT[13] = -1;
		LUT[15] = 1;
		LUT[17] = -1;
		LUT[19] = 1;
		LUT[21] = 1;
		LUT[23] = -1;
		LUT[25] = 3;
		LUT[27] = 1;
		LUT[29] = 1;
		LUT[31] = -1;
		LUT[33] = -3;
		LUT[35] = -1;
		LUT[37] = 3;
		LUT[39] = 1;
		LUT[41] = 1;
		LUT[43] = -1;
		LUT[45] = 3;
		LUT[47] = 1;
		LUT[49] = -1;
		LUT[51] = 1;

		LUT[53] = 1;
		LUT[55] = -1;
		LUT[57] = 3;
		LUT[59] = 1;
		LUT[61] = 1;
		LUT[63] = -1;
		LUT[65] = -3;
		LUT[67] = 3;
		LUT[69] = -1;
		LUT[71] = 1;
		LUT[73] = 1;
		LUT[75] = 3;
		LUT[77] = -1;
		LUT[79] = 1;
		LUT[81] = -1;
		LUT[83] = 1;
		LUT[85] = 1;
		LUT[87] = -1;
		LUT[89] = 3;
		LUT[91] = 1;
		LUT[93] = 1;
		LUT[95] = -1;
		LUT[97] = 1;
		LUT[99] = 3;
		LUT[101] = 3;
		LUT[103] = 1;

		LUT[105] = 5;
		LUT[107] = 3;
		LUT[109] = 3;
		LUT[111] = 1;
		LUT[113] = -1;
		LUT[115] = 1;
		LUT[117] = 1;
		LUT[119] = -1;
		LUT[121] = 3;
		LUT[123] = 1;
		LUT[125] = 1;
		LUT[127] = -1;
		LUT[129] = -7;
		LUT[131] = -1;
		LUT[133] = -1;
		LUT[135] = 1;
		LUT[137] = -3;
		LUT[139] = -1;
		LUT[141] = -1;
		LUT[143] = 1;
		LUT[145] = -1;
		LUT[147] = 1;
		LUT[149] = 1;
		LUT[151] = -1;
		LUT[153] = 3;
		LUT[155] = 1;

		LUT[157] = 1;
		LUT[159] = -1;
		LUT[161] = -3;
		LUT[163] = -1;
		LUT[165] = 3;
		LUT[167] = 1;
		LUT[169] = 1;
		LUT[171] = -1;
		LUT[173] = 3;
		LUT[175] = 1;
		LUT[177] = -1;
		LUT[179] = 1;
		LUT[181] = 1;
		LUT[183] = -1;
		LUT[185] = 3;
		LUT[187] = 1;
		LUT[189] = 1;
		LUT[191] = -1;
		LUT[193] = -3;
		LUT[195] = 3;
		LUT[197] = -1;
		LUT[199] = 1;
		LUT[201] = 1;
		LUT[203] = 3;
		LUT[205] = -1;
		LUT[207] = 1;

		LUT[209] = -1;
		LUT[211] = 1;
		LUT[213] = 1;
		LUT[215] = -1;
		LUT[217] = 3;
		LUT[219] = 1;
		LUT[221] = 1;
		LUT[223] = -1;
		LUT[225] = 1;
		LUT[227] = 3;
		LUT[229] = 3;
		LUT[231] = 1;
		LUT[233] = 5;
		LUT[235] = 3;
		LUT[237] = 3;
		LUT[239] = 1;
		LUT[241] = -1;
		LUT[243] = 1;
		LUT[245] = 1;
		LUT[247] = -1;
		LUT[249] = 3;
		LUT[251] = 1;
		LUT[253] = 1;
		LUT[255] = -1;
	}

	/**
	 * Check if a point is Euler invariant
	 *
	 * @param neighbors
	 *            neighbor pixels of the point
	 * @param LUT
	 *            Euler LUT
	 * @return true or false if the point is Euler invariant or not
	 */
	private boolean isEulerInvariant(final byte[] neighbors, final int[] LUT) {
		// Calculate Euler characteristic for each octant and sum up
		int eulerChar = 0;
		char n;
		// Octant SWU
		n = 1;
		if (neighbors[24] == 1)
			n |= 128;
		if (neighbors[25] == 1)
			n |= 64;
		if (neighbors[15] == 1)
			n |= 32;
		if (neighbors[16] == 1)
			n |= 16;
		if (neighbors[21] == 1)
			n |= 8;
		if (neighbors[22] == 1)
			n |= 4;
		if (neighbors[12] == 1)
			n |= 2;
		eulerChar += LUT[n];
		// Octant SEU
		n = 1;
		if (neighbors[26] == 1)
			n |= 128;
		if (neighbors[23] == 1)
			n |= 64;
		if (neighbors[17] == 1)
			n |= 32;
		if (neighbors[14] == 1)
			n |= 16;
		if (neighbors[25] == 1)
			n |= 8;
		if (neighbors[22] == 1)
			n |= 4;
		if (neighbors[16] == 1)
			n |= 2;
		eulerChar += LUT[n];
		// Octant NWU
		n = 1;
		if (neighbors[18] == 1)
			n |= 128;
		if (neighbors[21] == 1)
			n |= 64;
		if (neighbors[9] == 1)
			n |= 32;
		if (neighbors[12] == 1)
			n |= 16;
		if (neighbors[19] == 1)
			n |= 8;
		if (neighbors[22] == 1)
			n |= 4;
		if (neighbors[10] == 1)
			n |= 2;
		eulerChar += LUT[n];
		// Octant NEU
		n = 1;
		if (neighbors[20] == 1)
			n |= 128;
		if (neighbors[23] == 1)
			n |= 64;
		if (neighbors[19] == 1)
			n |= 32;
		if (neighbors[22] == 1)
			n |= 16;
		if (neighbors[11] == 1)
			n |= 8;
		if (neighbors[14] == 1)
			n |= 4;
		if (neighbors[10] == 1)
			n |= 2;
		eulerChar += LUT[n];
		// Octant SWB
		n = 1;
		if (neighbors[6] == 1)
			n |= 128;
		if (neighbors[15] == 1)
			n |= 64;
		if (neighbors[7] == 1)
			n |= 32;
		if (neighbors[16] == 1)
			n |= 16;
		if (neighbors[3] == 1)
			n |= 8;
		if (neighbors[12] == 1)
			n |= 4;
		if (neighbors[4] == 1)
			n |= 2;
		eulerChar += LUT[n];
		// Octant SEB
		n = 1;
		if (neighbors[8] == 1)
			n |= 128;
		if (neighbors[7] == 1)
			n |= 64;
		if (neighbors[17] == 1)
			n |= 32;
		if (neighbors[16] == 1)
			n |= 16;
		if (neighbors[5] == 1)
			n |= 8;
		if (neighbors[4] == 1)
			n |= 4;
		if (neighbors[14] == 1)
			n |= 2;
		eulerChar += LUT[n];
		// Octant NWB
		n = 1;
		if (neighbors[0] == 1)
			n |= 128;
		if (neighbors[9] == 1)
			n |= 64;
		if (neighbors[3] == 1)
			n |= 32;
		if (neighbors[12] == 1)
			n |= 16;
		if (neighbors[1] == 1)
			n |= 8;
		if (neighbors[10] == 1)
			n |= 4;
		if (neighbors[4] == 1)
			n |= 2;
		eulerChar += LUT[n];
		// Octant NEB
		n = 1;
		if (neighbors[2] == 1)
			n |= 128;
		if (neighbors[1] == 1)
			n |= 64;
		if (neighbors[11] == 1)
			n |= 32;
		if (neighbors[10] == 1)
			n |= 16;
		if (neighbors[5] == 1)
			n |= 8;
		if (neighbors[4] == 1)
			n |= 4;
		if (neighbors[14] == 1)
			n |= 2;
		eulerChar += LUT[n];
		if (eulerChar == 0)
			return true;
		else
			return false;
	}

	/*
	 * -----------------------------------------------------------------------
	 */
	/**
	 * Check if current point is a Simple Point. This method is named
	 * 'N(v)_labeling' in [Lee94]. Outputs the number of connected objects in a
	 * neighborhood of a point after this point would have been removed.
	 *
	 * @param neighbors
	 *            neighbor pixels of the point
	 * @return true or false if the point is simple or not
	 */
	private boolean isSimplePoint(final byte[] neighbors) {
		// copy neighbors for labeling
		final int cube[] = new int[26];
		int i = 0;
		for (i = 0; i < 13; i++) // i = 0..12 -> cube[0..12]
			cube[i] = neighbors[i];
		// i != 13 : ignore center pixel when counting (see [Lee94])
		for (i = 14; i < 27; i++) // i = 14..26 -> cube[13..25]
			cube[i - 1] = neighbors[i];
		// set initial label
		int label = 2;
		// for all points in the neighborhood
		for (i = 0; i < 26; i++) {
			if (cube[i] == 1) // voxel has not been labelled yet
			{
				// start recursion with any octant that contains the point i
				switch (i) {
				case 0:
				case 1:
				case 3:
				case 4:
				case 9:
				case 10:
				case 12:
					octreeLabeling(1, label, cube);
					break;
				case 2:
				case 5:
				case 11:
				case 13:
					octreeLabeling(2, label, cube);
					break;
				case 6:
				case 7:
				case 14:
				case 15:
					octreeLabeling(3, label, cube);
					break;
				case 8:
				case 16:
					octreeLabeling(4, label, cube);
					break;
				case 17:
				case 18:
				case 20:
				case 21:
					octreeLabeling(5, label, cube);
					break;
				case 19:
				case 22:
					octreeLabeling(6, label, cube);
					break;
				case 23:
				case 24:
					octreeLabeling(7, label, cube);
					break;
				case 25:
					octreeLabeling(8, label, cube);
					break;
				}
				label++;
				if (label - 2 >= 2) {
					return false;
				}
			}
		}
		// return label-2; in [Lee94] if the number of connected components
		// would be needed
		return true;
	}

	/*
	 * -----------------------------------------------------------------------
	 */
	/**
	 * This is a recursive method that calculates the number of connected
	 * components in the 3D neighborhood after the center pixel would have been
	 * removed.
	 *
	 * @param octant
	 * @param label
	 * @param cube
	 */
	private void octreeLabeling(final int octant, final int label, final int[] cube) {
		// check if there are points in the octant with value 1
		if (octant == 1) {
			// set points in this octant to current label
			// and recursive labeling of adjacent octants
			if (cube[0] == 1)
				cube[0] = label;
			if (cube[1] == 1) {
				cube[1] = label;
				octreeLabeling(2, label, cube);
			}
			if (cube[3] == 1) {
				cube[3] = label;
				octreeLabeling(3, label, cube);
			}
			if (cube[4] == 1) {
				cube[4] = label;
				octreeLabeling(2, label, cube);
				octreeLabeling(3, label, cube);
				octreeLabeling(4, label, cube);
			}
			if (cube[9] == 1) {
				cube[9] = label;
				octreeLabeling(5, label, cube);
			}
			if (cube[10] == 1) {
				cube[10] = label;
				octreeLabeling(2, label, cube);
				octreeLabeling(5, label, cube);
				octreeLabeling(6, label, cube);
			}
			if (cube[12] == 1) {
				cube[12] = label;
				octreeLabeling(3, label, cube);
				octreeLabeling(5, label, cube);
				octreeLabeling(7, label, cube);
			}
		}
		if (octant == 2) {
			if (cube[1] == 1) {
				cube[1] = label;
				octreeLabeling(1, label, cube);
			}
			if (cube[4] == 1) {
				cube[4] = label;
				octreeLabeling(1, label, cube);
				octreeLabeling(3, label, cube);
				octreeLabeling(4, label, cube);
			}
			if (cube[10] == 1) {
				cube[10] = label;
				octreeLabeling(1, label, cube);
				octreeLabeling(5, label, cube);
				octreeLabeling(6, label, cube);
			}
			if (cube[2] == 1)
				cube[2] = label;
			if (cube[5] == 1) {
				cube[5] = label;
				octreeLabeling(4, label, cube);
			}
			if (cube[11] == 1) {
				cube[11] = label;
				octreeLabeling(6, label, cube);
			}
			if (cube[13] == 1) {
				cube[13] = label;
				octreeLabeling(4, label, cube);
				octreeLabeling(6, label, cube);
				octreeLabeling(8, label, cube);
			}
		}
		if (octant == 3) {
			if (cube[3] == 1) {
				cube[3] = label;
				octreeLabeling(1, label, cube);
			}
			if (cube[4] == 1) {
				cube[4] = label;
				octreeLabeling(1, label, cube);
				octreeLabeling(2, label, cube);
				octreeLabeling(4, label, cube);
			}
			if (cube[12] == 1) {
				cube[12] = label;
				octreeLabeling(1, label, cube);
				octreeLabeling(5, label, cube);
				octreeLabeling(7, label, cube);
			}
			if (cube[6] == 1)
				cube[6] = label;
			if (cube[7] == 1) {
				cube[7] = label;
				octreeLabeling(4, label, cube);
			}
			if (cube[14] == 1) {
				cube[14] = label;
				octreeLabeling(7, label, cube);
			}
			if (cube[15] == 1) {
				cube[15] = label;
				octreeLabeling(4, label, cube);
				octreeLabeling(7, label, cube);
				octreeLabeling(8, label, cube);
			}
		}
		if (octant == 4) {
			if (cube[4] == 1) {
				cube[4] = label;
				octreeLabeling(1, label, cube);
				octreeLabeling(2, label, cube);
				octreeLabeling(3, label, cube);
			}
			if (cube[5] == 1) {
				cube[5] = label;
				octreeLabeling(2, label, cube);
			}
			if (cube[13] == 1) {
				cube[13] = label;
				octreeLabeling(2, label, cube);
				octreeLabeling(6, label, cube);
				octreeLabeling(8, label, cube);
			}
			if (cube[7] == 1) {
				cube[7] = label;
				octreeLabeling(3, label, cube);
			}
			if (cube[15] == 1) {
				cube[15] = label;
				octreeLabeling(3, label, cube);
				octreeLabeling(7, label, cube);
				octreeLabeling(8, label, cube);
			}
			if (cube[8] == 1)
				cube[8] = label;
			if (cube[16] == 1) {
				cube[16] = label;
				octreeLabeling(8, label, cube);
			}
		}
		if (octant == 5) {
			if (cube[9] == 1) {
				cube[9] = label;
				octreeLabeling(1, label, cube);
			}
			if (cube[10] == 1) {
				cube[10] = label;
				octreeLabeling(1, label, cube);
				octreeLabeling(2, label, cube);
				octreeLabeling(6, label, cube);
			}
			if (cube[12] == 1) {
				cube[12] = label;
				octreeLabeling(1, label, cube);
				octreeLabeling(3, label, cube);
				octreeLabeling(7, label, cube);
			}
			if (cube[17] == 1)
				cube[17] = label;
			if (cube[18] == 1) {
				cube[18] = label;
				octreeLabeling(6, label, cube);
			}
			if (cube[20] == 1) {
				cube[20] = label;
				octreeLabeling(7, label, cube);
			}
			if (cube[21] == 1) {
				cube[21] = label;
				octreeLabeling(6, label, cube);
				octreeLabeling(7, label, cube);
				octreeLabeling(8, label, cube);
			}
		}
		if (octant == 6) {
			if (cube[10] == 1) {
				cube[10] = label;
				octreeLabeling(1, label, cube);
				octreeLabeling(2, label, cube);
				octreeLabeling(5, label, cube);
			}
			if (cube[11] == 1) {
				cube[11] = label;
				octreeLabeling(2, label, cube);
			}
			if (cube[13] == 1) {
				cube[13] = label;
				octreeLabeling(2, label, cube);
				octreeLabeling(4, label, cube);
				octreeLabeling(8, label, cube);
			}
			if (cube[18] == 1) {
				cube[18] = label;
				octreeLabeling(5, label, cube);
			}
			if (cube[21] == 1) {
				cube[21] = label;
				octreeLabeling(5, label, cube);
				octreeLabeling(7, label, cube);
				octreeLabeling(8, label, cube);
			}
			if (cube[19] == 1)
				cube[19] = label;
			if (cube[22] == 1) {
				cube[22] = label;
				octreeLabeling(8, label, cube);
			}
		}
		if (octant == 7) {
			if (cube[12] == 1) {
				cube[12] = label;
				octreeLabeling(1, label, cube);
				octreeLabeling(3, label, cube);
				octreeLabeling(5, label, cube);
			}
			if (cube[14] == 1) {
				cube[14] = label;
				octreeLabeling(3, label, cube);
			}
			if (cube[15] == 1) {
				cube[15] = label;
				octreeLabeling(3, label, cube);
				octreeLabeling(4, label, cube);
				octreeLabeling(8, label, cube);
			}
			if (cube[20] == 1) {
				cube[20] = label;
				octreeLabeling(5, label, cube);
			}
			if (cube[21] == 1) {
				cube[21] = label;
				octreeLabeling(5, label, cube);
				octreeLabeling(6, label, cube);
				octreeLabeling(8, label, cube);
			}
			if (cube[23] == 1)
				cube[23] = label;
			if (cube[24] == 1) {
				cube[24] = label;
				octreeLabeling(8, label, cube);
			}
		}
		if (octant == 8) {
			if (cube[13] == 1) {
				cube[13] = label;
				octreeLabeling(2, label, cube);
				octreeLabeling(4, label, cube);
				octreeLabeling(6, label, cube);
			}
			if (cube[15] == 1) {
				cube[15] = label;
				octreeLabeling(3, label, cube);
				octreeLabeling(4, label, cube);
				octreeLabeling(7, label, cube);
			}
			if (cube[16] == 1) {
				cube[16] = label;
				octreeLabeling(4, label, cube);
			}
			if (cube[21] == 1) {
				cube[21] = label;
				octreeLabeling(5, label, cube);
				octreeLabeling(6, label, cube);
				octreeLabeling(7, label, cube);
			}
			if (cube[22] == 1) {
				cube[22] = label;
				octreeLabeling(6, label, cube);
			}
			if (cube[24] == 1) {
				cube[24] = label;
				octreeLabeling(7, label, cube);
			}
			if (cube[25] == 1)
				cube[25] = label;
		}

	}

	/*
	 * -----------------------------------------------------------------------
	 */
	/**
	 * Show plug-in information.
	 *
	 */
	void showAbout() {
		IJ.showMessage("About Skeletonize3D...",
				"This plug-in filter produces 3D thinning (skeletonization) of binary 3D images.\n");
	} /* end showAbout */
	/*
	 * -----------------------------------------------------------------------
	 */

} /* end Skeletonize3D_ */
