package org.doube.skeleton;

import org.doube.util.ImageCheck;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.plugin.PlugIn;

public class ShapeSkeletoniser implements PlugIn {

	private static final byte BLACK = (byte) 255;
	private static final byte WHITE = (byte) 0;

	public void run(String arg) {
		if (!ImageCheck.checkIJVersion())
			return;
		ImagePlus imp = IJ.getImage();
		ImageCheck ic = new ImageCheck();
		if (!ic.isBinary(imp)) {
			IJ.error("Shape Skeletoniser requires a binary image.");
			return;
		}
	}

	/**
	 * Determine if a point is a (26, 6) simple point. For <i>p</i> to be a
	 * simple point, 4 conditions must be met:
	 * <ol>
	 * <li><i>p</i> has at least one black 26-neighbour</li>
	 * <li><i>p</i> has at least one white 6-neighbour</li>
	 * <li>The set of black 26-neighbours of <i>p</i> is 26-connected</li>
	 * <li>The set of white 6-neighbours of <i>p</i> is 6-connected in the set
	 * of white 18-neighbours of <i>p</i></li>
	 * </ol>
	 * 
	 * @param stack
	 * @param x
	 * @param y
	 * @param z
	 * @param w
	 * @param h
	 * @param d
	 * @return
	 */
	private boolean isSimplePoint(ImageStack stack, int x, int y, int z, int w,
			int h, int d) {
		byte[] neighbours = getNeighborhood(stack, x, y, z, w, h, d);

		// Condition 2. check 6 neighbourhood for white pixels
		boolean hasWhite6Neighbour = false;
		for (int i = 0; i < 6; i++) {
			if (neighbours[sixNeighbours[i]] == WHITE) {
				hasWhite6Neighbour = true;
				break;
			}
		}
		if (!hasWhite6Neighbour)
			return false;

		// Condition 1. check 26 neighbourhood for black pixels
		boolean hasBlack26Neighbour = false;
		for (int i = 0; i < 27; i++) {
			if (i == 13)
				continue;
			if (neighbours[i] == BLACK) {
				hasBlack26Neighbour = true;
				break;
			}
		}
		if (!hasBlack26Neighbour)
			return false;

		// Condition 3. check that black 26-neighbours are part of the same
		// 26-connected set
		if (!isBlack26ConnectedSet(neighbours))
			return false;

		// Condition 4.
		if (!isWhite6ConnectedIn18Set(neighbours))
			return false;

		// all conditions are true so p is a simple point
		return true;
	}

	/**
	 * Check that the white elements in the 6 neighbourhood are 6-connected by
	 * the white elements in the 18 neighbourhood (face-to-face around edges)
	 * 
	 * @param neighbours
	 * @return
	 */
	private boolean isWhite6ConnectedIn18Set(byte[] neighbours) {
		// count white 6-points
		int nwp = 0;
		for (int i = 0; i < 6; i++) {
			if (neighbours[sixNeighbours[i]] == WHITE) {
				nwp++;
			}
		}

		if (nwp == 1)
			return true;

		// else check for connection of the 2 white 6-points in the white
		// 18-neighbours

		// TODO Auto-generated method stub
		return false;
	}

	/**
	 * Check the black elements of the 26 neighbourhood and return true if they
	 * form a single 26-connected set
	 * 
	 * @param neighbours
	 * @return
	 */
	private boolean isBlack26ConnectedSet(byte[] neighbours) {
		// TODO Auto-generated method stub
		return false;
	}

	/** LUT to find the 6 neighbours in a 27 neighbourhood array */
	private static final int[] sixNeighbours = { 4, 10, 12, 14, 16, 22 };

	// -------------------------------------------------------------------//
	// Utility methods from Ignacio Arganda Carreras' Skeletonize3D
	// If classes merge, following methods are safe to exclude

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
	 * @param w
	 *            image width
	 * @param h
	 *            image height
	 * @param d
	 *            image depth
	 * @return corresponding pixel (0 if out of image)
	 */
	private byte getPixel(ImageStack image, int x, int y, int z, int w, int h,
			int d) {
		if (x >= 0 && x < w && y >= 0 && y < h && z >= 0 && z < d)
			return ((byte[]) image.getPixels(z + 1))[x + y * w];
		else
			return 0;
	} /* end getPixel */

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
	 * @param w
	 *            image width
	 * @param h
	 *            image height
	 * @param d
	 *            image depth
	 * @return corresponding 27-pixels neighborhood (0 if out of image)
	 */
	private byte[] getNeighborhood(ImageStack image, int x, int y, int z,
			int w, int h, int d) {
		byte[] neighborhood = new byte[27];

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

}
