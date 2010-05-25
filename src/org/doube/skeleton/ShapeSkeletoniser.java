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

		// label 6-connected white 18-neighbours by particle
		int id = 1;
		int[] labels = new int[27];
		for (int i = 0; i < 18; i++) {
			// find the 18-neighbour in the 27 neighbourhood
			final int p = eighteenNeighbours[i];
			if (neighbours[p] == WHITE) {
				// assign an initial value
				labels[p] = id;
				int minId = id;

				// check 6-neighbours
				int[] neigh6 = conn6Neigh18LUT[p];
				for (int n = 0; n < neigh6.length; n++) {
					final int n6 = neigh6[n];
					if (neighbours[n6] == WHITE) {
						final int nId = labels[n6];
						if (nId != 0 && nId < minId) {
							minId = nId;
						}
					}
				}
				labels[p] = minId;
				if (minId == id)
					id++;
			}
		}
		// find minimum label in each 18-neighbour's 6-neighbourhood
		for (int i = 0; i < 18; i++) {
			final int p = eighteenNeighbours[i];
			if (neighbours[p] == WHITE) {
				int minId = labels[p];
				int[] neigh6 = conn6Neigh18LUT[p];
				for (int n = 0; n < neigh6.length; n++) {
					final int n6 = neigh6[n];
					if (neighbours[n6] == WHITE) {
						final int nId = labels[n6];
						if (nId < minId)
							minId = nId;
					}
				}
				for (int n = 0; n < neigh6.length; n++) {
					final int n6 = neigh6[n];
					final int nId = labels[n6];
					if (nId != minId && nId != 0) {
						for (int j = 0; j < 27; j++) {
							if (labels[j] == nId) {
								labels[j] = minId;
							}
						}
					}
				}
			}
		}
		int b = 0;
		for (int i = 0; i < 6; i++) {
			int a = labels[sixNeighbours[i]];
			if (a != 0) {
				if (b != 0 && a != b)
					return false;
				b = a;
			}
		}
		return true;
	}

	/**
	 * Check the black elements of the 26 neighbourhood and return true if they
	 * form a single 26-connected set
	 * 
	 * @param neighbours
	 * @return
	 */
	private boolean isBlack26ConnectedSet(byte[] neighbours) {
		// join 26 black neighbours into particles
		int id = 1;
		int[] labels = new int[27];
		for (int p = 0; p < 27; p++) {
			if (p == 13)
				continue;
			if (neighbours[p] == BLACK) {
				// assign an initial value
				labels[p] = id;
				int minId = id;

				// check 26-neighbours
				int[] neigh26 = conn26Neigh26LUT[p];
				final int nNeigh = neigh26.length;
				for (int n = 0; n < nNeigh; n++) {
					final int n26 = neigh26[n];
					if (neighbours[n26] == BLACK) {
						final int nId = labels[n26];
						if (nId != 0 && nId < minId) {
							minId = nId;
						}
					}
				}
				labels[p] = minId;
				if (minId == id)
					id++;
			}
		}
		// find minimum label in each 26-neighbour's 26-neighbourhood
		for (int p = 0; p < 27; p++) {
			if (p == 13)
				continue;
			if (neighbours[p] == BLACK) {
				int minId = labels[p];
				int[] neigh26 = conn26Neigh26LUT[p];
				final int nNeigh = neigh26.length;
				for (int n = 0; n < nNeigh; n++) {
					final int n26 = neigh26[n];
					if (neighbours[n26] == BLACK) {
						final int nId = labels[n26];
						if (nId < minId)
							minId = nId;
					}
				}
				for (int n = 0; n < nNeigh; n++) {
					final int n26 = neigh26[n];
					final int nId = labels[n26];
					if (nId != minId && nId != 0) {
						for (int j = 0; j < 27; j++) {
							if (labels[j] == nId) {
								labels[j] = minId;
							}
						}
					}
				}
			}
		}
		int b = 0;
		for (int i = 0; i < 27; i++) {
			if (i == 13)
				continue;
			int a = labels[i];
			if (a != 0) {
				if (b != 0 && a != b)
					return false;
				b = a;
			}
		}
		return true;
	}
	
	//-----------------------------------------------------------------//
	//Look-up tables
	
	/**
	 * List of 18-neighbour points that are 6 connected to the input point,
	 * which must also be an 18-neighbourhood point
	 */
	private static final int[][] conn6Neigh18LUT = { null, { 4, 10 }, null,
			{ 4, 12 }, { 1, 3, 5, 7 }, { 4, 14 }, null, { 4, 16 }, null,
			{ 10, 12 }, { 1, 9, 11, 19 }, { 10, 14 }, { 3, 9, 15, 21 }, null,
			{ 5, 11, 17, 23 }, { 12, 16 }, { 7, 15, 17, 25 }, { 14, 16 }, null,
			{ 10, 22 }, null, { 12, 22 }, { 19, 21, 23, 25 }, { 14, 22 }, null,
			{ 16, 22 }, null };

	/**
	 * List of 26-neighbour points that are 26 connected to the input point,
	 * which must also be an 26-neighbourhood point
	 */
	private static final int[][] conn26Neigh26LUT = { { 1, 3, 4, 9, 10, 12 }, // 0
			{ 0, 2, 3, 4, 5, 9, 10, 11, 12, 14 }, // 1
			{ 1, 4, 5, 10, 11, 14 },// 2
			{ 0, 1, 4, 6, 7, 9, 10, 12, 15, 16 },// 3
			{ 0, 1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 14, 15, 16, 17 },// 4
			{ 1, 2, 4, 7, 8, 10, 11, 14, 16, 17 }, // 5
			{ 3, 4, 7, 12, 15, 16 },// 6
			{ 3, 4, 5, 6, 8, 12, 14, 15, 16, 17 },// 7
			{ 4, 5, 7, 14, 16, 17 },// 8
			{ 0, 1, 3, 4, 10, 12, 18, 19, 21, 22 },// 9
			{ 0, 1, 2, 3, 4, 5, 9, 11, 12, 14, 18, 19, 20, 21, 22, 23 },// 10
			{ 1, 2, 4, 5, 10, 14, 19, 20, 22, 23 },// 11
			{ 0, 1, 3, 4, 6, 7, 9, 10, 15, 16, 18, 19, 21, 22, 24, 25 },// 12
			{},// 13
			{ 1, 2, 4, 5, 7, 8, 10, 11, 16, 17, 19, 20, 22, 23, 25, 26 },// 14
			{ 3, 4, 6, 7, 12, 16, 21, 22, 24, 25 },// 15
			{ 3, 4, 5, 6, 7, 8, 12, 14, 15, 17, 21, 22, 23, 24, 25, 26 },// 16
			{ 4, 5, 7, 8, 14, 16, 22, 23, 25, 26 },// 17
			{ 9, 10, 12, 19, 21, 22 },// 18
			{ 9, 10, 11, 12, 14, 18, 20, 21, 22, 23 },// 19
			{ 10, 11, 14, 19, 22, 23 },// 20
			{ 9, 10, 12, 15, 16, 18, 19, 22, 24, 25 },// 21
			{ 9, 10, 11, 12, 14, 15, 16, 17, 18, 19, 20, 21, 23, 24, 25, 26 }, // 22
			{ 10, 11, 14, 16, 17, 19, 20, 22, 25, 26 },// 23
			{ 12, 15, 16, 21, 22, 25 },// 24
			{ 12, 14, 15, 16, 17, 21, 22, 23, 24, 26 },// 25
			{ 14, 16, 17, 22, 23, 25 } // 26
	};

	/** LUT to find the 6 neighbours in a 27 neighbourhood array */
	private static final int[] sixNeighbours = { 4, 10, 12, 14, 16, 22 };

	/** LUT to find the 18 neighbours in a 27 neighbourhood array */
	private static final int[] eighteenNeighbours = { 1, 3, 4, 5, 7, 9, 10, 11,
			12, 14, 15, 16, 17, 19, 21, 22, 23, 25 };

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
