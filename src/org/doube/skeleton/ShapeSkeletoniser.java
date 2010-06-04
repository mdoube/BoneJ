package org.doube.skeleton;

import java.util.Arrays;

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
		ImagePlus outImp = shapeSkeletonize(imp);
		outImp.show();
	}

	/**
	 * Main skeletonisation method
	 * 
	 * @param imp
	 * @return
	 */
	public ImagePlus shapeSkeletonize(ImagePlus imp) {
		ImageCheck ic = new ImageCheck();
		if (!ic.isBinary(imp))
			throw new IllegalArgumentException();
		final int w = imp.getWidth();
		final int h = imp.getHeight();
		final int d = imp.getImageStackSize();
		short[][] workArray = initialiseWorkArray(imp);
		short threshold = Short.MIN_VALUE + 1;
		long nDeletable = Long.MAX_VALUE;
		short iteration = 1;

		// "scan" until there are no more deletable voxels
		// scan means go through all voxels in volume
		// iterate means traverse all surface points in a topological way and
		// may require multiple scans
		while (nDeletable > 0) {
			long dc = 0;
			for (int z = 0; z < d; z++) {
				for (int y = 0; y < h; y++) {
					final int yw = y * w;
					for (int x = 0; x < w; x++) {
						if (isDeletable(workArray, threshold, x, y, z, w, h, d)) {
							workArray[z][yw + x] = threshold;
							dc++;
						}
					}
				}
			}
			iteration++;
			threshold++;
			nDeletable = dc;
		}
		return imp;
	}

	/**
	 * Create work array, converting original white values to large negative
	 * number and black values to 0
	 * 
	 * @param imp
	 * @return
	 */
	private short[][] initialiseWorkArray(ImagePlus imp) {
		ImageStack stack = imp.getImageStack();
		final int w = stack.getWidth();
		final int h = stack.getHeight();
		final int d = stack.getSize();
		final int wh = w * h;
		short[][] workArray = new short[d][wh];
		for (int z = 0; z < d; z++) {
			byte[] pixels = (byte[]) stack.getPixels(z + 1);
			for (int i = 0; i < wh; i++) {
				if (pixels[i] == WHITE) {
					workArray[z][i] = Short.MIN_VALUE;
				} else {
					workArray[z][i] = 0;
				}
			}
		}
		return workArray;
	}

	/**
	 * Determine if a point (x, y, z) in the work array is marked
	 * 
	 * @param workArray
	 * @param x
	 * @param y
	 * @param z
	 * @param w
	 * @param h
	 * @param d
	 * @return
	 */
	private boolean isMarked(short[][] workArray, int x, int y, int z, int w,
			int h, int d) {
		if (workArray[z][y * w + x] > 0)
			return true;
		else
			return false;
	}

	/**
	 * Determine if a voxel is black (true) or white (false)
	 * 
	 * @param workArray
	 * @param x
	 * @param y
	 * @param z
	 * @param w
	 * @param h
	 * @param d
	 * @return
	 */
	private boolean isBlack(short[][] workArray, int x, int y, int z, int w,
			int h, int d) {
		if (workArray[z][y * w + x] >= 0)
			return true;
		else
			return false;
	}

	/**
	 * Determine if the point at (x, y, z) in the work array is deletable
	 * 
	 * @param workArray
	 * @param x
	 * @param y
	 * @param z
	 * @param w
	 * @param h
	 * @param d
	 * @return
	 */
	private boolean isDeletable(short[][] workArray, long threshold, int x,
			int y, int z, int w, int h, int d) {
		return false;
	}

	/**
	 * Determine if a black point p is s-open, which means at least one of its
	 * 6-neighbours is white prior to iteration
	 * 
	 * @param neighbours
	 * @return
	 */
	private boolean isSOpen(byte[] neighbours) {
		for (int i = 0; i < 6; i++) {
			if (neighbours[sixNeighbours[i]] == WHITE)
				return true;
		}
		return false;
	}

	/**
	 * Determine if a black point p is e-open p is e-open if and only if: p is
	 * not s open
	 * 
	 * @param x
	 *            coordinate of p
	 * @param y
	 *            coordinate of p
	 * @param z
	 *            coordinate of p
	 * @param neighbours
	 * @return
	 */
	private boolean isEOpen(byte[] neighbours, ImageStack stack, int x, int y,
			int z, int w, int h, int d) {
		if (isSOpen(neighbours))
			return false;
		final byte[] f1Points = getF1SPoints(stack, x, y, z, w, h, d);
		for (int i = 0; i < 12; i++) {
			final int n = ePoints[i];
			if (neighbours[n] == WHITE) {
				if (f1Points[edgeSixes[n][0]] == BLACK
						&& f1Points[edgeSixes[n][1]] == BLACK) {
					return true;
				}
			}
		}
		return false;
	}

	/**
	 * Check if a point is v-open
	 * 
	 * @param neighbours
	 * @param stack
	 * @param x
	 * @param y
	 * @param z
	 * @param w
	 * @param h
	 * @param d
	 * @return
	 */
	private boolean isVOpen(byte[] neighbours, ImageStack stack, int x, int y,
			int z, int w, int h, int d) {
		if (isSOpen(neighbours) || isEOpen(neighbours, stack, x, y, z, w, h, d))
			return false;

		byte[][] f1Points = getF1VPoints(stack, x, y, z, w, h, d);
		// check v-points
		checkNeighbours: for (int i = 0; i < 8; i++) {
			final int n = vPoints[i];
			// if white
			if (neighbours[n] == WHITE) {
				// check 3 f1s for blackness
				for (int f = 0; f < 3; f++) {
					if (f1Points[n][f] == WHITE) {
						continue checkNeighbours;
					}
				}
				// return true if all f1s are black
				return true;
			}
		}
		return false;
	}

	/**
	 * Get the f<sub>1</sub> of the s-points around p, which are the points that
	 * are 6-adjacent to the 6-neighbours, but which are not in p's
	 * neighbourhood.
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
	private byte[] getF1SPoints(ImageStack stack, int x, int y, int z, int w,
			int h, int d) {
		byte[] f1Points = new byte[27];
		f1Points[4] = getPixel(stack, x, y, z - 2, w, h, d);
		f1Points[10] = getPixel(stack, x, y - 2, z, w, h, d);
		f1Points[12] = getPixel(stack, x - 2, y, z, w, h, d);
		f1Points[14] = getPixel(stack, x + 2, y, z, w, h, d);
		f1Points[16] = getPixel(stack, x, y + 2, z, w, h, d);
		f1Points[22] = getPixel(stack, x, y, z + 2, w, h, d);
		return f1Points;
	}

	/**
	 * Get the f<sub>1</sub> of the v-points around p, which are the points that
	 * are 6-adjacent to the vertex neighbours, but which are not in p's
	 * neighbourhood.
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
	private byte[][] getF1VPoints(ImageStack stack, int x, int y, int z, int w,
			int h, int d) {
		byte[][] f1Points = new byte[27][3];
		f1Points[0][0] = getPixel(stack, x - 1, y - 1, z - 2, w, h, d);
		f1Points[0][1] = getPixel(stack, x - 1, y - 2, z - 1, w, h, d);
		f1Points[0][2] = getPixel(stack, x - 2, y - 1, z - 1, w, h, d);
		f1Points[2][0] = getPixel(stack, x + 1, y - 1, z - 2, w, h, d);
		f1Points[2][1] = getPixel(stack, x + 1, y - 2, z - 1, w, h, d);
		f1Points[2][2] = getPixel(stack, x + 2, y - 1, z - 1, w, h, d);
		f1Points[6][0] = getPixel(stack, x - 1, y + 1, z - 2, w, h, d);
		f1Points[6][1] = getPixel(stack, x - 1, y + 2, z - 1, w, h, d);
		f1Points[6][2] = getPixel(stack, x - 2, y + 1, z - 1, w, h, d);
		f1Points[8][0] = getPixel(stack, x + 1, y + 1, z - 2, w, h, d);
		f1Points[8][1] = getPixel(stack, x + 1, y + 2, z - 1, w, h, d);
		f1Points[8][2] = getPixel(stack, x + 2, y + 1, z - 1, w, h, d);
		f1Points[18][0] = getPixel(stack, x - 1, y - 1, z + 2, w, h, d);
		f1Points[18][1] = getPixel(stack, x - 1, y - 2, z + 1, w, h, d);
		f1Points[18][2] = getPixel(stack, x - 2, y - 1, z + 1, w, h, d);
		f1Points[20][0] = getPixel(stack, x + 1, y - 1, z + 2, w, h, d);
		f1Points[20][1] = getPixel(stack, x + 1, y - 2, z + 1, w, h, d);
		f1Points[20][2] = getPixel(stack, x + 2, y - 1, z + 1, w, h, d);
		f1Points[24][0] = getPixel(stack, x - 1, y + 1, z + 2, w, h, d);
		f1Points[24][1] = getPixel(stack, x - 1, y + 2, z + 1, w, h, d);
		f1Points[24][2] = getPixel(stack, x - 2, y + 1, z + 1, w, h, d);
		f1Points[26][0] = getPixel(stack, x + 1, y + 1, z + 2, w, h, d);
		f1Points[26][1] = getPixel(stack, x + 1, y + 2, z + 1, w, h, d);
		f1Points[26][2] = getPixel(stack, x + 2, y + 1, z + 1, w, h, d);
		return f1Points;
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
	 * @param neighbours
	 *            from Ignacio's getNeighborhood() method
	 * @return
	 */
	@SuppressWarnings("unused")
	private boolean isSimplePoint(byte[] neighbours) {

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

		// Condition 4. check that white 6 neighbours are 6-connected in the set
		// of white 18-neighbours
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

	/**
	 * Find the values of the voxels in the middle plane, which is a plane
	 * perpendicular to p and 2 opposite s-neighbours, p0 and p1
	 * 
	 * @param neighbours
	 * @param p0
	 *            zeroth s-neighbour, must have smaller neighbour index than p1
	 * @param p1
	 *            first s-neighbour, must have larger neighbour index than p0
	 * @return list of voxel values in the middle plane
	 * @throws IllegalArgumentException
	 *             if p0 and p1 are not a pair of opposite s-neighbours or p0 >=
	 *             p1
	 */
	private byte[] getMiddlePlane(byte[] neighbours, int p0, int p1) {
		byte[] middlePlane = new byte[8];
		int[] plane = new int[8];
		if (p0 == 4 && p1 == 22)
			plane = middlePlanes[0];
		else if (p0 == 10 && p1 == 16)
			plane = middlePlanes[1];
		else if (p0 == 12 && p1 == 14)
			plane = middlePlanes[2];
		else
			throw new IllegalArgumentException();

		for (int i = 0; i < 8; i++)
			middlePlane[i] = neighbours[plane[i]];

		return middlePlane;
	}

	// -----------------------------------------------------------------//
	// Look-up tables

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

	/**
	 * LUT to find the e-points, which are 18-adjacent but not 6-adjacent to p
	 * (i.e. share an edge with p)
	 */
	private static final int[] ePoints = { 1, 3, 5, 7, 9, 11, 15, 17, 19, 21,
			23, 25 };

	/**
	 * LUT to find the v-points, which are 26-adjacent but not 18-adjacent to p
	 * (i.e. share only a vertex with p)
	 */
	private static final int[] vPoints = { 0, 2, 6, 8, 18, 20, 24, 26 };

	/**
	 * LUT to find the pair of 6-neighbours of p that are six neighbours of each
	 * edge point
	 */
	private static final int[][] edgeSixes = { null,// 0
			{ 4, 10 }, // 1
			null,// 2
			{ 4, 12 }, // 3
			null,// 4
			{ 4, 14 },// 5
			null,// 6
			{ 4, 16 },// 7
			null,// 8
			{ 10, 12 },// 9
			null,// 10
			{ 10, 14 },// 11
			null,// 12
			null,// 13
			null,// 14
			{ 12, 16 },// 15
			null,// 16
			{ 14, 16 },// 17
			null,// 18
			{ 10, 22 },// 19
			null,// 20
			{ 12, 22 },// 21
			null,// 22
			{ 14, 22 },// 23
			null,// 24
			{ 16, 22 },// 25
			null // 26
	};

	/**
	 * LUT to find the 6 neighbours in a 27 neighbourhood array. s-points are
	 * 6-adjacent to p (i.e. share a face with p)
	 */
	private static final int[] sixNeighbours = { 4, 10, 12, 14, 16, 22 };

	/**
	 * LUT to find the 18 neighbours in a 27 neighbourhood array, these points
	 * share an edge and/or a face with p.
	 */
	private static final int[] eighteenNeighbours = { 1, 3, 4, 5, 7, 9, 10, 11,
			12, 14, 15, 16, 17, 19, 21, 22, 23, 25 };

	/**
	 * LUT of the 3 middle planes. 0 is the xy plane, 1 is the xz plane and 2 is
	 * the yz plane.
	 */
	private static final int[][] middlePlanes = {
			{ 9, 10, 11, 12, 14, 15, 16, 17 }, { 3, 4, 5, 12, 14, 21, 22, 23 },
			{ 1, 4, 7, 10, 16, 19, 22, 25 } };

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
