package org.doube.bonej;

/**
 * BranchLabeler plugin for ImageJ
 * Copyright 2010 Michael Doube 
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

import java.util.ArrayList;
import java.util.Iterator;

import org.doube.skeleton.AnalyzeSkeleton;
import org.doube.skeleton.Edge;
import org.doube.skeleton.Graph;
import org.doube.skeleton.Point;
import org.doube.skeleton.Skeletonize3D;
import org.doube.util.ImageCheck;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;

/**
 * Perform branch labelling (individual trabecula segmentation). Start with a
 * medial axis skeleton, label junctions and branches, then dilate each branch
 * until the original structure is filled.
 * 
 * @author Michael Doube
 * 
 */
public class BranchLabeler implements PlugIn {
	static int width, height, depth;

	public void run(String arg) {
		if (!ImageCheck.checkEnvironment())
			return;
		ImagePlus imp = IJ.getImage();
		ImageCheck ic = new ImageCheck();
		if (!ic.isBinary(imp)) {
			IJ.error("Branch Labeler requires a binary image.");
			return;
		}
		width = imp.getWidth();
		height = imp.getHeight();
		depth = imp.getImageStackSize();
		GenericDialog gd = new GenericDialog("Options");
		gd.addCheckbox("Prune_ends", false);
		gd.showDialog();
		if (gd.wasCanceled())
			return;
		boolean pruneEnds = gd.getNextBoolean();

		ImageStack stack = imp.getImageStack();
		// skeletonise original
		Skeletonize3D skel = new Skeletonize3D();
		ImagePlus skeleton = skel.getSkeleton(imp);

		// Process the skeleton and retrieve the Graphs for each tree
		AnalyzeSkeleton as = new AnalyzeSkeleton();
		as.run(AnalyzeSkeleton.NONE, pruneEnds, skeleton, imp, true, false);
		Graph[] graphs = as.getGraphs();
		skeleton.close();
		skeleton.flush();

		// The image stack where we will keep the labelled branches
		ImageStack its = new ImageStack(width, height, depth);
		for (int s = 1; s <= depth; s++) {
			float[] pixels = new float[width * height];
			for (int i = 0; i < pixels.length; i++)
				pixels[i] = (int) 0;
			its.setPixels(pixels, s);
			its.setSliceLabel("" + s, s);
		}

		// iteratively grow neighbourhoods from labelled branches
		// until all the space in imp is filled
		int labels = 2;
		final int nTrees = graphs.length;
		for (int t = 0; t < nTrees; t++) {
			Graph graph = graphs[t];

			// Initial labelling of branches
			final int nBranches = graph.getEdges().size();
			int nSeeds = 0;
			ArrayList<Edge> branches = graph.getEdges();
			for (int b = 0; b < nBranches; b++) {
				Edge branch = branches.get(b);
				ArrayList<Point> points = branch.getSlabs();
				nSeeds += points.size();
				final int nPoints = points.size();
				for (int i = 0; i < nPoints; i++) {
					Point p = points.get(i);
					setPixel(its, p.x, p.y, p.z, labels + b);
				}
			}

			while (nSeeds > 0) {
				nSeeds = countSlabs(branches);
				for (int b = 0; b < nBranches; b++) {
					Edge branch = branches.get(b);
					ArrayList<Point> seeds = branch.getSlabs();
					Iterator<Point> sit = seeds.listIterator();
					ArrayList<Point> nextSeeds = new ArrayList<Point>();
					while (sit.hasNext()) {
						Point seed = sit.next();
						// check seed's neighbours in original image
						// and labelled image
						for (int z = seed.z - 1; z <= seed.z + 1; z++) {
							for (int y = seed.y - 1; y <= seed.y + 1; y++) {
								for (int x = seed.x - 1; x <= seed.x + 1; x++) {
									final float label = getFloatPixel(its, x,
											y, z);
									final byte template = getBytePixel(stack,
											x, y, z);
									if (label == 0 && template != 0) {
										setPixel(its, x, y, z, labels + b);
										nextSeeds.add(new Point(x, y, z));
									}
								}
							}
						}
						// once checked, remove the seed
						sit.remove();
					}
					// replace the old seeds with new seeds for the next
					// iteration
					for (Point s : nextSeeds) {
						seeds.add(s);
					}
				}
			}
			labels += nBranches;
		}

		ImagePlus itsImp = new ImagePlus();
		itsImp.setTitle("ITS");
		itsImp.setStack(its);
		itsImp.show();
		itsImp.setDisplayRange(0, labels);
		IJ.run("Fire");
		// return an imp containing the labelled branches
	}

	private static int countSlabs(ArrayList<Edge> branches) {
		int nSlabs = 0;
		final int nBranches = branches.size();
		for (int b = 0; b < nBranches; b++) {
			nSlabs += branches.get(b).getSlabs().size();
		}
		return nSlabs;
	}

	/**
	 * Set pixel in 3D 32-bit image.
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
	private void setPixel(ImageStack image, int x, int y, int z, float value) {
		if (x >= 0 && x < width && y >= 0 && y < height && z >= 0 && z < depth)
			((float[]) image.getPixels(z + 1))[x + y * width] = value;
	} // end setPixel

	/**
	 * Get pixel in 3D byte (8-bit) image (0 border conditions)
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
	private static byte getBytePixel(final ImageStack image, final int x,
			final int y, final int z) {
		if (x >= 0 && x < width && y >= 0 && y < height && z >= 0 && z < depth)
			return ((byte[]) image.getPixels(z + 1))[x + y * width];
		else
			return 0;
	}

	/**
	 * Get pixel in 3D float (32-bit) image (0 border conditions)
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
	private static float getFloatPixel(final ImageStack image, final int x,
			final int y, final int z) {
		final int width = image.getWidth();
		final int height = image.getHeight();
		final int depth = image.getSize();
		if (x >= 0 && x < width && y >= 0 && y < height && z >= 0 && z < depth)
			return ((float[]) image.getPixels(z + 1))[x + y * width];
		else
			return 0;
	}
}
