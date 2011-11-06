package org.doube.util;

import java.awt.Rectangle;
import java.util.ArrayList;

import ij.plugin.PlugIn;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;
import ij.gui.Roi;
import ij.plugin.filter.ThresholdToSelection;
import ij.plugin.frame.RoiManager;

public class RoiInterpolator implements PlugIn {
	int[][] idt;
//	int w, h;
	
	ImagePlus image;
	int w, h, d;
	int current, total;


	public void run(String arg) {
		RoiManager roiman = RoiManager.getInstance();
		if (roiman == null){
			IJ.error("Please populate the ROI Manager with ROIs");
			return;
		}
		Roi[] rois = roiman.getRoisAsArray();
		int xmax = 0;
		int xmin = Integer.MAX_VALUE;
		int ymax = 0;
		int ymin = Integer.MAX_VALUE;
		int zmax = 1;
		int zmin = Integer.MAX_VALUE;
		ArrayList<Integer> templateSlices = new ArrayList<Integer>();
		for (Roi roi : rois){
			final int slice = roiman.getSliceNumber(roi.getName());
			templateSlices.add(new Integer(slice));
			if (slice == 0) //ignore non-slice associated ROIs
				continue;
			zmin = Math.min(slice, zmin);
			zmax = Math.max(slice, zmax);
			Rectangle bounds = roi.getBounds();
			xmin = Math.min(xmin, bounds.x);
			ymin = Math.min(ymin, bounds.y);
			xmax = Math.max(xmax, bounds.x + bounds.width);
			ymax = Math.max(ymax, bounds.y + bounds.height);
		}
		//create the binary stack
		final int nSlices = zmax - zmin + 1;
		ImageStack stack = new ImageStack(xmax, ymax);
		for (int s = 0; s < nSlices; s++){
			ByteProcessor bp = new ByteProcessor(xmax, ymax);
			bp.setColor(255);
			for (Roi roi : rois){
				final int slice = roiman.getSliceNumber(roi.getName());
				if (slice == zmin + s){
					bp.setRoi(roi);
					if (roi.getType() == Roi.RECTANGLE)
						bp.fill();
					else 
						bp.fill(roi);
				}
			}
			stack.addSlice(""+s, bp);
		}
		//do the binary interpolation
		run(stack);
		ImagePlus binary = new ImagePlus("interpolated", stack);

		//get the ROIs
		ThresholdToSelection ts = new ThresholdToSelection();
		ts.setup("", binary);
		for(int s = 0; s < nSlices; s++) {
			if (templateSlices.contains(new Integer(s + zmin)))
				continue;
			ImageProcessor bp = stack.getProcessor(s+1);
			int threshold = 255;
			bp.setThreshold(threshold, threshold, ImageProcessor.NO_LUT_UPDATE);
			Roi roi = ts.convert(bp);
			roi.setPosition(s + zmin);
			roiman.addRoi(roi);
		}
		IJ.showStatus("ROIs interpolated");
	}

//------------------ From Fiji --------------------------//
//http://fiji.sc/cgi-bin/gitweb.cgi?p=fiji.git;a=blob_plain;f=src-plugins/VIB-lib/vib/BinaryInterpolator.java;h=f6a610659ad624d13f94639bc5c0149712071f9f;hb=refs/heads/master

/*
 * This plugin takes a binary stack as input, where some slices are
 * labeled (i.e. contain white regions), and some are not. The unlabaled
 * regions are interpolated by weighting the signed integer distance
 * transformed labeled slices.
 */

	public void run(ImagePlus image, Roi[] rois) {
		w = image.getWidth();
		h = image.getHeight();
		ImageStack stack = new ImageStack(w, h);
		int firstIndex = -1, lastIndex = -1;
		for(int i = 1; i < rois.length; i++) {
			if(rois[i] != null) {
				firstIndex = (firstIndex == -1) ? i : firstIndex;
				lastIndex = i;
			}
		}
		if (firstIndex == -1) {
			IJ.error("There must be at least one selection in order to interpolate.");
			return;
		}

		for(int i = firstIndex; i <= lastIndex; i++) {
			ByteProcessor bp = new ByteProcessor(w, h);
			if(rois[i] != null) {
				bp.copyBits(rois[i].getMask(), 
							rois[i].getBounds().x, 
							rois[i].getBounds().y, 
							ij.process.Blitter.ADD);
			}
			stack.addSlice("", bp);
		}
		run(stack);
		ImagePlus roiImage = new ImagePlus("bla", stack);
		
		ThresholdToSelection ts = new ThresholdToSelection();
		ts.setup("", roiImage);
		for(int i = firstIndex; i <= lastIndex; i++) {
			ImageProcessor bp = stack.getProcessor(1);
			stack.deleteSlice(1);
			int threshold = 255;
			bp.setThreshold(threshold, threshold, ImageProcessor.NO_LUT_UPDATE);
			ts.run(bp);
			rois[i] = roiImage.getRoi();
		}
	}

	public void run(ImageStack stack) {	
		int sliceCount = stack.getSize();
		if (sliceCount < 3) {
			IJ.error("Too few slices to interpolate!");
			return;
		}

		IJ.showStatus("getting signed integer distance transform");
		w = stack.getWidth();
		h = stack.getHeight();
		idt = new int[sliceCount][];
		int first = sliceCount, last = -1;

		for (int z = 0; z < sliceCount; z++) {
			idt[z] = getIDT(stack.getProcessor(z + 1).getPixels());
			if (idt[z] != null) {
				if (z < first)
					first = z;
				last = z;
			}
		 }

		if (first == last || last < 0) {
			IJ.error("Not enough to interpolate");
			return;
		}

		IJ.showStatus("calculating weights");
		int current = 0, next = first;
		for (int z = first; z < last; z++) {
			if (z == next) {
				current = z;
				for (next = z + 1; idt[next] == null; next++);
				continue;
			}

			byte[] p =
				(byte[])stack.getProcessor(z + 1).getPixels();
			for (int i = 0; i < w * h; i++)
				if (0 <= idt[current][i] * (next - z)
						+ idt[next][i] * (z - current))
					p[i] = (byte)255;
			IJ.showProgress(z - first + 1, last - z);
		}
	}

	/*
	 * The following calculates the signed integer distance transform.
	 * Distance transform means that each pixel is assigned the distance
	 * to the boundary.
	 * IDT means that the distance is not the Euclidean, but the minimal
	 * sum of neighbour distances with 3 for horizontal and neighbours,
	 * and 4 for diagonal neighbours (in 3d, the 3d diagonal neighbour
	 * would be 5).
	 * Signed means that the outside pixels have a negative sign.
	 */
	class IDT {
		int[] result;

		IDT() {
			result = new int[w * h];
			int infinity = (w + h) * 9;

			for (int i = 0; i < result.length; i++)
				result[i] = infinity;
		}

		int init(byte[] p) {
			int count = 0;

			for (int j = 0; j < h; j++)
				for (int i = 0; i < w; i++) {
					int idx = i + w * j;
					if (isBoundary(p, i, j)) {
						result[idx] = 0;
						count++;
					} else if (isJustOutside(p, i, j))
						result[idx] = -1;
				}
			return count;
		}

		final void idt(int x, int y, int dx, int dy) {
			if (x + dx < 0 || y + dy < 0 ||
					x + dx >= w || y + dy >= h)
				return;
			int value = result[x + dx + w * (y + dy)];
			int distance = (dx == 0 || dy == 0 ? 3 : 4);
			value += distance * (value < 0 ? -1 : 1);
			if (Math.abs(result[x + w * y]) > Math.abs(value))
				result[x + w * y] = value;
		}

		void propagate() {
			for (int j = 0; j < h; j++)
				for (int i = 0; i < w; i++) {
					idt(i, j, -1, 0);
					idt(i, j, -1, -1);
					idt(i, j, 0, -1);
				}

			for (int j = h - 1; j >= 0; j--)
				for (int i = w - 1; i >= 0; i--) {
					idt(i, j, +1, 0);
					idt(i, j, +1, +1);
					idt(i, j, 0, +1);
				}

			for (int i = w - 1; i >= 0; i--)
				for (int j = h - 1; j >= 0; j--) {
					idt(i, j, +1, 0);
					idt(i, j, +1, +1);
					idt(i, j, 0, +1);
				}

			for (int i = 0; i < w; i++)
				for (int j = 0; j < h; j++) {
					idt(i, j, -1, 0);
					idt(i, j, -1, -1);
					idt(i, j, 0, -1);
				}
		}
	}

	int[] getIDT(Object pixels) {
		IDT idt = new IDT();
		if (idt.init((byte[])pixels) == 0)
			return null;
		idt.propagate();
		return idt.result;
	}

	final boolean isBoundary(byte[] pixels, int x, int y) {
		if (pixels[x + w * y] == 0)
			return false;
		if (x <= 0 || pixels[x - 1 + w * y] == 0)
			return true;
		if (x >= w - 1 || pixels[x + 1 + w * y] == 0)
			return true;
		if (y <= 0 || pixels[x + w * (y - 1)] == 0)
			return true;
		if (y >= h - 1 || pixels[x + w * (y + 1)] == 0)
			return true;
		if (x <= 0 || y <= 0 || pixels[x - 1 + w * (y - 1)] == 0)
			return true;
		if (x <= 0 || y >= h - 1 || pixels[x - 1 + w * (y + 1)] == 0)
			return true;
		if (x >= w - 1 || y <= 0 || pixels[x + 1 + w * (y - 1)] == 0)
			return true;
		if (x >= w - 1 || y >= h - 1 ||
				pixels[x + 1 + w * (y + 1)] == 0)
			return true;
		return false;
	}

	final boolean isJustOutside(byte[] pixels, int x, int y) {
		if (pixels[x + w * y] != 0)
			return false;
		if (x > 0 && pixels[x - 1 + w * y] != 0)
			return true;
		if (x < w - 1 && pixels[x + 1 + w * y] != 0)
			return true;
		if (y > 0 && pixels[x + w * (y - 1)] != 0)
			return true;
		if (y < h - 1 && pixels[x + w * (y + 1)] != 0)
			return true;
		if (x > 0 && y > 0 && pixels[x - 1 + w * (y - 1)] != 0)
			return true;
		if (x > 0 && y < h - 1 && pixels[x - 1 + w * (y + 1)] != 0)
			return true;
		if (x < w - 1 && y > 0 && pixels[x + 1 + w * (y - 1)] != 0)
			return true;
		if (x < w - 1 && y < h - 1 &&
				pixels[x + 1 + w * (y + 1)] != 0)
			return true;
		return false;
	}
	
	//--------------- From Fiji --------------------------
	//http://fiji.sc/cgi-bin/gitweb.cgi?p=fiji.git;a=blob_plain;f=src-plugins/Fiji_Plugins/fiji/process3d/EDT.java;h=eb433b2912a3a7955e0c56b91b009d9ff6fabb5e;hb=refs/heads/master

	/*
	 * The idea of the Euclidean Distance Transform is to get the
	 * distance of every outside pixel to the nearest outside pixel.
	 *
	 * We use the algorithm proposed in

		@TECHREPORT{Felzenszwalb04distancetransforms,
		    author = {Pedro F. Felzenszwalb and Daniel P. Huttenlocher},
		    title = {Distance transforms of sampled functions},
		    institution = {Cornell Computing and Information Science},
		    year = {2004}
		}

	 * Felzenszwalb & Huttenlocher's idea is to extend the concept to
	 * a broader one, namely to minimize not only the distance to an
	 * outside pixel, but to minimize the distance plus a value that
	 * depends on the outside pixel.
	 *
	 * In mathematical terms: we determine the minimum of the term
	 *
	 *	g(x) = min(d^2(x, y) + f(y) for all y)
	 *
	 * where y runs through all pixels and d^2 is the square of the
	 * Euclidean distance. For the Euclidean distance transform, f(y)
	 * is 0 for all outside pixels, and infinity for all inside
	 * pixels, and the result is the square root of g(x).
	 *
	 * The trick is to calculate g in one dimension, store the result,
	 * and use it as f(y) in the next dimension. Continue until you
	 * covered all dimensions.
	 *
	 * In order to find the minimum in one dimension (i.e. row by
	 * row), the following fact is exploited: for two different
	 * y1 < y2, (x - y1)^2 + f(y1) < (x - y2)^2 + f(y2) for x < s,
	 * where s is the intersection point of the two parabolae (there
	 * is the corner case where one parabola is always lower than
	 * the other one, in that case there is no intersection).
	 *
	 * Using this fact, for each row of n elements, a maximum number
	 * of n parabolae are constructed, adding them one by one for each
	 * y, adjusting the range of x for which this y yields the minimum,
	 * possibly overriding a number of previously added parabolae.
	 *
	 * At most n parabolae can be added, so the complexity is still
	 * linear.
	 *
	 * After this step, the list of parabolae is iterated to calculate
	 * the values for g(x).
	 */


		public ImagePlus compute(ImageStack stack) {
			w = stack.getWidth();
			h = stack.getHeight();
			d = stack.getSize();
			ImageStack result = new ImageStack(w, h, d);
			for (int i = 1; i <= d; i++)
				result.setPixels(new float[w * h], i);

			current = 0;
			total = w * h * d * 3;

			new Z(stack, result).compute();
			new Y(result).compute();
			new X(result).compute();

			return new ImagePlus("EDT", result);
		}

		abstract class EDTBase {
			int width;
			/*
			 * parabola k is defined by y[k] (v in the paper)
			 * and f[k] (f(v[k]) in the paper): (y, f) is the
			 * coordinate of the minimum of the parabola.
			 * z[k] determines the left bound of the interval
			 * in which the k-th parabola determines the lower
			 * envelope.
			 */

			int k;
			float[] f, z;
			int[] y;

			EDTBase(int rowWidth) {
				width = rowWidth;
				f = new float[width + 1];
				z = new float[width + 1];
				y = new int[width + 1];
			}

			final void computeRow() {
				// calculate the parabolae ("lower envelope")
				f[0] = Float.MAX_VALUE;
				y[0] = -1;
				z[0] = Float.MAX_VALUE;
				k = 0;
				float fx, s;
				for (int x = 0; x < width; x++) {
					fx = get(x);
					for (;;) {
						// calculate the intersection
						s = ((fx + x * x) - (f[k] + y[k] * y[k])) / 2 / (x - y[k]);
						if (s > z[k])
							break;
						if (--k < 0)
							break;
					}
					k++;
					y[k] = x;
					f[k] = fx;
					z[k] = s;
				}
				z[++k] = Float.MAX_VALUE;
				// calculate g(x)
				int i = 0;
				for (int x = 0; x < width; x++) {
					while (z[i + 1] < x)
						i++;
					set(x, (x - y[i]) * (x - y[i]) + f[i]);
				}
			}

			abstract float get(int column);

			abstract void set(int column, float value);

			final void compute() {
				while (nextRow()) {
					computeRow();
					if (total > 0) {
						current += width;
						IJ.showProgress(current, total);
					}
				}
			}

			abstract boolean nextRow();
		}

		class Z extends EDTBase {
			byte[][] inSlice;
			float[][] outSlice;
			int offset;

			Z(ImageStack in, ImageStack out) {
				super(d);
				inSlice = new byte[d][];
				outSlice = new float[d][];
				for (int i = 0; i < d; i++) {
					inSlice[i] = (byte[])in.getPixels(i + 1);
					outSlice[i] = (float[])out.getPixels(i + 1);
				}
				offset = -1;
			}

			final float get(int x) {
				return inSlice[x][offset] == 0 ? 0 : Float.MAX_VALUE;
			}

			final void set(int x, float value) {
				outSlice[x][offset] = value;
			}

			final boolean nextRow() {
				return ++offset < w * h;
			}
		}

		abstract class OneDimension extends EDTBase {
			ImageStack stack;
			float[] slice;
			int offset, lastOffset, rowStride, columnStride, sliceIndex;

			OneDimension(ImageStack out, boolean iterateX) {
				super(iterateX ? w : h);
				stack = out;
				columnStride = iterateX ? 1 : w;
				rowStride = iterateX ? w : 1;
				offset = w * h;
				lastOffset = rowStride * (iterateX ? h : w);
				sliceIndex = -1;
			}

			final float get(int x) {
				return slice[x * columnStride + offset];
			}

			final boolean nextRow() {
				offset += rowStride;
				if (offset >= lastOffset) {
					if (++sliceIndex >= d)
						return false;
					offset = 0;
					slice = (float[])stack.getPixels(sliceIndex + 1);
				}
				return true;
			}
		}

		class Y extends OneDimension {
			Y(ImageStack out) {
				super(out, false);
			}

			final void set(int x, float value) {
				slice[x * columnStride + offset] = value;
			}
		}

		class X extends OneDimension {
			X(ImageStack out) {
				super(out, true);
			}

			final void set(int x, float value) {
				slice[x * columnStride + offset] = (float)Math.sqrt(value);
			}
		}
}
