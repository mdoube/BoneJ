package org.doube.util;

import java.util.HashSet;
import java.util.Iterator;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.plugin.PlugIn;
import ij.process.ByteProcessor;


public class RayTracer implements PlugIn {

	@Override
	public void run(String arg) {
		// TODO Auto-generated method stub
		findCollisionPoints(IJ.getImage(), new int[]{0,0,0});
	}
	
	
	private int[][] findCollisionPoints(ImagePlus imp, int[] startPoint){
		
		final int w = imp.getWidth();
		final int h = imp.getHeight();
		final int d = imp.getStackSize();
		
		final int sX = startPoint[0];
		final int sY = startPoint[1];
		final int sZ = startPoint[2];
		
		ImageStack stack = imp.getImageStack();
		
		ByteProcessor[] pixels = new ByteProcessor[d]; 
		for (int i = 0; i < d; i++) {
			pixels[i] = (ByteProcessor) stack.getProcessor(i + 1);
		}
		
		
		int i = 0;
		
		//check whether the starting point is a collision point
		//might be problematic if starting from a surface point
		if (isBackground(pixels, sX, sY, sZ)) {
			return new int[][] {startPoint};
		}
		
		//stores the list of live sampling ray vectors and their current locations
		//each with the three centre coordinates followed by the three vector elements
		//and three elements giving the face identity
		HashSet<double[]> integerVectors = new HashSet<double[]>();
		
		//stores the list of collision points
		HashSet<int[]> collisionPoints = new HashSet<int[]>();

		//initial vectors
		for (int z = -1; z <= 1; z++) {
			for (int y = -1; y <= 1; y++) {
				for (int x = -1; x <= 1; x++) {
					if (x == 0 && y == 0 && z == 0)
						continue;
					final int xc = sX + x;
					final int yc = sY + y;
					final int zc = sZ + z;
					
					//last three coordinates are the parent face/edge/corner for later incrementing purposes
					//can determine which face (i.e. +- x, y, or z) or whether a face
					final double[] vector = new double[] {xc, yc, zc, x, y, z, x, y, z};
					integerVectors.add(vector);
				}
			}
		}
		
		i++;
		
		//loop until all the rays are killed by collisions
		while (integerVectors.size() > 0) {
			if (isPowerOfTwo(i))
	       integerVectors = spawn(integerVectors, startPoint);
		  
			checkRaysForCollisions(pixels, integerVectors, collisionPoints, w, h, d);
	    
			i++;
			
			integerVectors = incrementVectors(integerVectors);
		}
		
	  return collisionPoints.toArray(new int[collisionPoints.size()][3]);
	}
	
/**
 * Set the current position to the last position plus integer vector
 * 
 * @param integerVectors
 * @return updated list of vectors
 */
	private HashSet<double[]> incrementVectors(HashSet<double[]> integerVectors) {
		
		HashSet<double[]> nextPosition = new HashSet<double[]>();
		
		Iterator<double[]> iterator = integerVectors.iterator();
		while (iterator.hasNext()) {
			double[] vector = iterator.next();
			
			//set current position to integer vector plus last position
			vector[0] += vector[3];
			vector[1] += vector[4];
			vector[2] += vector[5];
			
			nextPosition.add(vector);
		}

		return nextPosition;
	}
	
	/**
	 * Bisect between existing rays within the current i plane
	 * 
	 * Do not spawn from corners or edges, only from faces
	 * 
	 * @param parentVectors
	 * @return 
	 */
	private HashSet<double[]> spawn(HashSet<double[]> parentVectors, int[] startPoint) {
		
		HashSet<double[]> childVectors = new HashSet<double[]>();
		
		Iterator<double[]> iterator = parentVectors.iterator();
		while (iterator.hasNext()) {
			double[] vector = iterator.next();
			
			//trim any imprecision to set ray to pixel grid
			vector[0] = Math.round(vector[0]);
			vector[1] = Math.round(vector[1]);
			vector[2] = Math.round(vector[2]);
			childVectors.add(vector);
			
			//create 4 child vectors as clones of the parent
			double[] child0 = vector.clone();
			double[] child1 = vector.clone();
			double[] child2 = vector.clone();
			double[] child3 = vector.clone();
			
			//get the directions to +- 0.5
			//if statement includes only faces
			//edges and corners are excluded
			//+- 0.5 in x & y
			if (vector[6] == 0 && vector[7] == 0) {
				child0[0] += 0.5;
				child1[0] -= 0.5;
				child2[1] += 0.5;
				child3[1] -= 0.5;
			}
			//+- 0.5 in x & z
			else if (vector[6] == 0 && vector[8] == 0) {
				child0[0] += 0.5;
				child1[0] -= 0.5;
				child2[2] += 0.5;
				child3[2] -= 0.5;
			}
			//+- 0.5 in y and z
			else if (vector[7] == 0 && vector[8] == 0) {
				child0[1] += 0.5;
				child1[1] -= 0.5;
				child2[2] += 0.5;
				child3[2] -= 0.5;
			}
			//handle edges
			else if (isEdge(vector)) {
				if (vector[6] == 0) {
					//spawn along the edge
					child0[0] += 0.5;
					child1[0] -= 0.5;
					//spawn across faces
					child2[1] -= vector[7] * 0.5;
					child3[2] -= vector[8] * 0.5;
				}
				else if (vector[7] == 0) {
				  //spawn along the edge
					child0[1] += 0.5;
					child1[1] -= 0.5;
					//spawn across faces
					child2[0] -= vector[6] * 0.5;
					child3[2] -= vector[8] * 0.5;
				}
				else if (vector[8] == 0) {
				  //spawn along the edge
					child0[2] += 0.5;
					child1[2] -= 0.5;
					//spawn across faces
					child2[0] -= vector[6] * 0.5;
					child3[1] -= vector[7] * 0.5;
				}
			}

			calculateIntegerVector(child0, startPoint);
			calculateIntegerVector(child1, startPoint);
			calculateIntegerVector(child2, startPoint);
			calculateIntegerVector(child3, startPoint);
						
			//add them to the childVectors HashSet
			childVectors.add(child0);
			childVectors.add(child1);
			childVectors.add(child2);
			childVectors.add(child3);
			
		}
		
		return childVectors;
	}

	/**
	 * Determine whether a vector is an edge
	 * Edges have one 0 and two +-1
	 * 
	 * Check that there is exactly 1 zero
	 *  
	 * @param vector
	 * @return
	 */
	private boolean isEdge(double[] vector) {
		int sum = 0;
		for (int i = 6; i < 9; i++)
			if (vector[i] == 0)
				sum++;
		
		return sum == 1;
	}


	private void calculateIntegerVector(double[] vector, int[] startPoint) {

		double l = 0;
		//in x & y plane, normal is z
		if (vector[6] == 0 && vector[7] == 0) {
			l = vector[2];
		}
		//in x & z plane, normal is y
		else if (vector[6] == 0 && vector[8] == 0) {
			l = vector[1];
		}
		//in y and z plane, normal is x
		else if (vector[7] == 0 && vector[8] == 0) {
			l = vector[0];
		}
		//x component
		vector[3] = (vector[0] - startPoint[0])/l;
		
		//y component
		vector[4] = (vector[1] - startPoint[1])/l;
		
		//z component
		vector[5] = (vector[2] - startPoint[2])/l;
	}


	/**
	 * Check each integer vector.
	 * 
	 * If it collides, remove it from integerVectors and add the 
	 * collision pixel int to collisionPoints.
	 * 
	 * If it is out of bounds remove it from integerVectors
	 * 
	 * @param integerVectors
	 * @param collisionPoints
	 */
	private void checkRaysForCollisions(ByteProcessor[] pixels, HashSet<double[]> integerVectors,
		HashSet<int[]> collisionPoints, final int w, final int h, final int d)
	{
		Iterator<double[]> iterator = integerVectors.iterator();
		while (iterator.hasNext()) {
			double[] vector = iterator.next();
			//floor double values to snap to pixel grid
			final int x = (int) Math.floor(vector[0]);
			final int y = (int) Math.floor(vector[1]);
			final int z = (int) Math.floor(vector[2]);

			if (isOutOfBounds(x, y, z, w, h, d)) {
				iterator.remove();
				continue;
			}
			if (isBackground(pixels, x, y, z)) {
				iterator.remove();
				collisionPoints.add(new int[] {x, y, z});
			}
		}
	}


	private boolean isOutOfBounds(final int x, final int y, final int z, 
		final int w, final int h, final int d) {
		
		return x < 0 || x >= w || y < 0 || y >= h || z < 0 || z >= d;
	}


	/**
	 * Check whether the given point is background
	 * @param pixels
	 * @param startPoint
	 * @return true if the point is background
	 */
	private boolean isBackground(ByteProcessor[] pixels,
		final int x, final int y, final int z) {
		
		return pixels[z].get(x, y) == 0;
	}


	private static boolean isPowerOfTwo(final int number) {
    return number > 0 && ((number & (number - 1)) == 0);
  }
}
