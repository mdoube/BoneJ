package org.doube.util;

import java.util.ArrayList;
import java.util.Arrays;
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

		ImagePlus imp = IJ.getImage();
		
		final int w = imp.getWidth();
		final int h = imp.getHeight();
		final int d = imp.getStackSize();
		
		final int sX = (int)Math.round(w / 2.0);
		final int sY = (int)Math.round(h / 2.0);
		final int sZ = (int)Math.round(d / 2.0);
		
		final long startTime = System.nanoTime();
		int[][] collisionPoints = findCollisionPoints(imp, new int[]{sX,sY,sZ});
		final long endTime = System.nanoTime();
		
		IJ.log("Found "+collisionPoints.length+" collision points");
		if (IJ.debugMode) {
		 for (int[] point : collisionPoints) {
		  	IJ.log("Collision point found at ("+point[0]+", "+point[1]+", "+point[2]+")");
		 }
		}
		IJ.log("Finished RayTracer in "+(endTime - startTime)/1000000.0+" ms");
		
		ByteProcessor[] bytes = new ByteProcessor[d];
		ImageStack stack = new ImageStack(w, h, d);
		
		for (int i = 0; i < d; i++)
			bytes[i] = new ByteProcessor(w, h);
		
		for (int[] point : collisionPoints)
			bytes[point[2]].set(point[0], point[1], 255);
		
		for (int i = 0; i < d; i++)
			stack.setProcessor(bytes[i], i+1);

		ImagePlus impOut = new ImagePlus(imp.getTitle()+"_Collision_Points", stack);
		impOut.show();
		
	}
	
	/**
	 * Given a start point, find the complete set of unique points
	 * resulting from collisions between rays from the start point and
	 * the image background
	 *  
	 * @param imp
	 * @param startPoint in int{x, y, z} format
	 * @return list of collision points 
	 */
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
		HashSet<ArrayList<Double>> integerVectors = new HashSet<ArrayList<Double>>();
		
		//stores the list of collision points
		HashSet<ArrayList<Integer>> collisionPoints = new HashSet<ArrayList<Integer>>();

		//initial vectors
		for (int z = -1; z <= 1; z++) {
			for (int y = -1; y <= 1; y++) {
				for (int x = -1; x <= 1; x++) {
					if (x == 0 && y == 0 && z == 0)
						continue;

					final int xc = sX + x;
					final int yc = sY + y;
					final int zc = sZ + z;
					
					if (IJ.debugMode)
						IJ.log("Set initial ray at ("+xc+", "+yc+", "+zc+", "+x+", "+y+", "+z+", "+x+", "+y+", "+z+")");
					
					//last three coordinates are the parent face/edge/corner for later incrementing purposes
					//can determine which face (i.e. +- x, y, or z) or whether a face
					ArrayList<Double> vector = new ArrayList<Double>(Arrays.asList(
						(double) xc, (double)yc, (double)zc, (double)x, (double)y,
						(double)z, (double)x, (double)y, (double)z));
					
					integerVectors.add(vector);
				}
			}
		}
		
		i++;
		
		//loop until all the rays are killed by collisions
		while (integerVectors.size() > 0) {
			if (IJ.debugMode)
				IJ.log("entering iteration number "+i);
			
			if (isPowerOfTwo(i))
	       integerVectors = spawn(integerVectors, startPoint);
			
			checkRaysForCollisions(pixels, integerVectors, collisionPoints, w, h, d);
	    
			integerVectors = incrementVectors(integerVectors);
			
			i++;
		}
		
		final int nPoints = collisionPoints.size();
		int[][] result = new int[nPoints][3];
		
		Iterator<ArrayList<Integer>> iterator = collisionPoints.iterator();
		for (int j = 0; j < nPoints; j++) {
			ArrayList<Integer> list = iterator.next();
			for (int k = 0; k < 3; k++) {
				result[j][k] = list.get(k);
			}
		}
		
		IJ.log("Finished finding collision points after "+i+" iterations");
	  return result;
	}
	
/**
 * Set the current position to the last position plus integer vector
 * 
 * @param integerVectors
 * @return updated list of vectors
 */
	private HashSet<ArrayList<Double>> incrementVectors(HashSet<ArrayList<Double>> integerVectors) {
		
		HashSet<ArrayList<Double>> nextPosition = new HashSet<ArrayList<Double>>(integerVectors.size());
		ArrayList<Double> vector = new ArrayList<Double>(9);
				
		Iterator<ArrayList<Double>> iterator = integerVectors.iterator();
		while (iterator.hasNext()) {
			vector = iterator.next();
			
			//current location
			final Double x = vector.get(0);
			final Double y = vector.get(1);
			final Double z = vector.get(2);
			//vector
			final Double dx = vector.get(3);
			final Double dy = vector.get(4);
			final Double dz = vector.get(5);
			//type (face, edge, corner)
			final Double px = vector.get(6);
			final Double py = vector.get(7);
			final Double pz = vector.get(8);

			//set current position to integer vector plus last position
			final ArrayList<Double> shiftedVector = new ArrayList<Double>(
					Arrays.asList( x + dx, y + dy, z + dz, dx, dy, dz, px, py, pz )
						);
				
			if (IJ.debugMode)
			 IJ.log("Vector incremented from "+vector.toString()+
				 " to "+shiftedVector.toString());
			
			nextPosition.add(shiftedVector);
		}

		return nextPosition;
	}
	
	/**
	 * Bisect between existing rays within the current i plane
	 * 
	 * Do not spawn from corners, only from faces and edges
	 * 
	 * @param parentVectors
	 * @return Set of parent vectors and their child vectors
	 */
	@SuppressWarnings("unchecked")
	private HashSet<ArrayList<Double>> spawn(HashSet<ArrayList<Double>> parentVectors, int[] startPoint) {
		
		if (IJ.debugMode) {
			IJ.log("Spawning");
			IJ.log("parentVectors has size = "+parentVectors.size());
		}
			
		HashSet<ArrayList<Double>> childVectors = new HashSet<ArrayList<Double>>();
		
		Iterator<ArrayList<Double>> iterator = parentVectors.iterator();
		while (iterator.hasNext()) {
			ArrayList<Double> vector = iterator.next();
			
			if (isCorner(vector))
				continue;
			
			final Double zero = new Double(0);
			
			if (IJ.debugMode) {
				IJ.log("Parent vector is "+vector.toString());
			}
			//trim any imprecision to set ray to pixel grid
			vector.set(0, (double) Math.round(vector.get(0)));
			vector.set(1, (double) Math.round(vector.get(1)));
			vector.set(2, (double) Math.round(vector.get(2)));
			childVectors.add(vector);
			
			if (IJ.debugMode)
				IJ.log("Rounded vector is "+vector.toString());
			
			//create 4 child vectors as clones of the parent
			ArrayList<Double> child0 = (ArrayList<Double>) vector.clone();
			ArrayList<Double> child1 = (ArrayList<Double>) vector.clone();
			ArrayList<Double> child2 = (ArrayList<Double>) vector.clone();
			ArrayList<Double> child3 = (ArrayList<Double>) vector.clone();
			
			//get the directions to +- 0.5
			//if statement includes only faces
			//edges and corners are excluded
			//+- 0.5 in x & y
			if (vector.get(6).equals(zero) && vector.get(7).equals(zero)) {
				child0.set(0, child0.get(0) + 0.5);
				child1.set(0, child1.get(0) - 0.5);
				child2.set(1, child2.get(1) + 0.5);
				child3.set(1, child3.get(1) - 0.5);
			}
			//+- 0.5 in x & z
			else if (vector.get(6).equals(zero) && vector.get(8).equals(zero)) {
				child0.set(0, child0.get(0) + 0.5);
				child1.set(0, child1.get(0) - 0.5);
				child2.set(2, child2.get(2) + 0.5);
				child3.set(2, child3.get(2) - 0.5);
			}
			//+- 0.5 in y and z
			else if (vector.get(7).equals(zero) && vector.get(8).equals(zero)) {
				child0.set(1, child0.get(1) + 0.5);
				child1.set(1, child1.get(1) - 0.5);
				child2.set(2, child2.get(2) + 0.5);
				child3.set(2, child3.get(2) - 0.5);
			}
			//handle edges
			else if (isEdge(vector)) {
				if (IJ.debugMode)
					IJ.log("Is edge");
				if (vector.get(6).equals(zero)) {
					//spawn along edge
					child0.set(0, child0.get(0) + 0.5);
					child1.set(0, child1.get(0) - 0.5);
					//spawn onto face
					child2.set(1, child2.get(1) - vector.get(7) * 0.5);
					child3.set(2, child3.get(2) - vector.get(8) * 0.5);
					//set vector type to face
					child2.set(7, 0.0);
					child3.set(8, 0.0);
				}
				else if (vector.get(7).equals(zero)) {
					child0.set(1, child0.get(1) + 0.5);
					child1.set(1, child1.get(1) - 0.5);
					child2.set(0, child2.get(0) - vector.get(6) * 0.5);
					child3.set(2, child3.get(2) - vector.get(8) * 0.5);
					child2.set(6, 0.0);
					child3.set(8, 0.0);
				}
				else if (vector.get(8).equals(zero)) {
					child0.set(2, child0.get(2) + 0.5);
					child1.set(2, child1.get(2) - 0.5);
					child2.set(0, child2.get(0) - vector.get(6) * 0.5);
					child3.set(1, child3.get(1) - vector.get(7) * 0.5);
					child2.set(6, 0.0);
					child3.set(7, 0.0);
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
			
			if (IJ.debugMode) {
				IJ.log("   |----- child0 vector is "+child0.toString());
				IJ.log("   |----- child1 vector is "+child1.toString());
				IJ.log("   |----- child2 vector is "+child2.toString());
				IJ.log("   |----- child3 vector is "+child3.toString());
			}
		}
		
		if (IJ.debugMode)
			IJ.log("childVectors has size = "+childVectors.size());
		return childVectors;
	}


	/**
	 * Determine whether a vector is an edge
	 * Edges have one 0 and two +-1
	 * 
	 * Check that there is exactly 1 zero
	 *  
	 * @param vector
	 * @return true if this pixel is on an edge ray
	 */
	private boolean isEdge(ArrayList<Double> vector) {
		int sum = 0;
		for (int i = 6; i < 9; i++)
			if (Double.compare(vector.get(i), 0) == 0)
				sum++;
		
		return sum == 1;
	}

	/**
	 * Determine whether a vector is a corner
	 * Corners have zero 0 and three +-1
 
	 * @param vector
	 * @return true if this pixel is on a corner ray
	 */
	private boolean isCorner(ArrayList<Double> vector) {
		int sum = 0;
		for (int i = 6; i < 9; i++)
			if (Double.compare(vector.get(i), 0) == 0)
				sum++;
		
		return sum == 0;
	}

	/**
	 * Calculate the vector that when added to the current location
	 * results in sampling on the next integer pixel coordinate
	 * 
	 * @param vector
	 * @param startPoint
	 */
	private void calculateIntegerVector(ArrayList<Double> vector, int[] startPoint) {
		final Double zero = new Double(0);
		double l = 0;
		//in x & y plane, normal is z
		if (vector.get(6).equals(zero) && vector.get(7).equals(zero)) {
			l = vector.get(2) - startPoint[2];
		}
		//in x & z plane, normal is y
		else if (vector.get(6).equals(zero) && vector.get(8).equals(zero)) {
			l = vector.get(1) - startPoint[1];
		}
		//in y and z plane, normal is x
		else if (vector.get(7).equals(zero) && vector.get(8).equals(zero)) {
			l = vector.get(0) - startPoint[0];
		}
		//handle edges
		else if (isEdge(vector)) {
			//edge is parallel to x
			if (vector.get(6).equals(zero)) {
				l = vector.get(2) - startPoint[2];
			}
			//edge is parallel to y
			else if (vector.get(7).equals(zero)) {
				l = vector.get(0) - startPoint[0];
			}
			//edge is parallel to z
			else if (vector.get(8).equals(zero)) {
				l = vector.get(0) - startPoint[0];
			}
		}
		l = Math.abs(l);
		//x component
		vector.set(3, (vector.get(0) - startPoint[0])/l);
		
		//y component
		vector.set(4, (vector.get(1) - startPoint[1])/l);
		
		//z component
		vector.set(5, (vector.get(2) - startPoint[2])/l);
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
	private void checkRaysForCollisions(ByteProcessor[] pixels, HashSet<ArrayList<Double>> integerVectors,
		HashSet<ArrayList<Integer>> collisionPoints, final int w, final int h, final int d)
	{
		if (IJ.debugMode) {
		 IJ.log("Checking rays for collisions");
		 IJ.log("integerVectors has size = "+integerVectors.size());
		 IJ.log("collisionPoints has size = "+collisionPoints.size());
		}
		Iterator<ArrayList<Double>> iterator = integerVectors.iterator();
		while (iterator.hasNext()) {
			ArrayList<Double> vector = iterator.next();
			//round double values to snap to pixel grid
			final int x = (int) Math.round(vector.get(0));
			final int y = (int) Math.round(vector.get(1));
			final int z = (int) Math.round(vector.get(2));

			if (isOutOfBounds(x, y, z, w, h, d)) {
				iterator.remove();
				continue;
			}
			if (isBackground(pixels, x, y, z)) {
				iterator.remove();
				collisionPoints.add(new ArrayList<Integer>(Arrays.asList(x, y, z)));
				if (IJ.debugMode)
					IJ.log("Added a collision point at ("+x+", "+y+", "+z+")");
			}
		}
		if (IJ.debugMode) {
			IJ.log("Finished checking rays for collisions");
			IJ.log("integerVectors has size = "+integerVectors.size());
			IJ.log("collisionPoints has size = "+collisionPoints.size());
		}
	}

	/**
	 * Check whether a coordinate (x, y, z) is outside the 
	 * image stack
	 * @param x
	 * @param y
	 * @param z
	 * @param w
	 * @param h
	 * @param d
	 * @return true if (x, y, z) is out of image bounds
	 */
	private boolean isOutOfBounds(final int x, final int y, final int z, 
		final int w, final int h, final int d) {
		
		return x < 0 || x >= w || y < 0 || y >= h || z < 0 || z >= d;
	}


	/**
	 * Check whether the given point is background
	 * @param pixels image pixel data
	 * @param x x-coordinate of the test point
	 * @param y y-coordinate of the test point
	 * @param z z-coordinate of the test point
	 * @return true if the point is background
	 */
	private boolean isBackground(ByteProcessor[] pixels,
		final int x, final int y, final int z) {
		
		return pixels[z].get(x, y) == 0;
	}

	/**
	 * Check whether an integer is a power of 2 (2<sup>n</sup>)
	 * @param number
	 * @return true if number is an integer power of 2 
	 */
	private static boolean isPowerOfTwo(final int number) {
    return number > 0 && ((number & (number - 1)) == 0);
  }
}
