package org.doube.util;

import org.doube.geometry.Ellipsoid;

/**
 * Methods for doing helpful things with arrays
 * 
 * @author Michael Doube
 * 
 */
public class ArrayHelper {

	/**
	 * Remove null values from an array
	 * 
	 * @param d
	 * @return array containing only non-null values. Element ordering remains
	 *         intact, null elements are replaced with no element, so the
	 *         resulting array is the length of the input array minus the number
	 *         of null elements.
	 */
	public static double[][] removeNulls(double[][] d) {
		final int l = d.length;
		int nullCount = 0;
		for (int i = 0; i < l; i++)
			if (d[i] == null)
				nullCount++;
		if (nullCount == 0)
			return d;
		final int nonNulls = l - nullCount;
		double[][] array = new double[nonNulls][];

		int j = 0;
		for (int i = 0; i < l; i++) {
			if (d[i] != null) {
				array[j] = d[i];
				j++;
			}
		}

		return array;
	}

	/**
	 * Remove null values from an array of Objects
	 * 
	 * @param o
	 *            1D Object[] array
	 * @return Array containing same elements in same order but without null
	 *         values, with length = o.length() - number of null values
	 */
	public static Object[] removeNulls(Object[] o) {
		final int l = o.length;
		int nullCount = 0;
		for (int i = 0; i < l; i++)
			if (o[i] == null)
				nullCount++;
		if (nullCount == 0)
			return o;
		final int nonNulls = l - nullCount;
		Object[] array = new Object[nonNulls];

		int j = 0;
		for (int i = 0; i < l; i++) {
			if (o[i] != null) {
				array[j] = o[i];
				j++;
			}
		}
		return array;
	}

	/**
	 * Transpose a square (rows and columns all the same length) double array so
	 * that row values become calumn values and vice versa
	 * 
	 * @param d
	 * @return
	 * @throws IllegalArgumentException
	 *             if arrays are not all of the same length
	 */
	public static double[][] transpose(double[][] d) {
		final int l = d.length;

		double[][] t = new double[l][l];
		for (int i = 0; i < l; i++) {
			if (d[i].length != l)
				throw new IllegalArgumentException();
			for (int j = 0; j < l; j++)
				t[i][j] = d[j][i];
		}
		return t;
	}

}
