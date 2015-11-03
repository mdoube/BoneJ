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
	public static double[][] removeNulls(final double[][] d) {
		final int l = d.length;
		int nullCount = 0;
		for (int i = 0; i < l; i++)
			if (d[i] == null)
				nullCount++;
		if (nullCount == 0)
			return d;
		final int nonNulls = l - nullCount;
		final double[][] array = new double[nonNulls][];

		int j = 0;
		for (int i = 0; i < l; i++) {
			if (d[i] != null) {
				array[j] = d[i];
				j++;
			}
		}

		return array;
	}

	public static Ellipsoid[] removeNulls(final Ellipsoid[] o) {
		final int l = o.length;
		int nullCount = 0;
		for (int i = 0; i < l; i++)
			if (o[i] == null)
				nullCount++;
		if (nullCount == 0)
			return o;
		final int nonNulls = l - nullCount;
		final Ellipsoid[] array = new Ellipsoid[nonNulls];

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
	 * that row values become column values and vice versa
	 *
	 * @param d
	 * @return
	 * @throws IllegalArgumentException
	 *             if arrays are not all of the same length
	 */
	public static double[][] transpose(final double[][] d) {
		final int l = d.length;

		final double[][] t = new double[l][l];
		for (int i = 0; i < l; i++) {
			final double[] di = d[i];
			if (di.length != l)
				throw new IllegalArgumentException();
			for (int j = 0; j < l; j++)
				t[j][i] = di[j];
		}
		return t;
	}

	/**
	 * Matrix-free version of Matrix.times
	 *
	 * Multiplies a by b in the matrix multiplication scheme c = ab
	 *
	 * @param a
	 * @param b
	 * @return c
	 */
	public static double[][] times(final double[][] a, final double[][] b) {
		final int am = a.length;
		final int an = a[0].length;
		final int bm = b.length;
		final int bn = b[0].length;
		if (bm != an) {
			throw new IllegalArgumentException("Matrix inner dimensions must agree.");
		}
		final double[][] c = new double[am][bn];
		final double[] bcolj = new double[an];
		for (int j = 0; j < bn; j++) {
			for (int k = 0; k < an; k++) {
				bcolj[k] = b[k][j];
			}
			for (int i = 0; i < am; i++) {
				final double[] arowi = a[i];
				double s = 0;
				for (int k = 0; k < an; k++) {
					s += arowi[k] * bcolj[k];
				}
				c[i][j] = s;
			}
		}
		return c;
	}
}
