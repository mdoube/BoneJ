package org.doube.util;

/**
 * ResultInserter plugin for ImageJ
 * Copyright 2009 2010 Michael Doube 
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

import ij.ImagePlus;
import ij.measure.ResultsTable;

/**
 * Receive results from analytic methods and insert them into the Results table
 * in a sensible way.
 * 
 * <p>
 * Each image gets a line; measurements of different types are added to the same
 * line; repeat measurements on same image go on a new line.
 * </p>
 * 
 * @author Michael Doube
 * 
 */
public class ResultInserter {
	public static final ResultInserter INSTANCE = new ResultInserter();
	private static ResultsTable rt;

	private ResultInserter() {
	}

	public static ResultInserter getInstance() {
		rt = ResultsTable.getResultsTable();
		String table = "Results";
		rt.show(table);
		return INSTANCE;
	}

	/**
	 * Finds the first available space for a result, avoiding lots of empty
	 * space when measurements of different types are made on the same image
	 * 
	 * @param imp
	 *            ImagePlus
	 * @param colHeading
	 *            column heading
	 * @param value
	 *            value to insert
	 */
	// TODO use a table other than the system Results table
	public void setResultInRow(ImagePlus imp, String colHeading, double value) {
		String title = imp.getTitle();

		// search for the first row that contains the image title
		// and contains no value for the heading
		for (int row = 0; row < rt.getCounter(); row++) {
			if (rt.getLabel(row) == null){
				rt.setLabel(title, row);
			}
			if (rt.getLabel(row).equals(title)) {
				// there could be no column called colHeading
				if (!rt.columnExists(rt.getColumnIndex(colHeading))) {
					// in which case, just insert the value
					rt.setValue(colHeading, row, value);
					// and fill all other newly created cells with NaN
					for (int r = 0; r < rt.getCounter(); r++) {
						if (r != row)
							rt.setValue(colHeading, r, Double.NaN);
					}
					return;
				} else {
					// but if there is, it might or might not have data in it
					Double currentValue = rt.getValue(colHeading, row);
					if (currentValue.equals(Double.NaN)) {
						rt.setValue(colHeading, row, value);
						return;
					} else {
						// look for another row with the right title
						continue;
					}
				}
			}
		}
		// we got to the end of the table without finding a space to insert
		// the value, so make a new row for it
		String label = "Image";
		rt.incrementCounter();
		rt.addLabel(label, title);
		rt.addValue(colHeading, value);
		// set all the other values in the row to Double.NaN
		// (IJ has set them to 0, which is unhelpful, and we can't use null)
		int row = rt.getCounter() - 1;
		for (int c = 0; c <= rt.getLastColumn(); c++) {
			if (!rt.getColumnHeading(c).equals(colHeading)
					&& !rt.getColumnHeading(c).equals(label)) {
				rt.setValue(c, row, Double.NaN);
			}
		}
		return;
	}

	/**
	 * Show the table
	 */
	public void updateTable() {
		String table = "Results";
		rt.show(table);
	}
}
