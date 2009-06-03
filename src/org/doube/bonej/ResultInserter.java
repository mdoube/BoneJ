package org.doube.bonej;

import ij.ImagePlus;
import ij.process.ImageProcessor;
import ij.measure.ResultsTable;
import ij.plugin.filter.PlugInFilter;

/**
 * Receive results from analytic methods and insert them into the 
 * Results table in a sensible way.
 * 
 * <p>Each image gets a line; measurements of different types 
 * are added to the same line; repeat measurements on same image
 * go on a new line.</p>
 * 
 * @author Michael Doube
 *
 */
public class ResultInserter implements PlugInFilter{
    @SuppressWarnings("unused")
    private ImagePlus imp;

    public int setup(String arg, ImagePlus imp){
	this.imp = imp;
	return DOES_ALL;
    }
    public void run(ImageProcessor ip){
	return;
    }
    /**
     * Finds the first available space for a result,
     * avoiding lots of empty space when measurements of different types
     * are made on the same image
     * 
     * @param imp ImagePlus
     * @param colHeading column heading
     * @param value value to insert
     */
    //TODO use a table other than the system Results table
    public void setResultInRow(ImagePlus imp, String colHeading, double value){
	ResultsTable rt = ResultsTable.getResultsTable();
	String title = imp.getTitle();
	String table = "Results";
	rt.show(table);

	//search for the first row that contains the image title
	//and contains no value for the heading
	for (int row = 0; row < rt.getCounter(); row++){
	    if (rt.getLabel(row).equals(title)){
		//there could be no column called colHeading
		if (!rt.columnExists(rt.getColumnIndex(colHeading))){
		    //in which case, just insert the value
		    rt.setValue(colHeading, row, value);
		    rt.show(table);
		    return;
		} else {
		    //but if there is, it might or might not have data in it
		    double currentValue =  rt.getValue(colHeading, row);
		    if(currentValue == 0){
			rt.setValue(colHeading, row, value);
			rt.show(table);
			return;
		    } else {
			//look for another row with the right title
			continue;
		    }
		}
	    }
	}
	//we got to the end of the table without finding a space to insert
	//the value, so make a new row for it
	rt.incrementCounter();
	rt.addLabel("Image", title);
	rt.addValue(colHeading, value);
	rt.show(table);
	return;
    }
}
