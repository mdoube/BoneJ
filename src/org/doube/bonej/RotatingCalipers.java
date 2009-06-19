package org.doube.bonej;

/**
 * RotatingCalipers plugin for ImageJ
 * Copyright 2009 Michael Doube 
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


import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ImageProcessor;
import ij.plugin.filter.PlugInFilter;
import ij.measure.Calibration;
import ij.measure.ResultsTable;
import ij.gui.*;

import org.doube.bonej.ResultInserter;

public class RotatingCalipers implements PlugInFilter {
    ImagePlus imp;
    protected ImageStack stack;
    public static final double PI = 3.141592653589793;
    public int setup(String arg, ImagePlus imp) {
	stack = imp.getStack();
	this.imp = imp;
	return DOES_ALL+ROI_REQUIRED+NO_CHANGES;
    }

    public void run(ImageProcessor ip) {
	double dMin = rotatingCalipers(this.imp);
	if (dMin > 0){
	    showResults(this.imp, dMin);
	} else {return;}
    }
    public double rotatingCalipers(ImagePlus imp){
//	double dmax = 0;
	double dMin = 0;
	IJ.run("Convex Hull");
	Roi roi = imp.getRoi();
	int perimx[] = ((PolygonRoi)roi).getXCoordinates();
	int perimy[] = ((PolygonRoi)roi).getYCoordinates();
	int podal = 0;
	int antipodal = 0;
	//	Compute the polygon's extreme points in the y direction
	int ymin = Integer.MAX_VALUE;
	int ymax = 0;
	int nPoints = perimx.length;
	for (int n = 0; n < nPoints; n++){
	    if (perimy[n] < ymin){
		ymin = perimy[n];
		podal = n;    		//n is always 0 in ImageJ, start point of convex hull-ed selection is minimum x value within set with minimum y
	    }
	    if (perimy[n] > ymax){
		ymax = perimy[n];
		antipodal = n;
	    }
	}

	//	create list of angles beween points
	double angles[]; angles = new double[nPoints];
	for (int n = 0; n < nPoints - 1; n++){
	    angles[n] = Math.atan2(perimy[n]-perimy[n+1] , perimx[n]-perimx[n+1]);
	}

	//	angle between last point and first point of selection
	angles[nPoints-1] = Math.atan2(perimy[nPoints-1]-perimy[0] , perimx[nPoints-1]-perimx[0]);

	//	correct any angles < -pi
	for (int n=0; n < nPoints; n++) {
	    if (angles[n] < -PI) angles[n] = angles[n] + 2*PI;
	}

	//	Construct two horizontal lines of support through ymin and ymax.
	//	Since this is already an anti-podal pair, compute the distance, and keep as maximum.
	dMin = ymax - ymin;
	int end = antipodal;
	double startTheta = 0; double endTheta = 0; double testAngle = 0; 
	while (antipodal < nPoints  && podal <= end){
	    int podincr = 0;
	    int antipodincr = 0;
	    //	find the start condition for theta
	    if (podal > 0) {
		startTheta = endTheta;
	    } else {
		if (angles[antipodal - 1] > 0) { 
		    testAngle = angles[antipodal - 1] - PI;
		} else {
		    testAngle = angles[antipodal - 1] + PI;
		}
		if (angles[nPoints -1] <= testAngle){
		    startTheta = angles[nPoints - 1];
		} else {
		    startTheta = testAngle;
		}
	    }

	    //	Rotate the lines until one is flush with an edge of the polygon.
	    //	find the end condition for theta
	    //	increment the podal or antipodal point depending on
	    //	which caliper blade touches the next point first

	    if (angles[antipodal] > 0) { 
		testAngle = angles[antipodal] - PI;
	    } else {
		testAngle = angles[antipodal] + PI;
	    }
	    if (antipodal <= nPoints - 1) {
		if (angles[podal] > testAngle && angles[podal] * testAngle > 0) {
		    endTheta = angles[podal];
		    podincr = 1;
		} else if (angles[podal] < testAngle) {
		    endTheta = testAngle;
		    antipodincr = 1;
		} else if (angles[podal] > testAngle && angles[podal] * testAngle < 0) {
		    endTheta = testAngle;
		    antipodincr = 1;
		} else {
		    endTheta = angles[podal];
		    podincr = 1;
		    antipodincr = 1;
		}
	    } else {
		endTheta = 0;
	    }

	    //	A new anti-podal pair is determined. Compute the new distance,
	    //	compare to old maximum, and update if necessary.
	    double thetaH = Math.atan2(perimy[antipodal] - perimy[podal] , perimx[antipodal] - perimx[podal]);
	    double dp = Math.sqrt(Math.pow(perimx[antipodal]-perimx[podal],2) + Math.pow(perimy[antipodal]-perimy[podal],2));
	    if (startTheta >= endTheta){
		double ds = dp*Math.abs(Math.cos(-thetaH + startTheta - PI/2));
		double de = dp*Math.abs(Math.cos(-thetaH + endTheta - PI/2));
		if (ds < dMin) dMin = ds;
		else if (de < dMin) dMin = de;
	    }	
	    else if (startTheta < endTheta && startTheta * endTheta < 0){
		double ds = dp*Math.abs(Math.cos(-thetaH + startTheta - PI/2));
		double de = dp*Math.abs(Math.cos(-thetaH + - 3*PI/2));
		if (ds < dMin) dMin = ds;
		else if (de < dMin) dMin = de;
		ds = dp*Math.abs(Math.cos(-thetaH + PI/2));
		de = dp*Math.abs(Math.cos(-thetaH + endTheta - PI/2));
		if (ds < dMin) dMin = ds;
		else if (de < dMin) dMin = de;
	    }
	    //	increment either or both podal or antipodal
	    podal = podal + podincr;
	    antipodal = antipodal + antipodincr;
	}
	return dMin;
    }
    
    public void showResults(ImagePlus imp, double dMin){
	ResultInserter ri = new ResultInserter();
	ri.setResultInRow(imp, "RC", dMin);
    }
}