/*
	This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

	N.B.  the above text was copied from http://www.gnu.org/licenses/gpl.html
	unmodified. I have not attached a copy of the GNU license to the source...

    Copyright (C) 2011 Timo Rantalainen
*/

package org.doube.bonej.pqct.selectroi;
import java.util.*;	//Vector, Collections
import java.lang.Math; //atan2
import java.awt.*;			//Polygon, Rectangle
import org.doube.bonej.pqct.io.*;	//image data
import ij.*;		//ImagePlus
import ij.gui.*;	//ImagePlus ROI
import ij.text.*; 	//Debugging ...
import ij.process.*;	//Debugging
@SuppressWarnings(value ={"serial","unchecked"}) //Unchecked for obtaining Vector<Object> as a returnvalue

public class SelectROI extends RoiSelector{
	//ImageJ constructor
	public SelectROI(ScaledImageData dataIn,ImageAndAnalysisDetails detailsIn, ImagePlus imp,double boneThreshold,boolean setRoi){
		super(dataIn,detailsIn, imp,boneThreshold,setRoi);
		//Select ROI
		
		/*Select ROI and set everything else than the roi to minimum*/
		cortexROI = new double[width*height];	//Make a new copy of the image with only the ROI remaining
		roiI = new Vector<Integer>();
		roiJ = new Vector<Integer>();
		cortexRoiI = new Vector<Integer>();
		cortexRoiJ = new Vector<Integer>();
		cortexAreaRoiI = new Vector<Integer>();
		cortexAreaRoiJ = new Vector<Integer>();
		boneMarrowRoiI = new Vector<Integer>();
		boneMarrowRoiJ = new Vector<Integer>();
		Roi ijROI = imp.getRoi();
		double[] tempScaledImage = (double[]) scaledImage.clone();
		if (ijROI != null && details.manualRoi){	/*Set pixels outside the manually selected ROI to zero*/
			/*Check whether pixel is within ROI, mark with bone threshold*/
			for (int j = 0;j< height;j++){
				for (int i = 0; i < width;i++){
					if (ijROI.contains(i,j)){
					}else{
						tempScaledImage[i+j*width] = minimum;
					}
				}
			}
			/*Check whether a polygon can be acquired and include polygon points too*/
			Polygon polygon = ijROI.getPolygon();
			if (polygon != null){
				for (int j = 0;j< polygon.npoints;j++){
					tempScaledImage[polygon.xpoints[j]+polygon.ypoints[j]*width] = scaledImage[polygon.xpoints[j]+polygon.ypoints[j]*width];
				}
			}
		}
		length					= new Vector<Integer> ();
		beginnings				= new Vector<Integer> ();
		iit						= new Vector<Integer> ();
		jiit					= new Vector<Integer> ();
		area	= new Vector<Integer> ();
		result 					= new byte[width*height];
		Vector<Object> boneMasks = getSieve(tempScaledImage,result,area,length,beginnings, iit, jiit,roiI,roiJ,boneThreshold,details.roiChoice,details.guessStacked,details.stacked,details.guessFlip,details.allowCleaving);
		sieve			= (byte[]) boneMasks.get(0);
		result	 		= (byte[]) boneMasks.get(1);
		iit 		 	= (Vector<Integer>) boneMasks.get(2);
		jiit 			= (Vector<Integer>) boneMasks.get(3);
		beginnings		= (Vector<Integer>) boneMasks.get(4);
		length			= (Vector<Integer>) boneMasks.get(5);
		area			= (Vector<Integer>) boneMasks.get(6);
		selection		= (Integer)	 boneMasks.get(7);
		/*Add the roi to the image*/
		if (setRoi){
			int[] xcoordinates = new int[roiI.size()];
			int[] ycoordinates = new int[roiJ.size()];
			for (int i = 0;i<roiI.size();++i){
				xcoordinates[i] = roiI.get(i);
				ycoordinates[i] = roiJ.get(i);
			}
			ijROI = new PolygonRoi(xcoordinates,ycoordinates,roiI.size(),Roi.POLYGON);
			imp.setRoi(ijROI);
		}
		
		for (int j = 0;j< height;j++){
			for (int i = 0; i < width;i++){
				if (scaledImage[i+j*width]<marrowThreshold && sieve[i+j*width] > 0){
					boneMarrowRoiI.add(i);
					boneMarrowRoiJ.add(j);
				}
				if (scaledImage[i+j*width]>=areaThreshold && sieve[i+j*width] > 0){
					cortexAreaRoiI.add(i);
					cortexAreaRoiJ.add(j);
				}
				if (scaledImage[i+j*width]>=BMDthreshold && sieve[i+j*width] > 0){
					cortexROI[i+j*width] = scaledImage[i+j*width];				
					cortexRoiI.add(i);
					cortexRoiJ.add(j);
				} else {
					cortexROI[i+j*width] = minimum;
				}
			}
		}
	}	
}
