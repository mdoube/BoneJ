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

public class SelectSoftROI extends RoiSelector{
	//ImageJ constructor
	public SelectSoftROI(ScaledImageData dataIn,ImageAndAnalysisDetails detailsIn, ImagePlus imp,double boneThreshold,boolean setRoi){
		this.scaledImageData = dataIn;
		this.imp = imp;
		details =detailsIn;
		scaledImage = (double[])dataIn.scaledImage.clone();
		softScaledImage = (double[])dataIn.softScaledImage.clone();
		pixelSpacing = dataIn.pixelSpacing;
		imageSavePath = details.imageSavePath;
		width =dataIn.width;
		height =dataIn.height;
		

		airThreshold = details.airThreshold;
		fatThreshold = details.fatThreshold;
		rotationThreshold = details.rotationThreshold;
		muscleThreshold = details.muscleThreshold;
		marrowThreshold = details.marrowThreshold;
		areaThreshold = details.areaThreshold;	//For cortical AREA analyses (CoA, SSI, I) + peeling distal pixels
		BMDthreshold = details.BMDthreshold;		//For cortical BMD analyses
		softThreshold = details.softThreshold;	//Thresholding soft tissues + marrow from bone
		this.boneThreshold = boneThreshold;
		minimum = dataIn.minimum;
		maximum = dataIn.maximum;
		/*A special function to check whether forearm has been measured palm up*/
		/*
		if (details.allowVerticalFlip){
			checkVerticalFlip();
		}
		*/
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
		
		//Soft tissue analysis
		softSieve = null;
		byte[] softResult = null;
		if (details.stOn){
			Vector<Integer> stLength		= new Vector<Integer>();
			Vector<Integer> stBeginnings	= new Vector<Integer>();
			Vector<Integer> stIit			= new Vector<Integer>();
			Vector<Integer> stJiit			= new Vector<Integer>();
			Vector<Integer> stRoiI			= new Vector<Integer>();
			Vector<Integer> stRoiJ			= new Vector<Integer>();
			Vector<Integer> stArea		= new Vector<Integer> ();
			softResult = new byte[width*height];
			
			
			/*Get rid of measurement tube used at the UKK institute*/
			byte[] sleeve = null;
			if (details.sleeveOn){
				sleeve = removeSleeve(softScaledImage,sleeve,25.0);
				int removed=0;
				for (int ii =0;ii<width*height;ii++){
					if(sleeve[ii]==1){
						softScaledImage[ii]=minimum;
						++removed;
					}
				}
			}
			
			Vector<Object> masks = getSieve(softScaledImage,softResult,stArea,stLength,stBeginnings, stIit, stJiit,stRoiI,stRoiJ,airThreshold,details.roiChoiceSt,details.guessStacked,details.stacked,false,true);
			softSieve		= (byte[]) masks.get(0);
			softResult	 	= (byte[]) masks.get(1);
			stIit 		 	= (Vector<Integer>) masks.get(2);
			stJiit 			= (Vector<Integer>) masks.get(3);
			stBeginnings	= (Vector<Integer>) masks.get(4);
			stLength		= (Vector<Integer>) masks.get(5);
			stArea			= (Vector<Integer>) masks.get(6);
			
			/*Erode three layers of pixels from the fat sieve to get rid of higher density layer (i.e. skin) 
			on top of fat to enable finding muscle border
			*/
			byte[] muscleSieve = (byte[]) softSieve.clone();
			double[] muscleImage = (double[]) softScaledImage.clone();
			byte[] subCutaneousFat = null;
			for (int i = 0;i< 3;++i){
				muscleSieve = erode(muscleSieve);
				if (i == 0){
					/*Add subcut fat sieve... Skin has already been removed by eroding one layer of pixels-> remove muscle later on*/
					subCutaneousFat = (byte[]) muscleSieve.clone();
				}
			}
			
			/*Remove everything other than the selected limb from the image*/
			for (int i = 0; i<muscleSieve.length;++i){
				if (muscleSieve[i] < 1){
					muscleImage[i] = minimum;
				}
			}
			/*Look for muscle outline*/
			Vector<Object> muscleMasks = getSieve(muscleImage,new byte[width*height],new Vector<Integer>(),new Vector<Integer>(),new Vector<Integer>(), new Vector<Integer>(), new Vector<Integer>(),new Vector<Integer>(),new Vector<Integer>(),details.muscleThreshold,"Bigger",details.guessStacked,details.stacked,false,false);
			muscleSieve		= (byte[]) muscleMasks.get(0);
			
			/*Wipe muscle area +1 layer of pixels away from subcut.*/
			byte[] tempMuscleSieve = (byte[]) muscleSieve.clone();
			dilate(tempMuscleSieve,(byte)1,(byte)0,(byte)2);
			for (int i = 0;i<tempMuscleSieve.length;++i){
				if (tempMuscleSieve[i] == 1){subCutaneousFat[i] = 0;}
			}
			/*create temp boneResult to wipe out bone and marrow*/
			byte[] boneResult = new byte[width*height];
			Vector<Object> masks2 = getSieve(softScaledImage,boneResult,new Vector<Integer>(),new Vector<Integer>(),new Vector<Integer>(), new Vector<Integer>(), new Vector<Integer>(),new Vector<Integer>(),new Vector<Integer>(),softThreshold,details.roiChoiceSt,details.guessStacked,details.stacked,false,false);
			boneResult	= (byte[]) masks2.get(1);
			for (int i = 0;i<softSieve.length;++i){
				if (softSieve[i] ==1 && softScaledImage[i] >= airThreshold && softScaledImage[i] < fatThreshold){
					softSieve[i] =2;	//Fat
				}
				if (muscleSieve[i] ==1 && softScaledImage[i] >= muscleThreshold && softScaledImage[i] < softThreshold){
					softSieve[i] = 3;	//Muscle
				}
				if (muscleSieve[i] ==1 && softScaledImage[i] >= airThreshold && softScaledImage[i] < muscleThreshold){
					softSieve[i] = 4;	//Intra/Intermuscular fat
				}
				if (subCutaneousFat[i] ==1 ){
					softSieve[i] = 5;	//Subcut fat
				}
				if (boneResult[i] ==1 ){
					softSieve[i] = 6;	//Bone & marrow
				}
			}
		}
	}
}
