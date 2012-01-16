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
		super(dataIn,detailsIn, imp,boneThreshold,setRoi);
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
