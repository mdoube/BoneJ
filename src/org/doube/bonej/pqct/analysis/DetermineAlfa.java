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

package org.doube.bonej.pqct.analysis;
import java.util.*;	//Vector, Collections
import java.lang.Math; //atan2
import org.doube.bonej.pqct.selectroi.*;	//ROI selection..
import org.doube.bonej.pqct.io.*;

public class DetermineAlfa{
	public int rotationIndex;
	public double alfa = 0;
	public double rotationCorrection = 0;
	public Vector<Integer> pind;
	public Vector<Integer> pindColor;
	ImageAndAnalysisDetails details;
	public DetermineAlfa(SelectROI roi,ImageAndAnalysisDetails details){
		this.details = details;
		//Calculate CSMIs and rotation angle to align maximal and minimal bending axes with X and Y axes
		rotationCorrection = (((double) details.sectorWidth)/2.0); 
		/*Rotation according to Imax/Imin for bone of interest or according to all bones*/
		if (details.rotationChoice.equals("All_Bones_Imax/Imin") || details.rotationChoice.equals("According_to_Imax/Imin")){
			double[] csmiValues = new double[3];
			if (details.rotationChoice.equals("All_Bones_Imax/Imin")){
				csmiValues = csmi(roi.result,roi.width,roi.height);
			}
			if (details.rotationChoice.equals("According_to_Imax/Imin")){
				byte[] tempCsmiSieve = new byte[roi.width*roi.height];
				for (int j = 0; j< roi.height;j++){
					for (int i = 0; i<roi.width;i++){
						if (roi.cortexROI[i+j*roi.width] >= details.areaThreshold){
							tempCsmiSieve[i+j*roi.width] = 1;
						} else {
							tempCsmiSieve[i+j*roi.width] = 0;
						}
					}
				}
				csmiValues = csmi(tempCsmiSieve,roi.width,roi.height);
			}
			double xmax = csmiValues[0];
			double ymax = csmiValues[1];
			double moment = csmiValues[2];
			double vali1,vali2;		
			//Calculate rotation required to align rotation axes
			if (ymax == xmax){
				alfa = 0;	//check that xmax does not equal ymax (can't divide with 0...			
			}else{
				alfa = Math.atan(2.0*moment/(ymax-xmax))/2.0;
				//Calculate the maximal and minimial cross-sectional moments of inertia
				vali1 = (ymax+xmax)/2+(ymax-xmax)/2*Math.cos(2*(-alfa))-moment*Math.sin(2*(-alfa));
				vali2 =(ymax+xmax)/2-(ymax-xmax)/2*Math.cos(2*(-alfa))+moment*Math.sin(2*(-alfa));
				//The according to Imax/Imin alfa may align rotation axis corresponding to maximal CSMI with either horizontal 
				//or vertical axis, whichever rotation is smaller...
				//Always rotate towards horizontal axis... maximal bending axis will be aligned with horizontal axis
				//Note that e.g. tibial mid-shaft rotation is completely different if only tibia or if both tibia and fibula
				//are consireder!!!
				if (vali1 > vali2){
					if (alfa < 0) {
						alfa =Math.PI/2.0+alfa;
					}else{
						alfa =alfa-Math.PI/2.0;
					}
				}
				
			}
		}

		/*Rotation according to the furthest point*/
		if (details.rotationChoice.equals("Furthest_point")){
			/*Calculate alfa from periosteal radii*/
			double[] marrowCenter = new double[2];
			for (int i = 0; i< roi.boneMarrowRoiI.size();i++){
				marrowCenter[0]+=(double)roi.boneMarrowRoiI.get(i);
				marrowCenter[1]+=(double)roi.boneMarrowRoiJ.get(i);
			}
			marrowCenter[0] /=(double)roi.boneMarrowRoiI.size();
			marrowCenter[1] /=(double)roi.boneMarrowRoiJ.size();
			double[] radii = new double[roi.roiI.size()];
			for (int i = 0; i<roi.roiI.size();++i){
				radii[i] = Math.sqrt(Math.pow(roi.roiI.get(i)-marrowCenter[0],2)+Math.pow(roi.roiJ.get(i)-marrowCenter[1],2));
			}
			double[] sumRadii = new double[radii.length];
			for (int i = 5;i<radii.length-6;++i){
				for (int j = -5;j<6;++j){
					sumRadii[i]+=radii[i+j];
				}
			}
			double[] sortRadii = (double[]) sumRadii.clone();
			Arrays.sort(sortRadii);
			int largest =0;
			while (sumRadii[largest] != sortRadii[sortRadii.length-1]){
				++largest;
			}
			double x,y;
			x = roi.roiI.get(largest)-marrowCenter[0];
			y = roi.roiJ.get(largest)-marrowCenter[1];
			alfa = Math.PI-Math.atan2(y,x);
		}

		/*Manual rotation*/
		if (details.manualRotation){
			alfa = details.manualAlfa;
		}
		
		/*Flip distribution*/
		if (details.flipDistribution){
			rotationCorrection = -rotationCorrection;
		}
		
		rotationIndex = (int) (alfa/Math.PI*180.0+rotationCorrection);
		
		//Calculate CSMIs and rotation angle to align maximal and minimal bending axes with X and Y axes
		pind = rotateIndex(rotationIndex);
		
		if(details.flipDistribution){
			pindColor = rotateIndex((int) (rotationIndex));
		}else{
			pindColor = rotateIndex(-rotationIndex);
		}

	}
	
	Vector<Integer> rotateIndex(int rotationAngle){
		int initialIndex = 0;
		Vector<Integer> rotateIndexVector = new Vector<Integer>();
		if (rotationAngle >= 0){
			initialIndex = 360-rotationAngle; 
		}else{
			initialIndex = -rotationAngle;
		}
		int inde;
		inde = initialIndex;
		while (inde<360){
			rotateIndexVector.add(inde);
			++inde;
		}
		inde=0;
		while (inde < initialIndex){
			rotateIndexVector.add(inde);
			++inde;
		}
		
		/*Flip rotateIndexVector, for e.g. comparing left to right*/
		if (details.flipDistribution){
			Collections.reverse(rotateIndexVector);
		}
		return rotateIndexVector;
	}
	
	double[] csmi(byte[] sieve,int width, int height){
		double[] cortexCenter = new double[2];
		int points = 0;
		Vector<Integer> bmcI = new Vector<Integer>();
		Vector<Integer> bmcJ = new Vector<Integer>();
		for (int j = 0; j< height;j++){
			for (int i = 0; i< width;i++){
				if (sieve[i+j*width] > 0){
					cortexCenter[0]+=(double) i;
					cortexCenter[1]+=(double) j;
					bmcI.add(i);
					bmcJ.add(j);
					++points;
				}
			}
		}
		cortexCenter[0] /=(double)points;
		cortexCenter[1] /=(double)points;

		double[] returnValues = new double[3];
		for (int i = 0;i<returnValues.length;++i){returnValues[i] = 0;}
		//Calculating cross-sectional moment of inertia in the original image orientation
		for (int i = 0;i<bmcI.size();i++){
			returnValues[0] +=((bmcI.get(i)-cortexCenter[0]))*((bmcI.get(i)-cortexCenter[0]));
			returnValues[1] +=((bmcJ.get(i)-cortexCenter[1]))*((bmcJ.get(i)-cortexCenter[1]));
			returnValues[2] +=((bmcI.get(i)-cortexCenter[0]))*((bmcJ.get(i)-cortexCenter[1]));
		}		
		return returnValues;
	}

}