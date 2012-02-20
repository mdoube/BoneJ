package org.doube.bonej.pqct.selectroi;
import java.util.*;	//Vector, Collections
import java.lang.Math; //atan2
import org.doube.bonej.pqct.io.*;	//image data
import ij.*;		//ImagePlus
public class SoftTissueSide{
	public SelectROI roi;
	public byte[] softTissueSideSieve;
	public byte softTissueSide;	/*2 = left, 1 = right*/
	public ImageAndAnalysisDetails details;
	public boolean changeSign;
	public SoftTissueSide(SelectROI roi,ImageAndAnalysisDetails details, boolean changeSign){
		this.roi = roi;
		this.details = details;
		this.changeSign = changeSign;
		softTissueSideSieve = fillSides(roi);	/*Sieve filled in with soft tissues split to left and right*/
		softTissueSide = checkSide(softTissueSideSieve, roi.width,roi.height);
		softTissueSideSieve = removeUnselected(softTissueSideSieve,softTissueSide);
		
	}
	
	byte[] removeUnselected(byte[] softTissueSideSieve,byte softTissueSide){
		for (int i = 0; i<softTissueSideSieve.length;++i){
			if (softTissueSideSieve[i] != softTissueSide && softTissueSideSieve[i] != 0 && softTissueSideSieve[i] < 3){
				softTissueSideSieve[i] = 0;
			}
		}
		return softTissueSideSieve;
	}
	
	byte checkSide(byte[] softTissueSideSieve, int width,int height){
		byte side;
		int[] whichSide = new int[2];
		for (int h = 0;h<height;++h){
			for (int w = 0;w<width;++w){
				if (softTissueSideSieve[w+h*width] == 1) whichSide[1]++;
				if (softTissueSideSieve[w+h*width] == 2) whichSide[0]++;
			}
		}
		side = (whichSide[0] > whichSide[1]) ? (byte) 2 : (byte) 1;
		return side;
	}
	
	/*
		Divide the soft tissue into two halves by the line, which connects the bone centres
		Find the rotation angle to vertical
		Check on which side a given pixel resides in, taking the rotation into account
	
	*/
	byte[] fillSides(SelectROI roi){
		
		int[] consideredBones = roi.twoLargestBonesDetectedEdges(roi.edges);	/*Get the two bones*/
		double[][] boneCoordinates = getCoordinates(roi,consideredBones);	/*Fill a sieve with the two bones*/
		double rotationAngle = getRotation(boneCoordinates);	/*Find the rotation angle to vertical*/
		/*	The rotation rotates the line connecting the bones to vertical */
		byte[] sieve = new byte[roi.width*roi.height];
		sieve[(int)boneCoordinates[0][0]+(int)boneCoordinates[0][1]*roi.width] = 1;
		double[] tempImage = (double[]) roi.scaledImage.clone();
		if (details.sleeveOn){
			/*Get rid of measurement tube used at the UKK institute*/
			byte[] sleeve = null;
			if (details.sleeveOn){
				sleeve = roi.removeSleeve(tempImage,sleeve,25.0);
				for (int ii =0;ii<roi.width*roi.height;ii++){
					if(sleeve[ii]==1){
						tempImage[ii]=roi.minimum;
					}
				}
			}		
		}
		/*Dilate sieve, into neighbouring pixels, until air is found*/
		int tempDil = 1;
		while (tempDil>0){
			tempDil=roi.dilateLimb(sieve,(byte)1,(byte)0,(byte)4,roi.details.airThreshold,tempImage);
		}
		byte[] returnSieve = fillSieveSides(roi,rotationAngle,sieve,boneCoordinates);
		
		return returnSieve;
	}
	
	/*Fill sieve sides*/
	byte[] fillSieveSides(SelectROI roi, double rotationAngle,byte[] sieve, double[][] boneCoordinates){
		if (details.manualRotation){
			rotationAngle = details.manualAlfa;
		}
		if (changeSign) {rotationAngle = -rotationAngle;}
		for (int h = 0;h<roi.height;++h){
			for (int w = 0;w<roi.width;++w){
				if(w*Math.cos(rotationAngle)-h*Math.sin(rotationAngle) < boneCoordinates[0][0]*Math.cos(rotationAngle)-boneCoordinates[0][1]*Math.sin(rotationAngle) &&
					sieve[w+h*roi.width] > 0
					){
					sieve[w+h*roi.width] = 2;
				}
			}
		}
		
		for (int b = 0; b<2 ;++b){
			for (int h = (int)boneCoordinates[b][1]-1;h<=(int)boneCoordinates[b][1]+1;++h){
				for (int w = (int)boneCoordinates[b][0]-1;w<=(int)boneCoordinates[b][0]+1;++w){
					sieve[w+h*roi.width] = (byte)(3+b);		
				}
			}
		}
		
		return sieve;
	}
	
	
	/*Find the rotation angle to vertical*/
	double getRotation(double[][] boneCoordinates){
		double x = 0;
		double y = 0;
		double[] alfa = new double[2];
		/*Search for minimal rotation, check both from larger to smaller and from smaller to larger!!*/
		/*From larger to smaller -> boneCoorinate[0] is larger bone*/
		x = boneCoordinates[1][0]-boneCoordinates[0][0];	//Use the selected bone as origin for rotation
		y = boneCoordinates[1][1]-boneCoordinates[0][1];	//Use the selected bone as origin for rotation
		alfa[0]= Math.PI/2.0-Math.atan2(y,x);
		/*From smaller to larger*/
		x = boneCoordinates[0][0]-boneCoordinates[1][0];	//Use the selected bone as origin for rotation
		y = boneCoordinates[0][1]-boneCoordinates[1][1];	//Use the selected bone as origin for rotation
		alfa[1]= Math.PI/2.0-Math.atan2(y,x);
		double rotation = (Math.abs(alfa[0]) < Math.abs(alfa[1])) ? alfa[0] : alfa[1];
		System.out.println("R angle "+rotation);
		return rotation;	
	}
	/*Fill a sieve with the two bones*/
	double[][] getCoordinates(SelectROI tempRoi, int[] consideredBones){
		double[][] coordinates = new double[consideredBones.length][2];
		for (int j = 0; j<consideredBones.length;++j){
			Vector<Integer> sRoiI = tempRoi.edges.get(consideredBones[j]).iit;
			Vector<Integer> sRoiJ =  tempRoi.edges.get(consideredBones[j]).jiit;
			byte[] returnedSieve = tempRoi.fillSieve(sRoiI, sRoiJ, tempRoi.width,tempRoi.height,tempRoi.scaledImage,tempRoi.details.rotationThreshold);
			int counter  = 0;
			coordinates[j][0] = 0;
			coordinates[j][1] = 0;
			for (int h = 0;h<roi.height;++h){
				for (int w = 0;w<roi.width;++w){ 
					if(returnedSieve[w+h*roi.width] ==1){
						coordinates[j][0]+=(double)w;
						coordinates[j][1]+=(double)h;
						counter++;
					}
				}
			}
			coordinates[j][0]/=(double)counter;
			coordinates[j][1]/=(double)counter;
		}
		return coordinates;
	}
	
}
