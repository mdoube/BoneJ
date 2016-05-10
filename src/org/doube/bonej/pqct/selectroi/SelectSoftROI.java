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
import org.doube.bonej.pqct.selectroi.liveWireEngine.*;	//LiveWire
import java.util.concurrent.ExecutionException;

@SuppressWarnings(value ={"serial","unchecked"}) //Unchecked for obtaining Vector<Object> as a returnvalue

public class SelectSoftROI extends RoiSelector{
	//ImageJ constructor
	public SelectSoftROI(ScaledImageData dataIn,ImageAndAnalysisDetails detailsIn, ImagePlus imp,double boneThreshold,boolean setRoi) throws ExecutionException{
		super(dataIn,detailsIn, imp,boneThreshold,setRoi);
		//Soft tissue analysis
		softSieve = null;
		byte[] softResult = null;
		if (details.stOn){

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


			/**Ignore data outside manually selected ROI, if manualRoi has been selected*/
			Roi ijROI = imp.getRoi();
			if (ijROI != null && details.manualRoi){	/*Set pixels outside the manually selected ROI to zero*/
			double[] tempScaledImage = Arrays.copyOf(softScaledImage,softScaledImage.length);
				/*Check whether pixel is within ROI, mark with bone threshold*/
				for (int j = 0;j< height;j++){
					for (int i = 0; i < width;i++){
						if (ijROI.contains(i,j)){
						}else{
							softScaledImage[i+j*width] = minimum;
						}
					}
				}
				/*Check whether a polygon can be acquired and include polygon points too*/
				Polygon polygon = ijROI.getPolygon();
				if (polygon != null){
					for (int j = 0;j< polygon.npoints;j++){
						softScaledImage[polygon.xpoints[j]+polygon.ypoints[j]*width] = tempScaledImage[polygon.xpoints[j]+polygon.ypoints[j]*width];
					}
				}
			}



			try{
				Vector<Object> masks = getSieve(softScaledImage,airThreshold,details.roiChoiceSt,details.guessStacked,details.stacked,false,true);
				softSieve						= (byte[]) masks.get(0);
				softResult					 	= (byte[]) masks.get(1);
				Vector<DetectedEdge> stEdges	= (Vector<DetectedEdge>) masks.get(2);

				/*Erode three layers of pixels from the fat sieve to get rid of higher density layer (i.e. skin)
				on top of fat to enable finding muscle border
				*/
				byte[] muscleSieve = Arrays.copyOf(softSieve,softSieve.length);
				double[] muscleImage = Arrays.copyOf(softScaledImage,softScaledImage.length);
				byte[] subCutaneousFat = null;
				//Remove skin by eroding three layers of pixels
				for (int i = 0;i< 3;++i){
					muscleSieve = erode(muscleSieve);
					//Changed 2016/01/08
					/*
					if (i == 0){
						//Add subcut fat sieve... Skin has already been removed by eroding one layer of pixels-> remove muscle later on
						subCutaneousFat = (byte[]) muscleSieve.clone();
					}
					*/
				}
				subCutaneousFat = Arrays.copyOf(muscleSieve,muscleSieve.length);
				/*Remove everything other than the selected limb from the image*/
				for (int i = 0; i<muscleSieve.length;++i){
					if (muscleSieve[i] < 1){
						muscleImage[i] = minimum;
					}
				}
				/*Look for muscle outline*/
				Vector<Object> muscleMasks = getSieve(muscleImage,details.muscleThreshold,"Bigger",details.guessStacked,details.stacked,false,false);
				//muscleSieve		= (byte[]) muscleMasks.get(0);
				Vector<DetectedEdge> muscleEdges = (Vector<DetectedEdge>) muscleMasks.get(2);
				Collections.sort(muscleEdges,Collections.reverseOrder());
				int tempMuscleArea=0;
				muscleSieve = new byte[softSieve.length];
				int areaToAdd=0;
				/*Include areas that contribute more than 1.0% on top of what is already included*/
				while (areaToAdd< muscleEdges.size() && tempMuscleArea*0.01 < muscleEdges.get(areaToAdd).area){
					byte[] tempMuscleSieve = fillSieve(muscleEdges.get(areaToAdd).iit, muscleEdges.get(areaToAdd).jiit,width,height,muscleImage,details.muscleThreshold);
					for (int i = 0; i<tempMuscleSieve.length;++i){
						if (tempMuscleSieve[i] > 0){muscleSieve[i] = tempMuscleSieve[i];}
					}
					tempMuscleArea+=muscleEdges.get(areaToAdd).area;
					areaToAdd++;
				}

				//Visualise muscleSieve
				
				ImagePlus tempImage = NewImage.createByteImage("MuscleSieve",width,height,1, NewImage.FILL_BLACK);
				byte[] rPixels = (byte[])tempImage.getProcessor().getPixels();
				for (int i = 0;i<muscleSieve.length;++i){
					for (int c = 0;c<width;++c){
						rPixels[i] = muscleSieve[i];
					}
				}
				tempImage.setDisplayRange(0,1);
				tempImage.show();
				
				/**Re-segment soft-tissues using livewire based on the muscleSieve
					1) bring rays back from image edges to centre of soft-tissue mask 1 deg apart
					2) use livewire on the 360 edge pixels
					3) rotate livewire init pixels around a few times to get the segmentation to go through subcut/intramuscular fat
				
				*/
				ArrayList<Double> ii = new ArrayList<Double>();
				ArrayList<Double> jj = new ArrayList<Double>();
				double[] softCentre = new double[2];
				for (int i = 0;i<width;++i){
					for (int j = 0; j<height;++j){
						if (muscleSieve[i+j*width] > 0){
							ii.add((double) i);
							jj.add((double) j);						
							softCentre[0]+=(double) i;
							softCentre[1]+=(double) j;
						}
					}
				}
				softCentre[0]/=(double)ii.size();
				softCentre[1]/=(double)jj.size();
				double maxR = Math.sqrt(Math.pow(max(Math.abs(softCentre[0]-width),softCentre[0]),2d)+Math.pow(max(Math.abs(softCentre[1]-height),softCentre[1]),2d));
				double[] rs = new double[360];
				double r,t;
				double[] theta = new double[360];
				int[][] edgeCoords = new int[360][2];
				//Get the extremes of muscle area with 1 deg increments in polar coordinates
				for (int i = 0;i<360;++i){
					r = maxR;
					t = ((double)i)/180d*Math.PI;
					while ( Math.round(r*Math.cos(t)+softCentre[0]) < 0 || Math.round(r*Math.sin(t)+softCentre[1]) < 0 ||
							Math.round(r*Math.cos(t)+softCentre[0]) >= width || Math.round(r*Math.sin(t)+softCentre[1]) >= height){
						r-=0.1;
					}
					while ( Math.round(r*Math.cos(t)+softCentre[0]) >= 0 && Math.round(r*Math.sin(t)+softCentre[1]) >= 0 &&
							Math.round(r*Math.cos(t)+softCentre[0]) < width && Math.round(r*Math.sin(t)+softCentre[1]) < height &&
							muscleSieve[(int) (Math.round(r*Math.cos(t)+softCentre[0])+Math.round(r*Math.sin(t)+softCentre[1])*width)] < 1){
						r-=0.1;
					}
					rs[i] = r;
					theta[i] = t;
					edgeCoords[i][0]=(int) (Math.round(rs[i]*Math.cos(theta[i])+softCentre[0]));
					edgeCoords[i][1]=(int) (Math.round(rs[i]*Math.sin(theta[i])+softCentre[1]));
				}

				double[][] pixels = new double[width][height];
				for (int rr = 0;rr<height;++rr){
					for (int c = 0;c<width;++c){
						pixels[c][rr] = (double) muscleSieve[c+rr*width];
					}
				}
				//Arraylists for edge, and livewire seed coordinates
				ArrayList<Integer> edgeii = new ArrayList<Integer>();
				ArrayList<Integer> edgejj = new ArrayList<Integer>();
				ArrayList<Integer> seedii = new ArrayList<Integer>();
				ArrayList<Integer> seedjj = new ArrayList<Integer>();
				//Create list of seed coordinates
				double tempR,tempR2;
				int maxInd;
				for (int i = 0;i<360;i=i+10){
					//Look for the furthest point in this bracket
					maxInd = i;
					tempR = Math.sqrt(Math.pow(edgeCoords[i][0]-softCentre[0],2d)+Math.pow(edgeCoords[i][1]-softCentre[1],2d));
					for (int j = i+1;j<i+10;++j){
						tempR2 = Math.sqrt(Math.pow(edgeCoords[j][0]-softCentre[0],2d)+Math.pow(edgeCoords[j][1]-softCentre[1],2d));
						if (tempR2 > tempR){
							maxInd =j;
							tempR = tempR2;
						}
					}
					seedii.add(edgeCoords[maxInd][0]);
					seedjj.add(edgeCoords[maxInd][1]);
				}
				
				//Loop livewire 6 times over here
				for (int l = 0; l<6;++l){//details.edgeDivisions;++l){
					edgeii.clear();
					edgejj.clear();
					LiveWireCosts lwc = new LiveWireCosts(pixels);
					int[][] fromSeedToCursor;
					for (int i = 0;i<seedii.size()-1; ++i){
						lwc.setSeed(seedii.get(i),seedjj.get(i));
						while ((fromSeedToCursor = lwc.returnPath(seedii.get(i+1),seedjj.get(i+1))) == null){
							try{Thread.sleep(1);}catch(Exception e){}
						}
						for (int iii = 0;iii< fromSeedToCursor.length;++iii){
							edgeii.add((int) fromSeedToCursor[iii][0]);
							edgejj.add((int) fromSeedToCursor[iii][1]);
							//IJ.log("Edge Length "+edgeii.size()+" x "+edgeii.get(edgeii.size()-1)+" y "+edgejj.get(edgejj.size()-1));
						}
						
					}
					/*Connect the last bit*/
					lwc.setSeed(seedii.get(seedii.size()-1),seedjj.get(seedjj.size()-1));
					
					while ((fromSeedToCursor = lwc.returnPath(seedii.get(0),seedjj.get(0))) == null){
						try{Thread.sleep(1);}catch(Exception e){}
					}
					for (int i = 0;i< fromSeedToCursor.length;++i){
						edgeii.add((int) fromSeedToCursor[i][0]);
						edgejj.add((int) fromSeedToCursor[i][1]);
					}
					//Set new seeds
					seedii.clear();
					seedjj.clear();
					double divisions = details.edgeDivisions;
					
					//Get appropriate seeds, select the furthest point from the centre
					for (int i = (int) (l*((edgeii.size()/divisions)/6.));i<(int) ((divisions-1d)/divisions*edgeii.size());i+=1d/divisions*edgeii.size()){
						//Look for the furthest point in this bracket
						maxInd = i;
						tempR = Math.sqrt(Math.pow(edgeii.get(i)-softCentre[0],2d)+Math.pow(edgejj.get(i)-softCentre[1],2d));
						for (int j = i+1;j<i+1d/divisions*edgeii.size();++j){
							tempR2 = Math.sqrt(Math.pow(edgeii.get(j)-softCentre[0],2d)+Math.pow(edgejj.get(j)-softCentre[1],2d));
							if (tempR2 > tempR){
								maxInd =j;
								tempR = tempR2;
							}
						}
						seedii.add(edgeii.get(maxInd));
						seedjj.add(edgejj.get(maxInd));
					}
				}
				
				//Fill in muscle mask with inter-muscular fat 
				muscleSieve = getByteMask(width,height,edgeii,edgejj);
				muscleSieve = dilateMuscleMask(muscleSieve,softScaledImage,width,height,muscleThreshold); //Dilate the sieve to include all muscle pixels
				/*Re-segmenting done*/
				
				//Wipe muscle area +3 layer of pixels away from subcut.
				byte[] tempMuscleSieve = Arrays.copyOf(muscleSieve,muscleSieve.length);
				//dilate(tempMuscleSieve,(byte)1,(byte)0,(byte)2);
				//dilate(tempMuscleSieve,(byte)1,(byte)0,(byte)2);
				//dilate(tempMuscleSieve,(byte)1,(byte)0,(byte)2);
				for (int i = 0;i<tempMuscleSieve.length;++i){
					if (tempMuscleSieve[i] == 1){subCutaneousFat[i] = 0;}
				}
				
				/*create temp boneResult to wipe out bone and marrow*/
				Vector<Object> masks2 = getSieve(softScaledImage,softThreshold,details.roiChoiceSt,details.guessStacked,details.stacked,false,false);
				byte[] boneResult	= (byte[]) masks2.get(1);
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
			}catch (ExecutionException err){
				throw err;
			}
		}
	}
	
	//Helper function to get seed points for livewire
	public Vector<Object> getSeedPoints(ArrayList<Integer> edgeii, ArrayList<Integer> edgejj, double[] softCentre,double divisions, double steps, double l){
		ArrayList<Integer> seedii = new ArrayList<Integer>();
		ArrayList<Integer> seedjj = new ArrayList<Integer>();
		double tempR,tempR2;
		int maxInd;
		//Get appropriate seeds, select the furthest point from the centre
		for (int i = (int) (l*((edgeii.size()/divisions)/steps));i<(int) ((divisions-1d)/divisions*edgeii.size());i+=1d/divisions*edgeii.size()){
			//Look for the furthest point in this bracket
			maxInd = i;
			tempR = Math.sqrt(Math.pow(edgeii.get(i)-softCentre[0],2d)+Math.pow(edgejj.get(i)-softCentre[1],2d));
			for (int j = i+1;j<i+1d/divisions*edgeii.size();++j){
				tempR2 = Math.sqrt(Math.pow(edgeii.get(j)-softCentre[0],2d)+Math.pow(edgejj.get(j)-softCentre[1],2d));
				if (tempR2 > tempR){
					maxInd =j;
					tempR = tempR2;
				}
			}
			seedii.add(edgeii.get(maxInd));
			seedjj.add(edgejj.get(maxInd));
		}
		Vector<Object> returnVal = new Vector<Object>();
		returnVal.add(seedii);
		returnVal.add(seedjj);
		return returnVal;		
	}
	
	double max(double a,double b){
		return a >= b ? a:b;
	}
	
	
	public byte[] dilateMuscleMask(byte[] mask,double[] softScaledImage,int width, int height,double threshold){
		ArrayList<Integer> initialI = new ArrayList<Integer>();
		ArrayList<Integer> initialJ= new ArrayList<Integer>();
		int i,j;
		for (i = 0;i<width;++i){
			for (j = 0;j<height;++j){
				if (mask[i+j*width]==1){
					initialI.add(i);
					initialJ.add(j);
				}
			}
		}
		
		while (initialI.size() >0 && initialI.get(initialI.size()-1) > 0 &&  initialI.get(initialI.size()-1) < width-1
			&& initialJ.get(initialJ.size()-1) > 0 && initialJ.get(initialJ.size()-1) < height-1){
			i =initialI.get(initialI.size()-1);
			j = initialJ.get(initialJ.size()-1);
			initialI.remove( initialI.size()-1);
			initialJ.remove( initialJ.size()-1);

			if (mask[i+j*width] == 0 && softScaledImage[i+j*width] >=threshold){
				mask[i+j*width] = 1;
			}

			if (mask[i-1+j*width] == 0 && softScaledImage[i+j*width] >=threshold) {
				initialI.add(i-1);
				initialJ.add(j);
			}

			if (mask[i+1+j*width] == 0 && softScaledImage[i+j*width] >=threshold) {
				initialI.add(i+1);
				initialJ.add(j);
			}
			
			if (mask[i+(j-1)*width] == 0 && softScaledImage[i+j*width] >=threshold) {
				initialI.add(i);
				initialJ.add(j-1);
			}
			
			if (mask[i+(j+1)*width] == 0 && softScaledImage[i+j*width] >=threshold) {
				initialI.add(i);
				initialJ.add(j+1);
			}

		}
		return mask;
	}
	
	byte[] getByteMask(int width,int height,ArrayList<Integer> edgeii,ArrayList<Integer> edgejj){
		byte[] mask=new byte[width*height];
		for (int i = 0; i<edgeii.size();++i){
			mask[edgeii.get(i)+edgejj.get(i)*width] = (byte) 1;
		}
		int[] fillInitCoords = findMaskFillInit(mask,width,height,edgeii,edgejj);
		if (fillInitCoords != null){
			return fillMask(fillInitCoords[0],fillInitCoords[1],mask,width,height);
		}else{
			return mask;
		}
	}
	
	byte[] fillMask(int i, int j, byte[] mask, int width,int height){
		ArrayList<Integer> initialI = new ArrayList<Integer>();
		ArrayList<Integer> initialJ= new ArrayList<Integer>();
		initialI.add(i);
		initialJ.add(j);
		while (initialI.size() >0 && initialI.get(initialI.size()-1) > 0 &&  initialI.get(initialI.size()-1) < width-1
			&& initialJ.get(initialJ.size()-1) > 0 && initialJ.get(initialJ.size()-1) < height-1){
			i =initialI.get(initialI.size()-1);
			j = initialJ.get(initialJ.size()-1);
			initialI.remove( initialI.size()-1);
			initialJ.remove( initialJ.size()-1);

			if (mask[i+j*width] == 0 ){
				mask[i+j*width] = 1;
			}

			if (mask[i-1+j*width] == 0) {
				initialI.add(i-1);
				initialJ.add(j);
			}

			if (mask[i+1+j*width] == 0) {
				initialI.add(i+1);
				initialJ.add(j);
			}
			
			if (mask[i+(j-1)*width] == 0) {
				initialI.add(i);
				initialJ.add(j-1);
			}
			
			if (mask[i+(j+1)*width] == 0) {
				initialI.add(i);
				initialJ.add(j+1);
			}

		}
		return mask;
	}
	
	int[] findMaskFillInit(byte[] mask,int width,int height,ArrayList<Integer> edgeii,ArrayList<Integer> edgejj){
		int[] returnCoordinates = new int[2];
		int[] steer = new int[2];
		for (int j = 0; j< edgeii.size()-1; ++j){
			returnCoordinates[0] = edgeii.get(j);
			returnCoordinates[1] = edgejj.get(j);
			double direction = Math.atan2(edgejj.get(j+1)-returnCoordinates[1],edgeii.get(j+1)-returnCoordinates[0]);
			direction+=Math.PI/4.0;
			for (int i = 0; i< 3; ++i){
				
				steer[0] = (int) Math.round(Math.cos(direction));
				steer[1]= (int) Math.round(Math.sin(direction));
				/*Handle OOB*/
				while ((returnCoordinates[0]+steer[0])<0 || (returnCoordinates[0]+steer[0])>=width ||
						(returnCoordinates[1]+steer[1])<0 || (returnCoordinates[1]+steer[1])>=height){
					direction+=Math.PI/4.0;
					steer[0] = (int) Math.round(Math.cos(direction));
					steer[1]= (int) Math.round(Math.sin(direction));
				}
				if (mask[returnCoordinates[0]+steer[0]+(returnCoordinates[1]+steer[1])*width] == 0){
					returnCoordinates[0] +=steer[0];
					returnCoordinates[1] +=steer[1];
					return returnCoordinates;
				}
				if (result[returnCoordinates[0]+steer[0]+(returnCoordinates[1]+steer[1])*width] == 1){
					break;
				}
				direction+=Math.PI/4.0;
			}
		}
		return null;
	}
	
}
