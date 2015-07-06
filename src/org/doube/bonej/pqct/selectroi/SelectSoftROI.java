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

@SuppressWarnings(value ={"serial","unchecked"}) //Unchecked for obtaining Vector<Object> as a returnvalue

public class SelectSoftROI extends RoiSelector{
	//ImageJ constructor
	public SelectSoftROI(ScaledImageData dataIn,ImageAndAnalysisDetails detailsIn, ImagePlus imp,double boneThreshold,boolean setRoi){
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
			double[] tempScaledImage = (double[]) softScaledImage.clone();
			if (ijROI != null && details.manualRoi){	/*Set pixels outside the manually selected ROI to zero*/
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




			Vector<Object> masks = getSieve(softScaledImage,airThreshold,details.roiChoiceSt,details.guessStacked,details.stacked,false,true);
			softSieve						= (byte[]) masks.get(0);
			softResult					 	= (byte[]) masks.get(1);
			Vector<DetectedEdge> stEdges	= (Vector<DetectedEdge>) masks.get(2);

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
			Vector<Object> muscleMasks = getSieve(muscleImage,details.muscleThreshold,"Bigger",details.guessStacked,details.stacked,false,false);
			//muscleSieve		= (byte[]) muscleMasks.get(0);
			Vector<DetectedEdge> muscleEdges = (Vector<DetectedEdge>) muscleMasks.get(2);
			Collections.sort(muscleEdges,Collections.reverseOrder());
			int tempMuscleArea=0;
			muscleSieve = new byte[softSieve.length];
			int areaToAdd=0;
			/*Include areas that contribute more than 1% on top of what is already included*/
			while (areaToAdd< muscleEdges.size() && tempMuscleArea*0.01 < muscleEdges.get(areaToAdd).area){
				byte[] tempMuscleSieve = fillSieve(muscleEdges.get(areaToAdd).iit, muscleEdges.get(areaToAdd).jiit,width,height,muscleImage,details.muscleThreshold);
				for (int i = 0; i<tempMuscleSieve.length;++i){
					if (tempMuscleSieve[i] > 0){muscleSieve[i] = tempMuscleSieve[i];}
				}
				tempMuscleArea+=muscleEdges.get(areaToAdd).area;
				areaToAdd++;
			}

			//muscleSieve		= (byte[]) muscleMasks.get(1);	/*Use all areas encircled as muscle (needed if there's fat between muscles)!!*/

			/*Wipe muscle area +1 layer of pixels away from subcut.*/
			byte[] tempMuscleSieve = (byte[]) muscleSieve.clone();
			dilate(tempMuscleSieve,(byte)1,(byte)0,(byte)2);
			for (int i = 0;i<tempMuscleSieve.length;++i){
				if (tempMuscleSieve[i] == 1){subCutaneousFat[i] = 0;}
			}
			
			/*Pop-up muscleSieve image*/
			ImagePlus ipVis = new ImagePlus("Visual");
			ipVis.setProcessor(new ByteProcessor(width,height,muscleSieve));
			new ImageConverter(ipVis).convertToRGB();
			
			//ipVis.show();
			
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
			//IJ.log("Into r");
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
			
						//Visualize segmentation
			ipVis.getProcessor().setColor(new Color(255,0,0));
			
			for (int i = 0;i<359;++i){
				ipVis.getProcessor().drawLine((int) (Math.round(rs[i]*Math.cos(theta[i])+softCentre[0])), (int) (Math.round(rs[i]*Math.sin(theta[i])+softCentre[1])), (int) (Math.round(rs[i+1]*Math.cos(theta[i+1])+softCentre[0])), (int) (Math.round(rs[i+1]*Math.sin(theta[i+1])+softCentre[1])));
			}
			//IMPLEMENT REPEATING LIVEWIRE AFTER ROTATION!!! Get seed points from livewire result...
			//Use every tenth as the init for livewire, rotate twice (5 deg each)
			double[][] pixels = new double[width][height];
			for (int rr = 0;rr<height;++rr){
				for (int c = 0;c<width;++c){
					pixels[c][rr] = (double) muscleSieve[c+rr*width];
				}
			}
			
			ArrayList<Integer> edgeii = new ArrayList<Integer>();
			ArrayList<Integer> edgejj = new ArrayList<Integer>();
			LiveWireCosts lwc = new LiveWireCosts(pixels);
			
			for (int i = 0;i<360-10; i=i+10){
				lwc.setSeed(edgeCoords[i][0],edgeCoords[i][1]);
				while (lwc.returnPath(edgeCoords[i+10][0],edgeCoords[i+10][1]) == null){
					try{Thread.sleep(10);}catch(Exception e){}
				}
				int[][] fromSeedToCursor = lwc.returnPath(edgeCoords[i+10][0],edgeCoords[i+10][1]);
				for (int iii = 0;iii< fromSeedToCursor.length;++iii){
					edgeii.add((int) fromSeedToCursor[iii][0]);
					edgejj.add((int) fromSeedToCursor[iii][1]);
					IJ.log("Edge Length "+edgeii.size()+" x "+edgeii.get(edgeii.size()-1)+" y "+edgejj.get(edgejj.size()-1));
				}
				
			}
			/*Connect the last bit*/
			lwc.setSeed(edgeii.get(edgeii.size()-1),edgejj.get(edgejj.size()-1));
			while (lwc.returnPath(edgeCoords[0][0],edgeCoords[0][1]) == null){
				try{Thread.sleep(10);}catch(Exception e){}
			}
			int[][] fromSeedToCursor = lwc.returnPath(edgeCoords[0][0],edgeCoords[0][1]);
			for (int i = 0;i< fromSeedToCursor.length;++i){
				edgeii.add((int) fromSeedToCursor[i][0]);
				edgejj.add((int) fromSeedToCursor[i][1]);
			}
			//Visualize liveWire result...
			ipVis.getProcessor().setColor(new Color(0,255,0));
			for (int i = 0;i<edgeii.size()-1;++i){
				ipVis.getProcessor().drawLine(edgeii.get(i),edgejj.get(i),edgeii.get(i+1),edgejj.get(i+1));
			}
			//IJ.log("Got past r");

			ipVis.setDisplayRange(0,1);
			ipVis.show();
			ipVis.repaintWindow();
			/*Re-segmenting done*/
			
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
		}
	}
	
	double max(double a,double b){
		return a >= b ? a:b;
	}
}
