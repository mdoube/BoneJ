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
import java.awt.image.*; //Creating the image...
import java.awt.*;			//Polygon, Rectangle
import javax.swing.*;   //for createImage
import org.doube.bonej.pqct.io.*;	//image data
import ij.*;		//ImagePlus
import ij.gui.*;	//ImagePlus ROI
import ij.text.*; 	//Debugging ...

@SuppressWarnings("serial")

/*	//For Debugging
	TextWindow checkWindow = new TextWindow(new String("DICOM info"),new String(""),800,400);
	checkWindow.append((String) imp.getProperty("Info"));*/
public class SelectROI extends JPanel{
	public ImageAndAnalysisDetails details;
	public double[] scaledImage;
	public double[] cortexROI;
	public double minimum;
	public double maximum;
	public Vector<Integer> iit;		//indexes for x-coordinates
	public Vector<Integer> jiit;	//indexes for y-coordinates
	public Vector<Integer> roiI;
	public Vector<Integer> roiJ;

	public Vector<Integer> boneMarrowRoiI;
	public Vector<Integer> boneMarrowRoiJ;
	public Vector<Integer> cortexRoiI;	//For BMD analyses
	public Vector<Integer> cortexRoiJ;	//For BMD analyses
	public Vector<Integer> cortexAreaRoiI;	//For AREA analyses
	public Vector<Integer> cortexAreaRoiJ;	//For AREA analyses
	public Vector<Integer> length;
	public Vector<Integer> beginnings;
	public Vector<Integer> lengthMarrow;
	public Vector<Integer> beginningsMarrow;
	public Vector<Double> marrowDensities;	//For storing mean BMD of marrow concentric rings
	public byte leg;
	public int height;
	public int width;

	public double marrowThreshold;
	public double airThreshold;
	public double fatThreshold;
	public double muscleThreshold;
	public double areaThreshold;	//For cortical AREA analyses (CoA, SSI, I) + peeling distal pixels
	public double BMDthreshold;		//For cortical BMD analyses
	public double softThreshold;	//Thresholding soft tissues + marrow from bone
	public double boneThreshold;	//Thresholding bone from the rest and cortical AREA analyses (CoA, SSI, I)
	
	public double pixelSpacing;
	public byte[] result;			//Will contain filled bones
	public byte[] sieve;
	public byte[] marrowSieve;
	public int[] longestEdge;	//For storing which traced edge is the longest (i.e. outlines the bone of interest)
	public int[] marrowLongestEdge;	//For storing which traced edge is the longest (i.e. outlines the bone of interest)

	public String imageSaveName;
	public String imageSavePath;
	ImagePlus imp;
	public int bmcAlfaIndex = 0;
	//ImageJ constructor
	public SelectROI(ScaledImageData dataIn,ImageAndAnalysisDetails detailsIn, ImagePlus imp){

		this.imp = imp;
		details =detailsIn;
		scaledImage = (double[])dataIn.scaledImage.clone();
		pixelSpacing = dataIn.pixelSpacing;
		imageSavePath = details.imageSavePath;
		width =dataIn.width;
		height =dataIn.height;
		

		airThreshold = details.airThreshold;
		fatThreshold = details.fatThreshold;
		muscleThreshold = details.muscleThreshold;
		marrowThreshold = details.marrowThreshold;
		areaThreshold = details.areaThreshold;	//For cortical AREA analyses (CoA, SSI, I) + peeling distal pixels
		BMDthreshold = details.BMDthreshold;		//For cortical BMD analyses
		softThreshold = details.softThreshold;	//Thresholding soft tissues + marrow from bone
		boneThreshold = details.boneThreshold;
		minimum = dataIn.minimum;
		maximum = dataIn.maximum;
		//Select ROI
		
		iit = new Vector<Integer>();
		jiit = new Vector<Integer>();
		length = new Vector<Integer>();
		beginnings = new Vector<Integer>();
		lengthMarrow = new Vector<Integer>();
		beginningsMarrow = new Vector<Integer>();
		marrowDensities = new Vector<Double>();
		longestEdge = new int[4];	//Larger than one to accommodate a second bone if required (e.g. longest = tibia, second longest = fibula)
		marrowLongestEdge = new int[4];	//Larger than one to accommodate a second bone if required (e.g. longest = tibia, second longest = fibula)
		result = new byte[width*height];

		/*Select ROI and set everything else than the roi to minimum*/
		cortexROI = new double[width*height];	//Make a new copy of the image with only the ROI remaining
		marrowSieve= new byte[width*height];
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
		findEdge(tempScaledImage,length,beginnings, iit, jiit,boneThreshold);	//Trace bone edges	
		
		//For Debugging
		TextWindow checkWindow = new TextWindow(new String("Length"),new String(""),400,400);
		checkWindow.append("Monta loytyi "+length.size());
		
		/*Select correct bone outline*/
		int selection = 0;
		if (details.roiChoice.equals(details.choiceLabels[0])){selection = selectRoiBiggestBone(length);}
		if (details.roiChoice.equals(details.choiceLabels[1])){selection = selectRoiSmallestBone(length);}
		if (details.roiChoice.equals(details.choiceLabels[2])){selection = selectRoiLeftMostBone(beginnings,iit);}
		if (details.roiChoice.equals(details.choiceLabels[3])){selection = selectRoiRightMostBone(beginnings,iit);}
		if (details.roiChoice.equals(details.choiceLabels[4])){selection = selectRoiTopMostBone(beginnings,jiit);}
		if (details.roiChoice.equals(details.choiceLabels[5])){selection = selectRoiBottomMostBone(beginnings,jiit);}
		if (details.roiChoice.equals(details.choiceLabels[6])){selection = selectRoiCentralBone(beginnings,length,iit,jiit,tempScaledImage,details.fatThreshold);}
		if (details.roiChoice.equals(details.choiceLabels[7])){selection = selectRoiPeripheralBone(beginnings,length,iit,jiit,tempScaledImage,details.fatThreshold);}
		
		/*
		TextWindow checkWindow = new TextWindow(new String("guessFlip"),new String(""),400,200);
		checkWindow.append((String) );
		*/
		/*Try to guess whether to flip the distribution*/
		if (details.guessFlip && details.stacked){
			details.flipDistribution = guessFlip(beginnings,jiit,selection);
		}
		if (details.guessFlip && !details.stacked){
			details.flipDistribution = guessFlip(beginnings,iit,selection);
		}
		
		/*fill roiI & roiJ*/
		for (int i = beginnings.get(selection);i < beginnings.get(selection)+length.get(selection);i++){
			roiI.add(iit.get(i));
			roiJ.add(jiit.get(i));
		}
		
		/*Cleaving function to separate bones attached with a narrow ridge. Useful e.g. for distal tibia*/
		/*
		if (details.allowCleaving){
			cleaveEdge(roiI,roiJ,3.0,6.0);
		}
		*/
		/*Add the roi to the image*/
		int[] xcoordinates = new int[roiI.size()];
		int[] ycoordinates = new int[roiJ.size()];
		for (int i = 0;i<roiI.size();++i){
			xcoordinates[i] = roiI.get(i);
			ycoordinates[i] = roiJ.get(i);
		}
		ijROI = new PolygonRoi(xcoordinates,ycoordinates,roiI.size(),Roi.POLYGON);
		imp.setRoi(ijROI);
		sieve= new byte[width*height];
		
		/*
		TextWindow checkWindow = new TextWindow(new String("RoiSize"),new String(""),400,200);
		checkWindow.append("roiI size "+roiI.size()+" roiJ "+roiJ.size());
		*/
		
		fillSieve(roiI, roiJ, sieve);

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
	
	boolean guessFlip(Vector<Integer> beginning,Vector<Integer> iit, int selection){
		Vector<Integer> temp = new Vector<Integer>();
		Vector<Integer> temp2 = new Vector<Integer>();
		for (int iii =0;iii<beginning.size();iii++){
			temp.add(iit.get(beginning.get(iii)));
			temp2.add(iit.get(beginning.get(iii)));
		}
		Collections.sort(temp);
		int counter=0;
		while (temp2.get(counter) !=temp.get(0)){
			++counter;
		}
		if (selection == counter){return false;}
		return true;
	}
	
	
	int selectRoiBiggestBone(Vector<Integer> length){
		Vector<Integer> temp = new Vector<Integer>();
		for (int iii =0;iii<length.size();++iii){
			temp.add(length.get(iii));
		}
		Collections.sort(temp);
		int counter=0;
		while (length.get(counter) !=temp.get(temp.size()-1)){
			++counter;
		}
		return counter;
	}
	
	int selectRoiSmallestBone(Vector<Integer> length){
		Vector<Integer> temp = new Vector<Integer>();
		for (int iii =0;iii<length.size();iii++){
			temp.add(length.get(iii));
		}
		Collections.sort(temp);
		int counter=0;
		while (length.get(counter) !=temp.get(0)){
			++counter;
		}
		return counter;
	}
	
	int selectRoiLeftMostBone(Vector<Integer> beginning,Vector<Integer> iit){
		Vector<Integer> temp = new Vector<Integer>();
		Vector<Integer> temp2 = new Vector<Integer>();
		for (int iii =0;iii<beginning.size();iii++){
			temp.add(iit.get(beginning.get(iii)));
			temp2.add(iit.get(beginning.get(iii)));
		}
		Collections.sort(temp);
		int counter=0;
		while (temp2.get(counter) !=temp.get(0)){
			++counter;
		}
		return counter;
	}
	
	int selectRoiRightMostBone(Vector<Integer> beginning,Vector<Integer> iit){
		Vector<Integer> temp = new Vector<Integer>();
		Vector<Integer> temp2 = new Vector<Integer>();
		for (int iii =0;iii<beginning.size();iii++){
			temp.add(iit.get(beginning.get(iii)));
			temp2.add(iit.get(beginning.get(iii)));
		}
		Collections.sort(temp);
		int counter=0;
		while (temp2.get(counter) !=temp.get(temp.size()-1)){
			++counter;
		}
		return counter;
	}
	
	int selectRoiTopMostBone(Vector<Integer> beginning,Vector<Integer> iit){
		Vector<Integer> temp = new Vector<Integer>();
		Vector<Integer> temp2 = new Vector<Integer>();
		for (int iii =0;iii<beginning.size();iii++){
			temp.add(iit.get(beginning.get(iii)));
			temp2.add(iit.get(beginning.get(iii)));
		}
		Collections.sort(temp);
		int counter=0;
		while (temp2.get(counter) !=temp.get(0)){
			++counter;
		}
		return counter;
	}
	
	int selectRoiBottomMostBone(Vector<Integer> beginning,Vector<Integer> iit){
		Vector<Integer> temp = new Vector<Integer>();
		Vector<Integer> temp2 = new Vector<Integer>();
		for (int iii =0;iii<beginning.size();iii++){
			temp.add(iit.get(beginning.get(iii)));
			temp2.add(iit.get(beginning.get(iii)));
		}
		Collections.sort(temp);
		int counter=0;
		while (temp2.get(counter) !=temp.get(temp.size()-1)){
			++counter;
		}
		return counter;
	}
	
	int selectRoiCentralBone(Vector<Integer> beginnings,Vector<Integer> length,Vector<Integer> iit,Vector<Integer> jiit,double[] tempScaledImage,double fatThreshold){
		double[] distanceFromCentreOfLimb = calcDistancesFromCentreOfLimb(beginnings,length,iit,jiit,tempScaledImage,fatThreshold);
		double[] temp = (double[]) distanceFromCentreOfLimb.clone();
		Arrays.sort(temp);
		int counter=0;
		while (distanceFromCentreOfLimb[counter] !=temp[0]){
			++counter;
		}
		return counter;
	}

	int selectRoiPeripheralBone(Vector<Integer> beginnings,Vector<Integer> length,Vector<Integer> iit,Vector<Integer> jiit,double[] tempScaledImage,double fatThreshold){
		double[] distanceFromCentreOfLimb = calcDistancesFromCentreOfLimb(beginnings,length,iit,jiit,tempScaledImage,fatThreshold);
		double[] temp = (double[]) distanceFromCentreOfLimb.clone();
		Arrays.sort(temp);
		int counter=0;
		while (distanceFromCentreOfLimb[counter] !=temp[temp.length-1]){
			++counter;
		}
		return counter;
	}
	
	public double[] calcDistancesFromCentreOfLimb(Vector<Integer> beginnings,Vector<Integer> length,Vector<Integer> iit,Vector<Integer> jiit,double[] tempScaledImage,double fatThreshold){
		double[] softPoints = new double[3];
		for (int j=0;j<3;++j){
				softPoints[j]=0;
			}
		Vector<double[]> bones = new Vector<double[]>();		
		for (int i=0;i<beginnings.size();++i){
			bones.add(new double[3]);
			for (int j=0;j<3;++j){
				bones.get(i)[j]=0;
			}
		}
		/*Find the centre of area of the limb*/
		for (int j = 0; j<height;++j){
			for (int i = 0; i<width;++i){
				if (tempScaledImage[i+j*width]>=fatThreshold){
				softPoints[0]+=i;
				softPoints[1]+=j;
				softPoints[2]+=1;
				}
			}
		}
		softPoints[0]/=softPoints[2];	/*X coordinate of the centre of area of the limb... (assuming just one limb)*/
		softPoints[1]/=softPoints[2];	/*Y coordinate of the centre of area of the limb... (assuming just one limb)*/
		/*Find the centres of circumference of the bones*/
		double[] distanceFromCentreOfLimb = new double[beginnings.size()];
		for (int i=0;i<beginnings.size();++i){
			for (int j = beginnings.get(i);j < beginnings.get(i)+length.get(i);j++){
				bones.get(i)[0]+=iit.get(j);
				bones.get(i)[1]+=jiit.get(j);
				bones.get(i)[2]+=1;
			}
			bones.get(i)[0]/=bones.get(i)[2];
			bones.get(i)[1]/=bones.get(i)[2];
			distanceFromCentreOfLimb[i] =Math.pow(softPoints[0]-bones.get(i)[0],2.0)+Math.pow(softPoints[1]-bones.get(i)[1],2.0); /*Square root omitted, as it does not affect the order...*/
		}
		return distanceFromCentreOfLimb;
	}
		
	/*For density distribution visualization*/
	public BufferedImage getMyImage(double[] imageIn,double[] marrowCenter,Vector<Integer> pind, double[] R, double[] R2, double[] Theta2, 
		int width, int height, double minimum, double maximum, Component imageCreator) {
		int[] image = new int[width*height];
		int pixel;
		for (int x = 0; x < width*height;x++) {
			pixel = (int) (((((double) (imageIn[x] -minimum))/((double)(maximum-minimum)))*255.0)); //Korjaa tama...
			image[x]= 255<<24 | pixel <<16| pixel <<8| pixel; 
		}
		 //Draw rotated radii
		for(int i = 0; i< 360;i++) {//45;i++) {//
			image[((int) (marrowCenter[0]+R[pind.get(i)]*Math.cos(Theta2[i])))+  ((int) (marrowCenter[1]+R[pind.get(i)]*Math.sin(Theta2[i])))*width]= 255<<24 | 255 <<16| 0 <<8| 255;
			image[(int) (marrowCenter[0]+R2[pind.get(i)]*Math.cos(Theta2[i]))+ ((int) (marrowCenter[1]+R2[pind.get(i)]*Math.sin(Theta2[i])))*width]=255<<24 | 0 <<16| 255 <<8| 255;
		}
		 Image imageToDraw = createImage(new MemoryImageSource(width,height,image,0,width));
		 imageToDraw= imageToDraw.getScaledInstance(1000, -1, Image.SCALE_SMOOTH);
		 BufferedImage bufferedImage = (BufferedImage) imageCreator.createImage(imageToDraw.getWidth(null), imageToDraw.getHeight(null));
		 Graphics2D gbuf = bufferedImage.createGraphics();
		 gbuf.drawImage(imageToDraw, 0, 0,null);
		 return bufferedImage;
	}

	/*For cortical analysis only visualization*/
	public BufferedImage getMyImage(double[] imageIn,byte[] sieve,int width, int height, double minimum, double maximum, Component imageCreator) {
		int[] image = new int[width*height];
		int pixel;
		for (int x = 0; x < width*height;x++) {
			pixel = (int) (((((double) (imageIn[x] -minimum))/((double)(maximum-minimum)))*255.0));
			image[x]= 255<<24 | pixel <<16| pixel <<8| pixel; 
			if (sieve[x] == 1){   //Ting roi area color with violet
				image[x]= 255<<24 | pixel <<16| 0 <<8| pixel; 
			}
		}
		Image imageToDraw = createImage(new MemoryImageSource(width,height,image,0,width));
		 imageToDraw= imageToDraw.getScaledInstance(1000, -1, Image.SCALE_SMOOTH);
		 BufferedImage bufferedImage = (BufferedImage) imageCreator.createImage(imageToDraw.getWidth(null), imageToDraw.getHeight(null));
		 Graphics2D gbuf = bufferedImage.createGraphics();
		 gbuf.drawImage(imageToDraw, 0, 0,null);
		 return bufferedImage;
	}
	
	void fillSieve(Vector<Integer> roiI, Vector<Integer> roiJ, byte[] sieveTemp){	
		//Fill the area enclosed by the traced edge contained in roiI,roiJ
		//beginning needs to be within the traced edge
			int z=0;
			int kai,kaj,i,j;
			kai = 0;
			kaj = 0;
			for(z = 0;z<roiI.size();z++){
				kai = kai+roiI.get(z);
				kaj = kaj+roiJ.get(z);
				sieveTemp[roiI.get(z)+roiJ.get(z)*width]=1;
			}
			i = kai/roiI.size();
			j = kaj/roiJ.size();
			
			Vector<Integer> initialI = new Vector<Integer>();
			Vector<Integer> initialJ = new Vector<Integer>();
			initialI.add(i);
			initialJ.add(j);
			while (initialI.size()>0){
				i =initialI.lastElement();
				j =initialJ.lastElement();
				initialI.remove(initialI.size()-1);;
				initialJ.remove(initialJ.size()-1);;
			
				if (sieveTemp[i+j*width] == 0){
					sieveTemp[i+j*width] = 1;
			
				}
				//check whether the neighbour to the left should be added to the que
				if (sieveTemp[i-1+j*width] == 0) {
				initialI.add(i-1);
				initialJ.add(j);
				}
				//check whether the neighbour to the right should be added to the que
				if (sieveTemp[i+1+j*width] == 0) {
				initialI.add(i+1);
				initialJ.add(j);
				}
				//check whether the neighbour below should be added to the que
				if (sieveTemp[i+(j-1)*width] == 0) {
				initialI.add(i);
				initialJ.add(j-1);
				}
				//check whether the neighbour above should be added to the que
				if (sieveTemp[i+(j+1)*width] == 0) {
				initialI.add(i);
				initialJ.add(j+1);
				}
			
			}
	}
	
	/*New algorithm 
			trace edge by advancing according to the previous direction
			if above threshold, turn to negative direction
			if below threshold, turn to positive direction
			Idea taken from some paper, couldn't locate it anymore
			The paper traced continent edges on map/satellite image
		*/
	void traceEdge(double[] scaledImage,byte[] result,double threshold,Vector<Integer> iit,Vector<Integer> jiit,int i,int j){
		double direction = 0; //begin by advancing right. Positive angles rotate the direction clockwise.
		double previousDirection;
		boolean done = false;
		int initI,initJ;
		initI = i;
		initJ = j;
		while(!done){
			int counter = 0;
			previousDirection = direction;
			if (scaledImage[i+((int) Math.round(Math.cos(direction)))+(j+((int) Math.round(Math.sin(direction))))*width] > threshold){//Rotate counter clockwise
				while((scaledImage[i+((int) Math.round(Math.cos(direction-Math.PI/4.0)))+(j+((int) Math.round(Math.sin(direction-Math.PI/4.0))))*width] > threshold 
				)
				&& counter < 8 
				){
					direction-=Math.PI/4.0;
					++counter;
					if (Math.abs(direction-previousDirection) >= 180){
						break;
					}
					
				}
			}else{//Rotate clockwise
				while((scaledImage[i+((int) Math.round(Math.cos(direction)))+(j+((int) Math.round(Math.sin(direction))))*width] < threshold 
				)
				&& counter < 8){
					direction+=Math.PI/4.0;
					++counter;
					if (Math.abs(direction-previousDirection) >= 180){
						break;
					}
				}
			}
			i += (int) Math.round(Math.cos(direction));
			j += (int) Math.round(Math.sin(direction));
			if ((i == initI && j == initJ) || counter > 7 || scaledImage[i+j*width]<threshold || result[i+j*width] ==1 || result[i+j*width] >3){
				done = true;
			}
			else{
				if (result[i+j*width] == 0){
					result[i+j*width] = 2;
				}else if (result[i+j*width] != 1){
					result[i+j*width]++;
				}
				iit.add(i);
				jiit.add(j);
			}
			direction -=Math.PI/2.0; //Keep steering counter clockwise not to miss single pixel structs...
		}	
		for (int ii = 0; ii< result.length;++ii){
			if(result[ii] > 1){result[ii]=1;}
		}
		
	}
	
	boolean resultFill(int i, int j){	
		boolean possible = true;
		Vector<Integer> initialI = new Vector<Integer>();
		Vector<Integer> initialJ= new Vector<Integer>();
		initialI.add(i);
		initialJ.add(j);
		while (initialI.size() >0 && initialI.lastElement() > 0 &&  initialI.lastElement() < width-1 && initialJ.lastElement() > 0 && initialJ.lastElement() < height-1){
			i =initialI.lastElement();
			j = initialJ.lastElement();
			initialI.remove( initialI.size()-1);
			initialJ.remove( initialJ.size()-1);

			if (result[i+j*width] == 0){
				result[i+j*width] = 1;
			}

			if (result[i-1+j*width] == 0) {
			initialI.add(i-1);
			initialJ.add(j);
			}

			if (result[i+1+j*width] == 0) {
			initialI.add(i+1);
			initialJ.add(j);
			}
			
			if (result[i+(j-1)*width] == 0) {
			initialI.add(i);
			initialJ.add(j-1);
			}
			
			if (result[i+(j+1)*width] == 0) {
			initialI.add(i);
			initialJ.add(j+1);
			}

		}

		if (initialI.size() > 0 || initialJ.size()>0) {possible = false;}		
		return possible;
	}
	
	void findEdge(double[] scaledImage,Vector<Integer> length, Vector<Integer> beginnings,Vector<Integer> iit, Vector<Integer> jiit,double threshold)
	{
		int i,j,tempI,tempJ;
		int len;
		i = 0;
		j = 0;

		while ((i < (width-1)) && (j < (height -1) )){
			while (j < height-1 && i < width && scaledImage[i+j*width] <threshold){
				i++;
				if (result[i+j*width] == 1){
					while (j < height-1 && result[i+j*width]>0){
						i++;
						if (i == width && j < height-2){
							i = 0;
							j++;
						}
						
					}
				}

				if (i == width){
					j++;
					if (j >= height-1) break;
					i = 0;
				}
			}
			tempI = i;
			tempJ = j;

			if (i >= width-1 && j >= height-1){
				break;	/*Go to end...*/
			}
			result[i+j*width] = 1;
			Vector<Integer> newIit = new Vector<Integer>();
			Vector<Integer> newJiit = new Vector<Integer>();
			newIit.add(i);
			newJiit.add(j);
			//iit.add(i);
			//jiit.add(j);
			/*Tracing algorithm*/
			traceEdge(scaledImage,result,threshold,newIit,newJiit,i,j);
			len = newIit.size();
			/*Tracing algorithm done...*/
			//Allow cleaving?
			Vector<Vector<Vector<Integer>>>  returnedVectors = null;
			if (details.allowCleaving){
				returnedVectors = cleaveEdge(newIit,newJiit,3.0,6.0);
				//Debugging
				
				TextWindow checkWindow = new TextWindow(new String("ReturnedVectors..."),new String(""),400,200);
				checkWindow.append("size "+returnedVectors.size());
				checkWindow.append("firstSize "+returnedVectors.get(0).get(0).size());
				checkWindow.append("firstIndexI "+returnedVectors.get(0).get(0).get(0)+ "firstIndexJ "+returnedVectors.get(0).get(1).get(0));
				
				for (int iii = 0;iii<returnedVectors.size();++iii){	/*Go through all returned edges*/
					/*Fill edge within result..*/
					for (int ii = 0; ii<returnedVectors.get(iii).get(0).size();++ii){
						iit.add(returnedVectors.get(iii).get(0).get(ii));
						jiit.add(returnedVectors.get(iii).get(1).get(ii));
					}
					len = returnedVectors.get(iii).get(0).size();
					checkWindow.append("length size "+length.size()+" vector length"+returnedVectors.get(iii).get(0).size());
					fillResultEdge(length,beginnings,iit,jiit,len);
					checkWindow.append("length size after "+length.size());
				}
				
			}else{
				/*Fill edge within result..*/
				for (int ii = 0; ii<newIit.size();++ii){
					iit.add(newIit.get(ii));
					jiit.add(newJiit.get(ii));
				}
				fillResultEdge(length,beginnings,iit,jiit,len);
			}
			
			//Find next empty spot
			i = tempI;
			j = tempJ;
			while (j < height && scaledImage[i+j*width] >=threshold){
				i++;
				if (i == width){
				i = 0;
				j++;
				}
			}

		}

	}
	
	void fillResultEdge(Vector<Integer> length, Vector<Integer> beginnings,Vector<Integer> iit, Vector<Integer> jiit,int len){
		if (len > 0){
			length.add(len);
			beginnings.add(iit.size()-len);
			int kai,kaj;
			kai = 0;
			kaj = 0;
			for(int zz = beginnings.lastElement() ;zz<(beginnings.lastElement()+length.lastElement()) ;zz++){
				kai = kai+ iit.get(zz);
				kaj = kaj+ jiit.get(zz);
			}
			kai = kai/length.lastElement() ;
			kaj = kaj/length.lastElement() ;
			while(result[kai+kaj*width]> 1){
				kai = kai+1;
				kaj = kaj+1;
			}
						
			boolean possible = true;
			int jj =kaj;
			int ii;
			for (ii = kai;ii<width;ii++){
				if (result[ii+jj*width]> 0){break;}
			}
			if (ii>=width-1){possible = false;}
			
			for (ii = kai;ii>0;ii--){
				if (result[ii+jj*width]> 0){break;}
			}
			if (ii<=1){possible = false;}
			
			ii = kai;
			for (jj = kaj;jj<height;jj++){
				if (result[ii+jj*width]> 0){break;}
			}
			if (jj>=height-1){possible = false;}
			
			for (jj = kaj;jj>0;jj--){
				if (result[ii+jj*width]> 0){break;}
			}
			if (jj<=1){possible = false;}
			
			if(result[kai+kaj*width]==1){possible = false;}
			
			//Debugging
			TextWindow checkWindow = new TextWindow(new String("Possible"),new String(""),400,200);
			checkWindow.append("Possible "+possible);
			
			if (possible){
				possible = resultFill(kai,kaj);
				
				checkWindow.append("Possible "+possible);
				if (!possible){
					//Remove "extra ii and jii
					for (int po = 0;po <length.lastElement() ;po++){
						iit.remove(iit.size()-1);
						jiit.remove(jiit.size()-1);
					}
					length.remove(length.size()-1);
					beginnings.remove(beginnings.size()-1);
				}
			}else{
				for (int po = 0;po <length.lastElement() ;po++){
					iit.remove(iit.size()-1);
					jiit.remove(jiit.size()-1);
				}
				length.remove(length.size()-1);
				beginnings.remove(beginnings.size()-1);
			}
		}
	}
	
	
	/*Cleaving is made by looking at the ratios of
	distances between two points along the edge and the shortest distance 
	between the points. If the maximum of the  ratio is big enough, the 
	highest ratio points will be connected with a straigth
	line and the edge with higher indices will be removed. E.g. 
	for a circle, the maximum ratio is (pi/2)/d ~= 1.57 and for square
	it is 2/sqrt(2) = sqrt(2) ~= 1.41.*/	
	Vector<Vector<Vector<Integer>>> cleaveEdge(Vector<Integer> fatRoiI,Vector<Integer> fatRoiJ,double minRatio,double minLength){
		double distanceAlongTheEdge = 0;
		double distance = 0;
		double ratio;
		double minEdge = (double) fatRoiI.size()/minLength;
		int[] cleavingIndices = new int[2];
		boolean nextLoop = true;
		Vector<Vector<Vector<Integer>>> returnVectorVectorPointer = new Vector<Vector<Vector<Integer>>>();
		while (nextLoop){
			double highestRatio = minRatio-0.1;
			/*Go through all point pairs*/
			for (int i=0;i<fatRoiI.size()-11;++i){
				for (int j=i+10;j<fatRoiI.size();++j){
					distance = Math.sqrt(Math.pow((double) (fatRoiI.get(j)-fatRoiI.get(i)),2.0)+Math.pow((double) (fatRoiJ.get(j)-fatRoiJ.get(i)),2.0));
					distanceAlongTheEdge = min((double)(j-i),(double) fatRoiI.size()-j+i);
					ratio = distanceAlongTheEdge/distance;
					if (ratio>highestRatio && distanceAlongTheEdge > minEdge){
						highestRatio = ratio;
						cleavingIndices[0] = i;
						cleavingIndices[1] = j;
					}

				}
			}
			/*If ratio is high enough, cleave at the highest ratio point pair*/
			if (highestRatio >= minRatio){
				returnVectorVectorPointer.add(cleave(fatRoiI,fatRoiJ,cleavingIndices));
			} else {
				nextLoop = false;
			}
		}
		/*Insert the last retained part to first index.*/
		Vector<Vector<Integer>> returnVectorPair = new Vector<Vector<Integer>>();
		returnVectorPair.add(new Vector<Integer>());
		returnVectorPair.add(new Vector<Integer>());
		returnVectorPair.get(0).addAll(fatRoiI);
		returnVectorPair.get(1).addAll(fatRoiJ);
		if (returnVectorVectorPointer.size() < 1){
			returnVectorVectorPointer.add(returnVectorPair);
		}else{
			returnVectorVectorPointer.insertElementAt(returnVectorPair, 0);
		}
		return returnVectorVectorPointer;
	}
	
	/*
		Remove the extra part from vectors
		and replace with a straight line
	*/
	Vector<Vector<Integer>> cleave(Vector<Integer> fatRoiI,Vector<Integer> fatRoiJ,int[] cleavingIndices){
		int initialLength = fatRoiI.size();
		int initI = fatRoiI.get(cleavingIndices[0]);
		int initJ = fatRoiJ.get(cleavingIndices[0]);
		int targetI = fatRoiI.get(cleavingIndices[1]);
		int targetJ = fatRoiJ.get(cleavingIndices[1]);
		/*remove cleaved elements*/
		int replacementI = fatRoiI.get(cleavingIndices[0]);
		int replacementJ = fatRoiJ.get(cleavingIndices[0]);
		Vector<Integer> cleavedI = new Vector<Integer>(fatRoiI.subList(cleavingIndices[0]+1,cleavingIndices[1]+1)); /*the elements to be cleaved*/
		Vector<Integer> cleavedJ = new Vector<Integer>(fatRoiJ.subList(cleavingIndices[0]+1,cleavingIndices[1]+1)); /*the elements to be cleaved*/
		for (int i = cleavingIndices[0]; i <cleavingIndices[1];++i){
			fatRoiI.removeElementAt(cleavingIndices[0]);	/*Remove the elements to be cleaved*/
			fatRoiJ.removeElementAt(cleavingIndices[0]);	/*Remove the elements to be cleaved*/
		}
		/*Insert replacement line*/
		double replacementLength = (double)(cleavingIndices[1]-cleavingIndices[0]);
		double repILength = (double)(targetI-initI);
		double repJLength = (double)(targetJ-initJ);
		double relativeLength;
		Vector<Integer> insertionI = new Vector<Integer>();
		Vector<Integer> insertionJ = new Vector<Integer>();
		insertionI.add(replacementI);
		insertionJ.add(replacementJ);
		for (int k = cleavingIndices[0];k<cleavingIndices[1];++k){
			relativeLength = ((double)k)-((double)cleavingIndices[0]);
			replacementI = ((int) (repILength*(relativeLength/replacementLength)))+initI;
			replacementJ = ((int) (repJLength*(relativeLength/replacementLength)))+initJ;
			if (replacementI !=insertionI.lastElement() || replacementJ !=insertionJ.lastElement()){
				insertionI.add(replacementI);
				insertionJ.add(replacementJ);
			}
		}
		fatRoiI.addAll(cleavingIndices[0],insertionI);
		fatRoiJ.addAll(cleavingIndices[0],insertionJ);
		Collections.reverse(insertionI);
		Collections.reverse(insertionJ);
		cleavedI.addAll(0,insertionI);
		cleavedJ.addAll(0,insertionJ);
		Vector<Vector<Integer>> returnVectorPair = new Vector<Vector<Integer>>();
		returnVectorPair.add(new Vector<Integer>());
		returnVectorPair.add(new Vector<Integer>());
		returnVectorPair.get(0).addAll(cleavedI);
		returnVectorPair.get(1).addAll(cleavedJ);
		return  returnVectorPair;
	}
	
	double min(double a,double b){
		return (a < b) ? a : b;
	}
}
