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
import org.doube.jama.*;		//linear equation group solver http://math.nist.gov/javanumerics/jama/
import java.awt.image.*; //Creating the image...
import java.awt.*;			//Polygon, Rectangle
import java.io.*;				//File IO
import javax.imageio.*;		//Saving the image
import javax.swing.*;   //for createImage
import org.doube.bonej.pqct.io.*;	//image data
import ij.*;		//ImagePlus
import ij.gui.*;	//ImagePlus ROI
public class SelectROI extends JPanel{
	ImageAndAnalysisDetails details;
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
	public byte[] result;
	public byte[] sieve;
	public byte[] marrowSieve;
	public byte[] marrowKernelSieve;
	public byte[] legSieve;	//Sieve for calculating fat% contains the whole leg pixels
	public int[] longestEdge;	//For storing which traced edge is the longest (i.e. outlines the bone of interest)
	public int[] marrowLongestEdge;	//For storing which traced edge is the longest (i.e. outlines the bone of interest)

	public String imageSaveName;
	public String imageSavePath;
	ImagePlus imp;
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
		//System.out.println("Edgeen\n");
		//System.out.println("Jalka etsimaan");
		result = new byte[width*height];

				/*Select ROI and set everything else than the roi to minimum*/
		//System.out.println("soft found");
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
		if (ijROI == null){	/*No ROI defined, detect one automatically*/
			/*Make a switch for automated roi detections here...*/
			selectRoiBiggestBone();
		}else{ /*Manually selected ROI*/
			double[] tempScaledImage = (double[]) scaledImage.clone();
			/*Check whether pixel is within ROI, mark with bone threshold*/
			for (int j = 0;j< height;j++){
				for (int i = 0; i < width;i++){
					if (ijROI.contains(i,j)){
					}else{
						tempScaledImage[i+j*width] = minimum;
					}
				}
			}
			findEdge(tempScaledImage,length,beginnings, iit, jiit,longestEdge,boneThreshold);	//Delineate the manual ROI	
		}
		/*fill roiI & roiJ*/
		for (int i = beginnings.get(longestEdge[0]);i < beginnings.get(longestEdge[0])+length.get(longestEdge[0]);i++){
			roiI.add(iit.get(i));
			roiJ.add(jiit.get(i));
		}
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

	void selectRoiBiggestBone(){
		findEdge(scaledImage,length,beginnings, iit, jiit,longestEdge,boneThreshold);	//Bone area analysis
	}
	
	public BufferedImage getMyImage(double[] imageIn,double[] marrowCenter,Vector<Integer> pind, double[] R, double[] R2, double[] Theta2, 
		int width, int height, double minimum, double maximum, Component imageCreator) {
		int[] image = new int[width*height];
		int pixel;
		for (int x = 0; x < width*height;x++) {
			pixel = (int) (((((double) (imageIn[x] -minimum))/((double)(maximum-minimum)))*255.0)); //Korjaa tama...
			image[x]= 255<<24 | pixel <<16| pixel <<8| pixel; 
		}
		 //Draw rotated radii
		for(int i = 0; i< 360;i++) {
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
	void traceEdge(double[] scaledImage,byte[] result,double threshold,Vector<Integer> iit,Vector<Integer> jiit,int[] len,int i,int j){
		double direction = 0; //begin by advancing right. Positive angles rotate the direction clockwise.
		double previousDirection;
		boolean done = false;
		int initI,initJ;
		initI = i;
		initJ = j;
		//System.out.println("I "+i+" J "+j+"InI "+initI+" J "+initJ+" dir "+(direction/Math.PI*180.0)+" D "+ scaledImage[i+j*width]+" r "+result[i+j*width]);
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
						//System.out.println("Negative break");
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
						//System.out.println("Positive break");
						break;
					}
				}
			}
			//if (result[i+j*width] == 1 || initialI.size()<1){ /*Allow returning via already used point*/
			i += (int) Math.round(Math.cos(direction));
			j += (int) Math.round(Math.sin(direction));
			//if ((i == initI && j == initJ) || counter > 7 || i < 1 || j < 1 || i>=width-1|| j>=height-1 || scaledImage[i+j*width]<threshold){
			if ((i == initI && j == initJ) || counter > 7 || scaledImage[i+j*width]<threshold || result[i+j*width] ==1 || result[i+j*width] >3){
				//System.out.println("TEST "+(i == initI && j == initJ)+" "+(counter > 7)+" "+(scaledImage[i+j*width]<threshold));
				//System.out.println("Done "+i+" "+j+" "+counter+" "+scaledImage[i+j*width]);
				done = true;
			}
			else{
				if (result[i+j*width] == 0){
					result[i+j*width] = 2;
				}else if (result[i+j*width] != 1){
					result[i+j*width]++;
				}
				iit.add(new Integer(i));
				jiit.add(new Integer(j));
				len[0]++;
			}
			direction -=Math.PI/2.0; //Keep steering counter clockwise not to miss single pixel structs...
		}
		
		for (int ii = 0; ii< result.length;++ii){
			if(result[ii] > 1){result[ii]=1;}
		}
		
	}
	
	void findEdge_leg(double[] scaledImage,Vector<Integer> length, Vector<Integer> beginnings,Vector<Integer> iit, Vector<Integer> jiit, int[] longestEdge,double threshold)
	{
		int newRow;
		int newCol;
		int i,j,ii,jj,tempI,tempJ;
		
		Vector<Integer> initialI = new Vector<Integer>();
		Vector<Integer> initialJ= new Vector<Integer>(); 
		
		Vector<Integer> previousI = new Vector<Integer>();
		Vector<Integer> previousJ= new Vector<Integer>();
		Vector<Integer> lengths = new Vector<Integer>();
		boolean edgeEnd = false;
		int size;
		size = width*height;
		int[][] searchOrder = {{1,0,-1,-1,-1,0,1,1,1},{1,1,1,0,-1,-1,-1,0,1}};
		
		int[] len = new int[1];
		i = 0;
		j = 0;

		//System.out.println("Edge tracing begins "+i+" "+j);
		while ((i < (width-1)) && (j < (height -1) )){
			len[0] = 0;
			//printf("Prior to i %d j %d\n",i,j);
			previousI.clear();
			previousJ.clear();
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
			//System.out.println("temp I\n");
			tempI = i;
			tempJ = j;
			//Save previous i and j
			previousI.add(new Integer(i));
			previousJ.add(new Integer(j));
			//System.out.println("temp I "+(Integer) previousI.lastElement()+" J"+(Integer) previousJ.lastElement());
			//return if we're at the end of the scaledImage
			
			if (i >= width-1 && j >= height-1){
				break;	/*Go to end...*/
				//return;
			}
			//
			//Edge found
			result[i+j*width] = 1;
			iit.add(new Integer(i));
			jiit.add(new Integer(j));
			len[0]++;
			/*Tracing algorithm*/
			traceEdge(scaledImage,result,threshold,iit,jiit,len,i,j);
			/*Tracing algorithm done...*/

			if (len[0] > 0){
				length.add(new Integer(len[0]));
				lengths.add(new Integer(len[0]));
				if (length.size() < 2){
					beginnings.add(new Integer(0));
				}else{
					beginnings.add(new Integer(iit.size()-len[0]));
				}
				
				
				int kai,kaj;
				kai = 0;
				kaj = 0;
				for(int zz = beginnings.lastElement() ;zz<(beginnings.lastElement()+length.lastElement()) ;zz++){
					/*kai = kai+(Integer) iit.get(zz);
					kaj = kaj+(Integer) jiit.get(zz);*/
					kai = kai+ iit.get(zz);
					kaj = kaj+ jiit.get(zz);
				}
				//printf("Kai %d Kaj %d\n",kai,kaj);
				kai = kai/length.lastElement() ;
				kaj = kaj/length.lastElement() ;
				//System.out.println("kai "+kai+" kaj"+kaj);
				//In case the centre of the traced edge is on top of the trace
				while(result[kai+kaj*width]> 1){
					kai = kai+1;
					kaj = kaj+1;
				}
				
				//Check whether the trace can be filled
			
				boolean possible = true;
				
				jj =kaj;
				for (ii = kai;ii<width;ii++){
					if (result[ii+jj*width]> 0) break;
				}
				
				if (ii>=width-1) possible = false;
				for (ii = kai;ii>0;ii--){
					if (result[ii+jj*width]> 0) break;
				}
				
				if (ii<=1) possible = false;
				i = kai;
				for (jj = kaj;jj<height;jj++){
					if (result[ii+jj*width]> 0) break;
				}
				
				if (jj>=height-1) possible = false;
				for (jj = kaj;jj>0;jj--){
					if (result[ii+jj*width]> 0) break;
				}
				
				if (jj<=1) possible = false;
				
				if(result[kai+kaj*width]==1){
					possible = false;

				}
				//System.out.println(possible);

				if (possible){
					
					possible = resultFill(kai,kaj);
					
					if (!possible){
						//Remove "extra ii and jii
						for (int po = 0;po <length.lastElement() ;po++){
							iit.remove(iit.size()-1);
							jiit.remove(jiit.size()-1);
						}
						length.remove(length.size()-1);
						beginnings.remove(beginnings.size()-1);
					// return;
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
			//System.out.println("length "+length.size()+" beg "+beginnings.size());
		}
		whichLeg(lengths);
		
		longestEdge[0] = 0;
		longestEdge[1] = 0;
		Vector<Integer> temp = new Vector<Integer>();
		for (int iii =0;iii<length.size();iii++){
			temp.add(length.get(iii));
		}
		Collections.sort(temp);
		int counter=0;
		while (length.get(counter) !=temp.get(temp.size()-1)){
			++counter;
		}
		longestEdge[0]=counter;
		if (length.size() > 1){
			counter=0;
			while (length.get(counter) !=temp.get(temp.size()-2)){
				++counter;
			}
			longestEdge[1]=counter;
		}
		
	}
	
	void whichLeg(Vector<Integer> lengths)
	{
		leg = 0;
		int[] pitit = new int[2];
		
		

			//Find out which is tibia -> Look for the largest...
		int jar = 0;
		for (int i = 0; i<2;i++){
			int maxima = 0;
			for (int is = 0; is < lengths.size();is++){
				//if ((Integer) lengths.get(is) > (Integer) lengths.get(maxima)){maxima = is;}
				if (lengths.get(is) > lengths.get(maxima)){maxima = is;}
			}
			pitit[jar] = maxima;
			lengths.setElementAt(0,maxima);
			++jar;
		}
		if (pitit[0] < pitit[1]){
			leg = 0;  //false =right
		}else{
			leg = 1;  //true = left
		}	
	}
	boolean resultFill(int i, int j){	
		boolean possible = true;
		Vector<Integer> initialI = new Vector<Integer>();
		Vector<Integer> initialJ= new Vector<Integer>();
		initialI.add(new Integer(i));
		initialJ.add(new Integer(j));

		////System.out.println("I "+i+" J "+j);
		while (initialI.size() >0 && initialI.lastElement() > 0 &&  initialI.lastElement() < width-1 && initialJ.lastElement() > 0 && initialJ.lastElement() < height-1){
			//System.out.println("I "+i+" J "+j);
			i =initialI.lastElement();
			j = initialJ.lastElement();
			initialI.remove( initialI.size()-1);
			initialJ.remove( initialJ.size()-1);
			//System.out.println("2I "+i+" J "+j);
			if (result[i+j*width] == 0){
				result[i+j*width] = 1;
			}

			if (result[i-1+j*width] == 0) {
			initialI.add(new Integer(i-1));
			initialJ.add(new Integer(j));
			}

			if (result[i+1+j*width] == 0) {
			initialI.add(new Integer(i+1));
			initialJ.add(new Integer(j));
			}
			
			if (result[i+(j-1)*width] == 0) {
			initialI.add(new Integer(i));
			initialJ.add(new Integer(j-1));
			}
			
			if (result[i+(j+1)*width] == 0) {
			initialI.add(new Integer(i));
			initialJ.add(new Integer(j+1));
			}

		}

		if (initialI.size() > 0 || initialJ.size()>0) {possible = false;}
		
		return possible;
	}
	
	void findEdge(double[] scaledImage,Vector<Integer> length, Vector<Integer> beginnings,Vector<Integer> iit, Vector<Integer> jiit, int[] longestEdge,double threshold)
	{
		int newRow;
		int newCol;
		int i,j,ii,jj,tempI,tempJ;
		
		Vector<Integer> initialI = new Vector<Integer>();	
		Vector<Integer> initialJ= new Vector<Integer>(); 
		

		Vector<Integer> previousI = new Vector<Integer>();
		Vector<Integer> previousJ= new Vector<Integer>();
		Vector<Integer> lengths = new Vector<Integer>();
		boolean edgeEnd = false;
		int size;
		size = width*height;
		int[][] searchOrder = {{1,0,-1,-1,-1,0,1,1,1},{1,1,1,0,-1,-1,-1,0,1}};
		
		int[] len = new int[1];
		i = 0;
		j = 0;

		while ((i < (width-1)) && (j < (height -1) )){
			len[0] = 0;
			previousI.clear();
			previousJ.clear();
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

			previousI.add(new Integer(i));
			previousJ.add(new Integer(j));
			if (i >= width-1 && j >= height-1){
				break;	/*Go to end...*/
			}
			result[i+j*width] = 1;
			iit.add(new Integer(i));
			jiit.add(new Integer(j));
			len[0]++;

			/*Tracing algorithm*/
			traceEdge(scaledImage,result,threshold,iit,jiit,len,i,j);
			/*Tracing algorithm done...*/
			
			if (len[0] > 0){
				length.add(new Integer(len[0]));
				lengths.add(new Integer(len[0]));
				if (length.size() < 2){
					beginnings.add(new Integer(0));
				}else{
					beginnings.add(new Integer(iit.size()-len[0]));
				}
				
				
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
					jj =kaj;
				for (ii = kai;ii<width;ii++){
					if (result[ii+jj*width]> 0) break;
				}
				
				if (ii>=width-1) possible = false;
				for (ii = kai;ii>0;ii--){
					if (result[ii+jj*width]> 0) break;
				}
				
				if (ii<=1) possible = false;
				i = kai;
				for (jj = kaj;jj<height;jj++){
					if (result[ii+jj*width]> 0) break;
				}
				
				if (jj>=height-1) possible = false;
				for (jj = kaj;jj>0;jj--){
					if (result[ii+jj*width]> 0) break;
				}
				
				if (jj<=1) possible = false;
				
				if(result[kai+kaj*width]==1){
					possible = false;

				}
				

				if (possible){
					possible = resultFill(kai,kaj);
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

		
		longestEdge[0] = 0;
		longestEdge[1] = 0;
		Vector<Integer> temp = new Vector<Integer>();
		temp.addAll(length);
		Collections.sort(temp);
		int counter=0;
		while (length.get(counter) !=temp.get(temp.size()-1)){
			++counter;
		}
		longestEdge[0]=counter;
		if (length.size() > 1){
			counter=0;
			while (length.get(counter) !=temp.get(temp.size()-2)){
				++counter;
			}
			longestEdge[1]=counter;
			/*In case of dicom image having both the left and the right leg, switch left leg tibia to longestEdge[0]*/
			if (iit.get(beginnings.get(longestEdge[0])) < iit.get(beginnings.get(longestEdge[1]))){ //the one with smaller column is left...

			} else {
				int tempo = longestEdge[0];
				longestEdge[0] = longestEdge[1];
				longestEdge[1] = tempo;
			}
		}

		
		/*Search for Figula in DICOM having both feet*/
		if (length.size() > 3){
			counter=0;
			while (length.get(counter) !=temp.get(temp.size()-3)){
				++counter;
			}
			longestEdge[2]=counter; /*First fibula*/
			counter=0;
			while (length.get(counter) !=temp.get(temp.size()-4)){
				++counter;
			}
			longestEdge[3]=counter; /*Second fibula*/
			/*In case of dicom image having both the left and the right leg, switch left leg tibia to longestEdge[0]*/
			if (iit.get(beginnings.get(longestEdge[2])) < iit.get(beginnings.get(longestEdge[3]))){ //the one with smaller column is left...
				longestEdge[1] = longestEdge[2];
			} else {
				longestEdge[1] = longestEdge[3];
			}
		}

	}
	
	/*
	Calculate Spline 
	linear equation group solver http://math.nist.gov/javanumerics/jama/
	cubic spline http://www.physics.arizona.edu/~restrepo/475A/Notes/sourcea/node35.html
	*/

	public static Vector<Integer> calculateSpline(Vector<Integer> interpolants, int numberOfPointsBetweenInterpolants){
		double[][] a = new double[interpolants.size()][interpolants.size()];
		double[] b = new double[interpolants.size()];
		double [] x = new double[interpolants.size()];
		double[] bb = new double[interpolants.size()-1];
		double [] d = new double[interpolants.size()-1];
		for (int i = 0;i<interpolants.size();i++){	//Create xs
			x[i] = (double) i;
		}
		a[0][0] = 1;
		a[interpolants.size()-1][interpolants.size()-1] = 1;
		b[0] = 0;
		b[interpolants.size()-1] = 0;
		for (int i = 1; i< interpolants.size()-1;i++){
			a[i][i-1] =x[i]-x[i-1];
			a[i][i]   =2.0*((x[i]-x[i-1])+(x[i+1]-x[i]));
			a[i][i+1] =x[i+1]-x[i];
			b[i] = (3/(x[i+1]-x[i]))*((double) interpolants.get(i+1)-(double) interpolants.get(i))-(3/(x[i]-x[i-1]))*((double) interpolants.get(i)-(double) interpolants.get(i-1));
		}
		Matrix A = new Matrix(a);
		Matrix B = new Matrix(b,interpolants.size());
		Matrix c = A.solve(B);
		for (int i = 1; i < interpolants.size();i++){
			bb[i-1] = 1/(x[i]-x[i-1])*((double) interpolants.get(i)-(double) interpolants.get(i-1))-(x[i]-x[i-1])/3*(2*c.get(i-1,0)+c.get(i,0));
			d[i-1] = (c.get(i,0)-c.get(i-1,0))/(3*(x[i]-x[i-1]));
		}
		Vector<Integer> result = new Vector<Integer>();
		double inter;
		for (int i = 0;i<interpolants.size()-1;i++){
		   for (int j = 0;j<numberOfPointsBetweenInterpolants;j++){
			   inter = (x[i+1]-x[i])*(double)j/(double)numberOfPointsBetweenInterpolants;
			   result.add((int) ((double) interpolants.get(i)+bb[i]*inter+c.get(i,0)*Math.pow(inter,2.0)+d[i]*Math.pow(inter,3.0)));
		   }
		}
		return result;
	}

}
