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

package SelectRoi;
import java.util.*;	//Vector, Collections
import java.lang.Math; //atan2
import Jama.*;		//linear equation group solver http://math.nist.gov/javanumerics/jama/
import java.awt.image.*; //Creating the image...
import java.awt.*;
import java.awt.Graphics2D;
import java.io.*;				//File IO
import javax.imageio.*;		//Saving the image
import javax.swing.*;   //for createImage
import ImageReading.*;	//image data
import DrawImage.*;		//Drawing and saving images

public class SelectROI extends JPanel implements Runnable{
	ImageAndAnalysisDetails details;
	public double[] scaledImage;
	public double[] softScaledImage;
	public double[] cortexROI;
	public double minimum;
	public double maximum;
	public Vector<Integer> iit;		//indexes for x-coordinates
	public Vector<Integer> jiit;	//indexes for y-coordinates
	public Vector<Integer> roiI;
	public Vector<Integer> roiJ;
	public Vector<Integer> softIit;
	public Vector<Integer> softJiit;
	public Vector<Integer> softRoiI;
	public Vector<Integer> softRoiJ;
	public Vector<Integer> marrowIit;
	public Vector<Integer> marrowJiit;
	public Vector<Integer> marrowRoiI;
	public Vector<Integer> marrowRoiJ;
	public Vector<Integer> medullaryMarrowRoiI;
	public Vector<Integer> medullaryMarrowRoiJ;
	public Vector<Integer> boneMarrowRoiI;
	public Vector<Integer> boneMarrowRoiJ;
	public Vector<Integer> cortexRoiI;	//For BMD analyses
	public Vector<Integer> cortexRoiJ;	//For BMD analyses
	public Vector<Integer> cortexAreaRoiI;	//For AREA analyses
	public Vector<Integer> cortexAreaRoiJ;	//For AREA analyses
	public Vector<Integer> length;
	public Vector<Integer> beginnings;
	public Vector<Integer> lengthSoft;
	public Vector<Integer> beginningsSoft;
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
	public byte[] softSieve;
	public byte[] marrowSieve;
	public byte[] marrowKernelSieve;
	public byte[] muscleSieve;
	public byte[] legSieve;	//Sieve for calculating fat% contains the whole leg pixels
	public byte[] stSieve;		//Soft Tissue Sieve, What will be 1 on this one is the muscle = fat. Bones will be removed
	public int[] longestEdge;	//For storing which traced edge is the longest (i.e. outlines the bone of interest)
	public int[] softLongestEdge;	//For storing which traced edge is the longest (i.e. outlines the bone of interest)
	public int[] marrowLongestEdge;	//For storing which traced edge is the longest (i.e. outlines the bone of interest)

	DrawImage showFigure;
	public String imageSaveName;
	public String imageSavePath;

	public boolean imageJ;	//used to indicate to not try to print out images while using the object from ImageJ
	//ImageJ constructor
	public SelectROI(ScaledImageData dataIn,ImageAndAnalysisDetails detailsIn){
		imageJ = true; //Used to indicate that the object has been constructed with imageJ
		details =detailsIn;
		scaledImage = (double[])dataIn.scaledImage.clone();
		softScaledImage = (double[])dataIn.softScaledImage.clone();
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
		softIit = new Vector<Integer>();
		softJiit = new Vector<Integer>();
		marrowIit = new Vector<Integer>();
		marrowJiit = new Vector<Integer>();
		length = new Vector<Integer>();
		beginnings = new Vector<Integer>();
		lengthSoft = new Vector<Integer>();
		beginningsSoft = new Vector<Integer>();
		lengthMarrow = new Vector<Integer>();
		beginningsMarrow = new Vector<Integer>();
		marrowDensities = new Vector<Double>();
		longestEdge = new int[4];	//Larger than one to accommodate a second bone if required (e.g. longest = tibia, second longest = fibula)
		softLongestEdge = new int[4];	//Larger than one to accommodate a second bone if required (e.g. longest = tibia, second longest = fibula)
		marrowLongestEdge = new int[4];	//Larger than one to accommodate a second bone if required (e.g. longest = tibia, second longest = fibula)
		//System.out.println("Edgeen\n");
		//System.out.println("Jalka etsimaan");
		result = new byte[width*height];
		if (details.dicomOn == true){
			findEdge(scaledImage,length,beginnings, iit, jiit,longestEdge,boneThreshold);	//Bone area analysis
			leg = 1;  //true = vasen
			//System.out.println("Bone found "+marrowThreshold);
			
			if (details.mOn == true){
				result = new byte[width*height];
				findEdge(scaledImage,lengthMarrow,beginningsMarrow, marrowIit,marrowJiit,marrowLongestEdge,marrowThreshold);	//marrow analysis
			}
			
			//System.out.println("marrow found");
			if (details.stOn ==true){
				result = new byte[width*height];
				findEdge(softScaledImage,lengthSoft,beginningsSoft, softIit,softJiit,softLongestEdge,softThreshold);	//Soft tissue analysis
			}
		}else{
			findEdge_leg(scaledImage,length,beginnings, iit, jiit,longestEdge,boneThreshold);	//Bone area analysis
			//System.out.println("Bone found "+marrowThreshold);
			
			if (details.mOn == true){
				result = new byte[width*height];
				findEdge_leg(scaledImage,lengthMarrow,beginningsMarrow, marrowIit,marrowJiit,marrowLongestEdge,marrowThreshold);	//marrow analysis
			}
			
			//System.out.println("marrow found");
			if (details.stOn ==true){
				result = new byte[width*height];
				findEdge_leg(softScaledImage,lengthSoft,beginningsSoft, softIit,softJiit,softLongestEdge,softThreshold);	//Soft tissue analysis
			}
		}
		/*Select ROI and set everything else than the roi to minimum*/
		//System.out.println("soft found");
		cortexROI = new double[width*height];	//Make a new copy of the image with only the ROI remaining
		sieve= new byte[width*height];
		softSieve= new byte[width*height];
		marrowSieve= new byte[width*height];
		stSieve= new byte[width*height];
		legSieve = new byte[width*height]; //Includes the whole cross-section of the leg
		muscleSieve= new byte[width*height];
		roiI = new Vector<Integer>();
		roiJ = new Vector<Integer>();
		softRoiI = new Vector<Integer>();
		softRoiJ = new Vector<Integer>();
		marrowRoiI = new Vector<Integer>();
		marrowRoiJ = new Vector<Integer>();
		cortexRoiI = new Vector<Integer>();
		cortexRoiJ = new Vector<Integer>();
		cortexAreaRoiI = new Vector<Integer>();
		cortexAreaRoiJ = new Vector<Integer>();
		boneMarrowRoiI = new Vector<Integer>();
		boneMarrowRoiJ = new Vector<Integer>();
		medullaryMarrowRoiI = new Vector<Integer>();
		medullaryMarrowRoiJ = new Vector<Integer>();
		//System.out.println("Luominen valmis ");
	}
	
	
	public SelectROI(ScaledImageData dataIn,ImageAndAnalysisDetails detailsIn, DrawImage showFigureIn,String imageSaveNameIn){
		
		details =detailsIn;
		scaledImage = (double[])dataIn.scaledImage.clone();
		softScaledImage = (double[])dataIn.softScaledImage.clone();
		showFigure = showFigureIn;
		pixelSpacing = dataIn.pixelSpacing;
		imageSavePath = details.imageSavePath;
		imageSaveName = imageSaveNameIn;
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
		softIit = new Vector<Integer>();
		softJiit = new Vector<Integer>();
		marrowIit = new Vector<Integer>();
		marrowJiit = new Vector<Integer>();
		length = new Vector<Integer>();
		beginnings = new Vector<Integer>();
		lengthSoft = new Vector<Integer>();
		beginningsSoft = new Vector<Integer>();
		lengthMarrow = new Vector<Integer>();
		beginningsMarrow = new Vector<Integer>();
		marrowDensities = new Vector<Double>();
		longestEdge = new int[4];	//Larger than one to accommodate a second bone if required (e.g. longest = tibia, second longest = fibula)
		softLongestEdge = new int[4];	//Larger than one to accommodate a second bone if required (e.g. longest = tibia, second longest = fibula)
		marrowLongestEdge = new int[4];	//Larger than one to accommodate a second bone if required (e.g. longest = tibia, second longest = fibula)
		//System.out.println("Edgeen\n");
		//System.out.println("Jalka etsimaan");
		result = new byte[width*height];
		if (details.dicomOn == true){
			findEdge(scaledImage,length,beginnings, iit, jiit,longestEdge,boneThreshold);	//Bone area analysis
			leg = 1;  //true = vasen
			//System.out.println("Bone found "+marrowThreshold);
			/*
			if (details.mOn == true){
				result = new byte[width*height];
				findEdge(scaledImage,lengthMarrow,beginningsMarrow, marrowIit,marrowJiit,marrowLongestEdge,marrowThreshold);	//marrow analysis
			}
			*/
			//System.out.println("marrow found");
			if (details.stOn ==true){
				result = new byte[width*height];
				findEdge(softScaledImage,lengthSoft,beginningsSoft, softIit,softJiit,softLongestEdge,softThreshold);	//Soft tissue analysis
			}
		}else{
			findEdge_leg(scaledImage,length,beginnings, iit, jiit,longestEdge,boneThreshold);	//Bone area analysis
			//System.out.println("Bone found "+marrowThreshold);
			/*
			if (details.mOn == true){
				result = new byte[width*height];
				findEdge_leg(scaledImage,lengthMarrow,beginningsMarrow, marrowIit,marrowJiit,marrowLongestEdge,marrowThreshold);	//marrow analysis
			}
			*/
			//System.out.println("marrow found");
			if (details.stOn ==true){
				result = new byte[width*height];
				findEdge_leg(softScaledImage,lengthSoft,beginningsSoft, softIit,softJiit,softLongestEdge,softThreshold);	//Soft tissue analysis
			}
		}
		/*Select ROI and set everything else than the roi to minimum*/
		//System.out.println("soft found");
		cortexROI = new double[width*height];	//Make a new copy of the image with only the ROI remaining
		sieve= new byte[width*height];
		softSieve= new byte[width*height];
		marrowSieve= new byte[width*height];
		stSieve= new byte[width*height];
		legSieve = new byte[width*height]; //Includes the whole cross-section of the leg
		muscleSieve= new byte[width*height];
		roiI = new Vector<Integer>();
		roiJ = new Vector<Integer>();
		softRoiI = new Vector<Integer>();
		softRoiJ = new Vector<Integer>();
		marrowRoiI = new Vector<Integer>();
		marrowRoiJ = new Vector<Integer>();
		cortexRoiI = new Vector<Integer>();
		cortexRoiJ = new Vector<Integer>();
		cortexAreaRoiI = new Vector<Integer>();
		cortexAreaRoiJ = new Vector<Integer>();
		boneMarrowRoiI = new Vector<Integer>();
		boneMarrowRoiJ = new Vector<Integer>();
		medullaryMarrowRoiI = new Vector<Integer>();
		medullaryMarrowRoiJ = new Vector<Integer>();
		//System.out.println("Luominen valmis ");
	}
	public void run(){
		//System.out.println("beg "+beginnings.get(longestEdge[0])+" le "+length.get(longestEdge[0]));
		for (int i = beginnings.get(longestEdge[0]);i < beginnings.get(longestEdge[0])+length.get(longestEdge[0]);i++){
			roiI.add(iit.get(i));
			roiJ.add(jiit.get(i));
		}
		if (details.stOn == true){
			for (int i = beginningsSoft.get(softLongestEdge[0]);i < beginningsSoft.get(softLongestEdge[0])+lengthSoft.get(softLongestEdge[0]);i++){
				softRoiI.add(softIit.get(i));
				softRoiJ.add(softJiit.get(i));
			}
		}
		/*
		if (details.mOn == true){
			for (int i = beginningsMarrow.get(marrowLongestEdge[0]);i < beginningsMarrow.get(marrowLongestEdge[0])+lengthMarrow.get(marrowLongestEdge[0]);i++){
				marrowRoiI.add(marrowIit.get(i));
				marrowRoiJ.add(marrowJiit.get(i));
			}
		}
		*/
		double[] tempImage = (double[])scaledImage.clone();
		for (int i = 0;i<roiI.size();i++){
			tempImage[roiI.get(i)+roiJ.get(i)*width] = maximum;
		}
		
		if (imageJ ==false){
			showFigure.drawImage(tempImage,width, height,minimum,maximum);
			showFigure.paintImmediately(0,0,500,500);
		}
		if (details.mRoiDet ==true){	//Delineate ROI manually
			System.out.println("Digitize coordinates with left mouse button, right click once on the figure when you're done.");
			System.out.println("If ROI is OK, just right click once on the figure.");
			showFigure.notReady = true;
			showFigure.coordx.clear();
			showFigure.coordy.clear();
			while (showFigure.notReady == true){
				try{Thread.sleep(100); }catch (Exception err){System.err.println("Didn't finnish digitizing: "+err.getMessage());}
			}
			if(showFigure.coordx.size() > 3){
				showFigure.coordx.add(showFigure.coordx.firstElement());	//Close the ROI
				showFigure.coordy.add(showFigure.coordy.firstElement());	//Close the ROI
				roiI = calculateSpline(showFigure.coordx, 1000);
				roiJ = calculateSpline(showFigure.coordy, 1000);
				tempImage = (double[])scaledImage.clone();
            
				for (int i = 0;i<roiI.size();i++){
					tempImage[roiI.get(i)+roiJ.get(i)*width] = maximum;
				}
				showFigure.drawImageSpline(tempImage,width, height,minimum,maximum,roiI,roiJ);
				showFigure.paintImmediately(0,0,500,500);
			}
		}

		
		

		fillSieve(roiI, roiJ, sieve);
		/*
		if (details.mOn){
			fillSieve(marrowRoiI, marrowRoiJ, marrowSieve);
		}
		*/
		if (details.stOn){
			fillSieve(softRoiI, softRoiJ, softSieve);
		}
		for (int j = 0;j< height;j++){
			for (int i = 0; i < width;i++){
				/*
				if (details.mOn){
					if (scaledImage[i+j*width]<marrowThreshold && marrowSieve[i+j*width] > 0){
						medullaryMarrowRoiI.add(i);
						medullaryMarrowRoiJ.add(j);
					}
				}
				*/
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
		marrowSieve = new byte[width*height];	//marrowSieve needs to be cleared here...
		//ADD bone marrow ROI, erode marrow area by two layers of pixels
		if (details.mOn){
			/*
			for (int i=0;i<medullaryMarrowRoiI.size();i++){
				marrowSieve[medullaryMarrowRoiI.get(i)+width*medullaryMarrowRoiJ.get(i)]=1;
			}
			*/
			for (int i=0;i<boneMarrowRoiI.size();i++){
				marrowSieve[boneMarrowRoiI.get(i)+width*boneMarrowRoiJ.get(i)]=1;
			}
			
			//System.out.println("MarrowSize "+medullaryMarrowRoiI.size());
			//How many layers are to be eroded (try 0, 1, 2,3,4)
			
			/* //No eroding for now...
			for (int ee =0;ee<2;ee++){
				erode(marrowSieve);
			}
			*/
			int eroded =0;
			int marCSA=0;
			for (int i=0;i<height*width;i++){
				if (marrowSieve[i] == 1){marCSA++;}
			}
			//Marrow densities in concentric layers..
			eroded =0;
			byte [] marrowConcentric = (byte[])marrowSieve.clone();	//marrowSieve needs to be cleared here...
			while (eroded < marCSA){//while (eroded < boneMarrowRoiI.size()){
				eroded+=erodeConcentric(marrowConcentric, scaledImage, marrowDensities);
			}
			System.out.println("MarrowRings "+marrowDensities.size());
			//Marrow densities in concentric layers done
			//Add dividing the marrow area to half
			marrowKernelSieve = (byte[])marrowSieve.clone();	//marrowSieve needs to be cleared here...
			eroded =0;
			while (eroded < marCSA/2){
				eroded+=erodeKernel(marrowKernelSieve);
			}
			//Marrow area divided...
			//Bone marrow analysis done
		}
		if(details.stOn==true){
			//Add soft Tissue ROI
			for (int i=0;i<height*width;i++){
				stSieve[i] = 0;	
				legSieve[i] = 0;	
				if (softScaledImage[i] > airThreshold){// && softScaledImage[i] < muscleThreshold){
					stSieve[i] = 1;	
					legSieve[i] = 1;
				}
			}
			//System.out.println("stSieveTehty");
			softTissueRemoveBones(lengthSoft,beginningsSoft, softIit,softJiit,softLongestEdge);
			
			//Just for UKK institute data, get rid of the sleeve!!!
			if (details.ukkOn){
				UKK_special(softScaledImage,20.0);
			}
			//UKK special ends
			for (int i=0;i<height*width;i++){
				muscleSieve[i] = 0;	
				if (stSieve[i] ==1 && softScaledImage[i] >= fatThreshold && softScaledImage[i] < muscleThreshold){
					muscleSieve[i] = 1;	
				}
				//Remove muscle from stSieve
				if (stSieve[i] ==1 && softScaledImage[i] < airThreshold || softScaledImage[i] >= fatThreshold){
					stSieve[i] = 0;	
				}
			}
			

			//Soft Tissue ROI done..
		}
		//DRAW FIGURES
		if (details.imOn ==true && details.dOn ==false){
			try{
				 //Save ROI selection image
				 BufferedImage bi = getMyImage(scaledImage,roiI,roiJ,marrowSieve,marrowKernelSieve,stSieve,muscleSieve,legSieve); // retrieve image
				// System.out.println("Saatiin bi tehtya");
					File imageOutputFile = new File(details.imageSavePath+"/"+ imageSaveName+".png");
				// System.out.println("Saatiin oputFile tehtya");
					ImageIO.write(bi, "png", imageOutputFile);
				// System.out.println("Saatiin tulostettua");
			}catch (Exception err){System.err.println("Couldn't write image out: "+err.getMessage());}
		}
		if (details.imOn ==true && details.stOn == true){
			try{
				 //Save ROI selection image
				 BufferedImage bi = getMyImage(softScaledImage,roiI,roiJ,marrowSieve,marrowKernelSieve,stSieve,muscleSieve,legSieve); // retrieve image
				// System.out.println("Saatiin bi tehtya");
					File imageOutputFile = new File(details.imageSavePath+"/"+ imageSaveName+"soft.png");
				// System.out.println("Saatiin oputFile tehtya");
					ImageIO.write(bi, "png", imageOutputFile);
				 //System.out.println("Saatiin tulostettua");
			}catch (Exception err){System.err.println("Couldn't write image out: "+err.getMessage());}
		}
		//Remove extra pixels from the external surface of the bone
		if (details.dBoneSite == true && details.eOn == true){
			showFigure.drawImage(cortexROI,width, height,minimum,maximum);
			showFigure.paintImmediately(0,0,500,500);
			try{Thread.sleep(500); }catch (Exception err){System.err.println("Didn't finnish digitizing: "+err.getMessage());}
			while (peelDistal() > 0){System.out.println("Eroded one round");}
			showFigure.drawImage(cortexROI,width, height,minimum,maximum);
			showFigure.paintImmediately(0,0,500,500);
			try{Thread.sleep(500); }catch (Exception err){System.err.println("Didn't finnish digitizing: "+err.getMessage());}
		}
		
	}
	
	int erodeConcentric(byte[] data, double[] densities, Vector<Double> marrowDensities){
		int removed=0;
		double value =0;
		//Erode algorithm
		//Modified from the best dilate by one solution taken from http://ostermiller.org/dilate_and_erode.html
		for (int i=0; i<height; i++){
			for (int j=0; j<width; j++){
				if (data[i*width+j] > 0){
					if (i>0 && data[(i-1)*width+j]==0 ||
						j>0 && data[(i)*width+j-1]==0 ||
						i+1<height && data[(i+1)*width+j]==0 ||
						j+1<width && data[(i)*width+j+1]==0)
					{//Erode the pixel if any of the neighborhood pixels is background
						data[i*width+j] = 0-1;
						++removed;
						value+=densities[i*width+j];
					}	
						
				}
			}
		}
		for (int i=0; i<width*height; i++){
			if (data[i] < 0){
				data[i] = 0;
			}
		}
		value/=(double)removed;
		marrowDensities.add(value);
		return removed;
	}
	
	int erodeKernel(byte[] data){
		int removed=0;
		//Erode algorithm
		//Modified from the best dilate by one solution taken from http://ostermiller.org/dilate_and_erode.html
		for (int i=0; i<height; i++){
			for (int j=0; j<width; j++){
				if (data[i*width+j] > 0){
					if (i>0 && data[(i-1)*width+j]==0 ||
						j>0 && data[(i)*width+j-1]==0 ||
						i+1<height && data[(i+1)*width+j]==0 ||
						j+1<width && data[(i)*width+j+1]==0)
					{//Erode the pixel if any of the neighborhood pixels is background
						data[i*width+j] = 0-1;
						++removed;
					}	
						
				}
			}
		}
		for (int i=0; i<width*height; i++){
			if (data[i] < 0){
				data[i] = 0;
			}
		}
		return removed;
	}
	
	
	void erode(byte[] data){
		//Erode algorithm
		//Modified from the best dilate by one solution taken from http://ostermiller.org/dilate_and_erode.html
		for (int i=0; i<height; i++){
			for (int j=0; j<width; j++){
				if (data[i*width+j] > 0){
					if (i>0 && data[(i-1)*width+j]==0 ||
						j>0 && data[(i)*width+j-1]==0 ||
						i+1<height && data[(i+1)*width+j]==0 ||
						j+1<width && data[(i)*width+j+1]==0)
						{data[i*width+j] = 0-1;}	//Erode the pixel if any of the neighborhood pixels is background
				}
			}
		}
		for (int i=0; i<width*height; i++){
			if (data[i] < 0){
				data[i] = 0;
			}
		}
	}

	void softTissueRemoveBones(Vector<Integer> length, Vector<Integer> beginnings,Vector<Integer> iit, Vector<Integer> jiit, int[] longestEdge){	
		//Fill the area enclosed by the traced edge contained in roiI,roiJ
		//beginning needs to be within the traced edge
		int kai,kaj,i,j;
		//First Bone
		//Set traced edge to 2 and determine the seed point for the flood fill
		kai = 0;
		kaj = 0;
		for (int ii = beginnings.get(longestEdge[0]);ii < beginnings.get(longestEdge[0])+length.get(longestEdge[0]);ii++){
			kai+=iit.get(ii);
			kaj+=jiit.get(ii);
			stSieve[iit.get(ii)+width*jiit.get(ii)]=2;
		}
		i = kai/length.get(longestEdge[0]);
		j = kaj/length.get(longestEdge[0]);
		softTissueRemoveBonesFill(i,j);
		//Second Bone
		//System.out.println("Beginning "+beginnings.size());
		if (beginnings.size() > 1){
			kai = 0;
			kaj = 0;
			for (int ii = beginnings.get(longestEdge[1]);ii < beginnings.get(longestEdge[1])+length.get(longestEdge[1]);ii++){
				kai+=iit.get(ii);
				kaj+=jiit.get(ii);
				stSieve[iit.get(ii)+width*jiit.get(ii)]=2;
			}
			//i = kai/length.get(longestEdge[1]);
			//j = kaj/length.get(longestEdge[1]);
			i = iit.get(beginnings.get(longestEdge[1]));
			j = jiit.get(beginnings.get(longestEdge[1]));
			int[][] searchOrder = {{1,1,1,0,-1,-1,-1,0,1},{1,0,-1,-1,-1,0,1,1,1}};
			int z = 0;
			while (z < 9 && (softScaledImage[i+searchOrder[0][z]+(j+searchOrder[1][z])*width] <softThreshold) || stSieve[i+searchOrder[0][z]+(j+searchOrder[1][z])*width] == 2){
				++z;
			}
			
			softTissueRemoveBonesFill(i+searchOrder[0][z],j+searchOrder[1][z]);	
		}
		
		dilate(stSieve,((byte)2),((byte)1),((byte)3));	//To Remove remaining bone after applying a threshold of 690!
		dilate(stSieve,((byte)2),((byte)1),((byte)3));	//To Remove remaining bone after applying a threshold of 690!
			//printf("sihti i %d j %d\n",i,j);

	}
	
	public void dilate(byte[] data,byte dilateVal,byte min, byte temp){
	//Dilate algorithm
	// Best dilate by one solution taken from http://ostermiller.org/dilate_and_erode.html
	
		for (int i=0; i<height; i++){
			for (int j=0; j<width; j++){
				if (data[i*width+j] ==dilateVal){
					if (i>0 && data[(i-1)*width+j]==min) {data[(i-1)*width+j] = temp;}
					if (j>0 && data[(i)*width+j-1]==min) {data[(i)*width+j-1] = temp;}
					if (i+1<height && data[(i+1)*width+j]==min) {data[(i+1)*width+j] = temp;}
					if (j+1<width && data[(i)*width+j+1]==min) {data[(i)*width+j+1] = temp;}
				}
			}
		}
		for (int i=0; i<width*height; i++){
			if (data[i] == temp){
				data[i] = dilateVal;	//Set to proper value here...
			}
		}
	}
	
	
	void softTissueRemoveBonesFill(int i, int j){
		Vector<Integer> initialI = new Vector<Integer>();
		Vector<Integer> initialJ = new Vector<Integer>();
		initialI.add(i);
		initialJ.add(j);
		int montakoPoistettiin=0;
		while (initialI.size()>0){
			i =initialI.lastElement();
			j =initialJ.lastElement();
			initialI.remove(initialI.size()-1);;
			initialJ.remove(initialJ.size()-1);;
			//Do we paint the cell
			if (stSieve[i+j*width] == 1){
				stSieve[i+j*width] = 2;
				++montakoPoistettiin;
			}
			
			if (stSieve[i-1+j*width] == 1) {
				initialI.add(i-1);
				initialJ.add(j);
			}
			
			if (stSieve[i+1+j*width] == 1) {
				initialI.add(i+1);
				initialJ.add(j);
			}
			
			if (stSieve[i+(j-1)*width] == 1) {
				initialI.add(i);
				initialJ.add(j-1);
			}
			
			if (stSieve[i+(j+1)*width] == 1) {
				initialI.add(i);
				initialJ.add(j+1);
			}
		}
		
	}
	
	
	int peelDistal(){
		int pixelsPeeled = 0;
		//Use Erode algorithm for peeling extra pixels off...
		//Modified from the best dilate by one solution taken from http://ostermiller.org/dilate_and_erode.html
		for (int i=0; i<height; i++){
			for (int j=0; j<width; j++){
				if (sieve[i*width+j] > 0){
					if ((i>0 && sieve[(i-1)*width+j]==0 && cortexROI[(i)*width+j] < areaThreshold) ||
						(j>0 && sieve[(i)*width+j-1]==0 && cortexROI[(i)*width+j] < areaThreshold) ||
						(i+1<height && sieve[(i+1)*width+j]==0 && cortexROI[(i)*width+j] < areaThreshold) ||
						(j+1<width && sieve[(i)*width+j+1]==0 )&& cortexROI[(i)*width+j] < areaThreshold)
						{
							cortexROI[i*width+j] = minimum-1;
							sieve [i*width+j] = -1;
							pixelsPeeled++;
						}	//Erode the pixel if any of the neighborhood pixels is background
				}
			}
		}
		for (int i=0; i<width*height; i++){
			if (cortexROI[i] < minimum){
				cortexROI[i] = minimum;
				sieve[i] = 0;
			}
		}
		return pixelsPeeled;
	}

     //  private static final int DEFAULT_IMAGE_TYPE = BufferedImage.TYPE_INT_RGB;

		public BufferedImage getMyImage(double[] imageIn,Vector<Integer> xx,Vector<Integer> yy, byte[] marrowSieve,byte[] marrowKernelSieve, byte[] stSieve, byte[] muscleSieve, byte[] legSieve) {
			 int[] image = new int[width*height];
			 int pixel;
			 for (int x = 0; x < width*height;x++) {
				pixel = (int) (((((double) (imageIn[x] -minimum))/((double)(maximum-minimum)))*255.0));
				image[x]= 255<<24 | pixel <<16| pixel <<8| pixel; 

				
				if (marrowSieve[x] == 1){   //Change marrow area color to green
					image[x]= 255<<24 | 0 <<16| pixel+50<<8| 0; 
					if (marrowKernelSieve[x] == 1){   //Change marrow area color to green
						image[x]= 255<<24 | pixel+50 <<16| pixel+50<<8| 0; 
					}
				}
				
				if (stSieve[x] == 1){   //Change soft Tissue area color to red
					image[x]= 255<<24 | 0 <<16| 0 <<8| pixel; 
				}
				if (muscleSieve[x] == 1){   //Change muscle area color to violet
					image[x]= 255<<24 | pixel <<16| 0 <<8| 0; 
				}
				
			 }

			 Image imageToDraw = createImage(new MemoryImageSource(width,height,image,0,width));
			 imageToDraw= imageToDraw.getScaledInstance(1000, -1, Image.SCALE_SMOOTH);
			 //System.out.println("Scaled the image to draw");
			 BufferedImage bufferedImage = (BufferedImage) showFigure.createImage(imageToDraw.getWidth(null), imageToDraw.getHeight(null));
			 //System.out.println("BI done "+bufferedImage);
			 Graphics2D gbuf = bufferedImage.createGraphics();
			 //System.out.println("Graphics done");
			 gbuf.drawImage(imageToDraw, 0, 0,null);
			 //System.out.println("Drawn");
			 return bufferedImage;
		}

	
		public BufferedImage getMyImage(double[] imageIn,double[] marrowCenter,Vector<Integer> pind, double[] R, double[] R2, double[] Theta2, 
				byte[] marrowSieve,byte[] marrowKernelSieve, byte[] stSieve, byte[] muscleSieve, byte[] legSieve,
				int width, int height, double minimum, double maximum) {
			 int[] image = new int[width*height];
			 int pixel;
			 for (int x = 0; x < width*height;x++) {
				pixel = (int) (((((double) (imageIn[x] -minimum))/((double)(maximum-minimum)))*255.0)); 
				image[x]= 255<<24 | pixel <<16| pixel <<8| pixel; 
				/*
				if (legSieve[x] == 1){   //Change shank area color to blue
					image[x]= 255<<24 | 0 <<16| 0 <<8| pixel; 
				}
				*/
				
				if (marrowSieve[x] == 1){   //Change marrow area color to green
					image[x]= 255<<24 | 0 <<16| pixel+50<<8| 0; 
					if (marrowKernelSieve[x] == 1){   //Change marrow area color to green
						image[x]= 255<<24 | pixel+50 <<16| pixel+50<<8| 0; 
					}
				}
				
				if (stSieve[x] == 1){   //Change soft Tissue area color to red
					image[x]= 255<<24 | 0 <<16| 0 <<8| pixel; 
				}
				if (muscleSieve[x] == 1){   //Change muscle area color to violet
					image[x]= 255<<24 | pixel <<16| 0 <<8| 0; 
				}
				
			 }
			 //Draw rotated radii
			for(int i = 0; i< 360;i++) {
				image[((int) (marrowCenter[0]+R[pind.get(i)]*Math.cos(Theta2[i])))+  ((int) (marrowCenter[1]+R[pind.get(i)]*Math.sin(Theta2[i])))*width]= 255<<24 | 255 <<16| 0 <<8| 255;
				image[(int) (marrowCenter[0]+R2[pind.get(i)]*Math.cos(Theta2[i]))+ ((int) (marrowCenter[1]+R2[pind.get(i)]*Math.sin(Theta2[i])))*width]=255<<24 | 0 <<16| 255 <<8| 255;
			}
			 

			 Image imageToDraw = createImage(new MemoryImageSource(width,height,image,0,width));
			 imageToDraw= imageToDraw.getScaledInstance(1000, -1, Image.SCALE_SMOOTH);
			 BufferedImage bufferedImage = (BufferedImage) showFigure.createImage(imageToDraw.getWidth(null), imageToDraw.getHeight(null));
			 Graphics2D gbuf = bufferedImage.createGraphics();
			 gbuf.drawImage(imageToDraw, 0, 0,null);

			 return bufferedImage;
		}
		
		public BufferedImage getMyImage(double[] imageIn,double[] marrowCenter,Vector<Integer> pind, double[] R, double[] R2, double[] Theta2, 
				byte[] marrowSieve,byte[] marrowKernelSieve, byte[] stSieve, byte[] muscleSieve, byte[] legSieve,
				int width, int height, double minimum, double maximum, Component imageCreator) {
			 int[] image = new int[width*height];
			 int pixel;
			 for (int x = 0; x < width*height;x++) {
				pixel = (int) (((((double) (imageIn[x] -minimum))/((double)(maximum-minimum)))*255.0)); //Korjaa tama...
				image[x]= 255<<24 | pixel <<16| pixel <<8| pixel; 

				
				if (marrowSieve[x] == 1){   //Change marrow area color to green
					image[x]= 255<<24 | 0 <<16| pixel+50<<8| 0; 
					if (marrowKernelSieve[x] == 1){   //Change marrow area color to green
						image[x]= 255<<24 | pixel+50 <<16| pixel+50<<8| 0; 
					}
				}
				
				if (stSieve[x] == 1){   //Change soft Tissue area color to red
					image[x]= 255<<24 | 0 <<16| 0 <<8| pixel; 
				}
				if (muscleSieve[x] == 1){   //Change muscle area color to violet
					image[x]= 255<<24 | pixel <<16| 0 <<8| 0; 
				}
				
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
	
	public void UKK_special(double[] scaledImage,double sleeveThreshold){
		int i,j;
		i=10;
		j=10;
		while ((j < height-1 && i < width && scaledImage[i+j*width] <sleeveThreshold) || scaledImage[i+j*width] ==0){
			i++;
			if (i == width){
				j++;
				if (j >= height-1) break;
				i = 0;
			}
			//System.out.println("X Y BMD"+i+" "+j+" "+scaledImage[i+j*width]);
		}
		//System.out.println("Sleeve coords "+i+" "+j);
		//Sleeve found
		byte[] sleeve = new byte[width*height];
		Vector<Integer> initialI = new Vector<Integer>();
		Vector<Integer> initialJ= new Vector<Integer>();
		initialI.add(i);
		initialJ.add(j);
		while (initialI.size() >0 && initialI.lastElement() > 0 &&  initialI.lastElement() < width-1 && initialJ.lastElement() > 0 && initialJ.lastElement() < height-1){
			i =initialI.lastElement();
			j = initialJ.lastElement();
			initialI.remove( initialI.size()-1);
			initialJ.remove( initialJ.size()-1);
			if (scaledImage[i+j*width] > sleeveThreshold && sleeve[i+j*width] ==0){
				sleeve[i+j*width] = 1;
			}

			if (scaledImage[i-1+j*width] > sleeveThreshold && sleeve[i-1+j*width] ==0) {
				initialI.add(new Integer(i-1));
				initialJ.add(new Integer(j));
			}

			if (scaledImage[i+1+j*width] > sleeveThreshold && sleeve[i+1+j*width] ==0) {
				initialI.add(new Integer(i+1));
				initialJ.add(new Integer(j));
			}

			if (scaledImage[i+(j-1)*width] > sleeveThreshold && sleeve[i+(j-1)*width] ==0) {
				initialI.add(new Integer(i));
				initialJ.add(new Integer(j-1));
			}

			if (scaledImage[i+(j+1)*width] > sleeveThreshold && sleeve[i+(j+1)*width] ==0) {
				initialI.add(new Integer(i));
				initialJ.add(new Integer(j+1));
			}

		}
		dilate(sleeve,(byte)1,(byte)0,(byte)2);
		dilate(sleeve,(byte)1,(byte)0,(byte)2);
		
		int removed=0;
		for (int ii =0;ii<width*height;ii++){
			if(sleeve[ii]==1){
				stSieve[ii]=0;
				legSieve[ii] = 0;
				++removed;
			}
		}
		//System.out.println("UKK_special "+removed);
		
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
		
		int len;
		i = 0;
		j = 0;

		//System.out.println("Edge tracing begins "+i+" "+j);
		while ((i < (width-1)) && (j < (height -1) )){
			len = 0;
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
			len++;
			//Loop traces the edge until the initial point is found.

			char z;
			boolean done = false;
			while(!done){
				initialI.clear();
				initialJ.clear();
				
				
				
				//Check the surrouding pixels. If last iit and jiit will be empty.
				//System.out.println("Direction");
				for (z=1;z < 9;z++){ 
					if (scaledImage[i+searchOrder[0][z]+(j+searchOrder[1][z])*width] >= threshold && scaledImage[i+searchOrder[0][z-1]+(j+searchOrder[1][z-1])*width] < threshold && result[i+searchOrder[0][z]+(j+searchOrder[1][z])*width] < 1) {
						initialI.add(new Integer(i+searchOrder[0][z]));
						initialJ.add(new Integer(j+searchOrder[1][z]));		
					}
				}

				//Check which of the possible candicaes changes direction the least
				if (previousI.size() > 1 && initialI.size() >1){
					//printf("checking angle i %d j %d\n",initialI[0],initialJ[0]);
					Vector<Double> angle = new Vector<Double>();
					double[] v1 = new double[2];
					double[] v2 = new double[2];
					int selection = 0;
					int zz;
					double smallestAngle = Math.PI*2;
					for (zz = 0;zz < initialI.size()-1;zz++);{
						v1[0] = (Integer) previousI.elementAt(0)-(Integer) previousI.elementAt(1);
						v1[1] = (Integer) previousJ.elementAt(0)-(Integer) previousJ.elementAt(1);
						v2[0] = (Integer) initialI.elementAt(zz)-(Integer) previousI.elementAt(1);
						v2[1] = (Integer) initialJ.elementAt(zz)-(Integer) previousJ.elementAt(1);
						angle.add(new Double(Math.atan2(v1[0],v1[1])-Math.atan2(v2[0],v2[1])));
						if (angle.lastElement()  < smallestAngle) {
							selection = zz;
						}
					}
					
					i = (Integer) initialI.elementAt(selection);
					j = (Integer) initialJ.elementAt(selection);
					
				} else if (initialI.size() > 0){
					i = initialI.firstElement();
					j = initialJ.firstElement();
				}
				//System.out.println("Direction defined");

				
				//Stop if no candidate pixels left
				if (result[i+j*width] == 1 || initialI.size()<1){
					done = true;
				}
				else{
					result[i+j*width] = 1;
					iit.add(new Integer(i));
					jiit.add(new Integer(j));
					previousI.add(new Integer(i));
					previousJ.add(new Integer(j));
					if (previousI.size()>2){
						previousI.remove(0) ;
						previousJ.remove(0);
					}
					len++;
				}
			}
			if (len > 0){
				length.add(new Integer(len));
				lengths.add(new Integer(len));
				if (length.size() < 2){
					beginnings.add(new Integer(0));
				}else{
					beginnings.add(new Integer(iit.size()-len));
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
		
		int len;
		i = 0;
		j = 0;

		while ((i < (width-1)) && (j < (height -1) )){
			len = 0;
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
			len++;
			char z;
			boolean done = false;
			while(!done){
				initialI.clear();
				initialJ.clear();
				for (z=1;z < 9;z++){ 
					if (scaledImage[i+searchOrder[0][z]+(j+searchOrder[1][z])*width] >= threshold && scaledImage[i+searchOrder[0][z-1]+(j+searchOrder[1][z-1])*width] < threshold && result[i+searchOrder[0][z]+(j+searchOrder[1][z])*width] < 1) {
						initialI.add(new Integer(i+searchOrder[0][z]));
						initialJ.add(new Integer(j+searchOrder[1][z]));		
					}
				}
				if (previousI.size() > 1 && initialI.size() >1){
					Vector<Double> angle = new Vector<Double>();
					double[] v1 = new double[2];
					double[] v2 = new double[2];
					int valinta = 0;
					int zz;
					double smallestAngle = Math.PI*2;
					for (zz = 0;zz < initialI.size()-1;zz++);{
						v1[0] = (Integer) previousI.elementAt(0)-(Integer) previousI.elementAt(1);
						v1[1] = (Integer) previousJ.elementAt(0)-(Integer) previousJ.elementAt(1);
						v2[0] = (Integer) initialI.elementAt(zz)-(Integer) previousI.elementAt(1);
						v2[1] = (Integer) initialJ.elementAt(zz)-(Integer) previousJ.elementAt(1);
						angle.add(new Double(Math.atan2(v1[0],v1[1])-Math.atan2(v2[0],v2[1])));
						if (angle.lastElement()  < smallestAngle) {
							valinta = zz;
						}
					}
					
					i = (Integer) initialI.elementAt(valinta);
					j = (Integer) initialJ.elementAt(valinta);
				} else if (initialI.size() > 0){
					i = initialI.firstElement();
					j = initialJ.firstElement();
				}

				if (result[i+j*width] == 1 || initialI.size()<1){
					done = true;
				}
				else{
					result[i+j*width] = 1;
					iit.add(new Integer(i));
					jiit.add(new Integer(j));
					previousI.add(new Integer(i));
					previousJ.add(new Integer(j));
					if (previousI.size()>2){
						previousI.remove(0) ;
						previousJ.remove(0);
					}
					len++;
				}
			}
			if (len > 0){
				length.add(new Integer(len));
				lengths.add(new Integer(len));
				if (length.size() < 2){
					beginnings.add(new Integer(0));
				}else{
					beginnings.add(new Integer(iit.size()-len));
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


			


		//System.out.println("Loppu found funktio");
		
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
	public static double[] interpolateMarrow(Vector<Double> interpolants){
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
		
		double[] result = new double[10];
		double inter;
		result[0]=interpolants.get(0);
		result[9]=interpolants.get(interpolants.size()-1);
		//System.out.println("Interpolate ... size "+(interpolants.size()-1));
		for (int i = 1; i<9; i++){
			inter = ((double) (i))/9.0*((double)(interpolants.size()-1));
			//System.out.println(i+" "+inter+" "+((int) Math.floor(inter)) +" ");
			result[i] = interpolants.get((int) Math.floor(inter))+bb[(int) Math.floor(inter)]*(inter-Math.floor(inter))+c.get((int) Math.floor(inter),0)*Math.pow(inter-Math.floor(inter),2.0)+d[(int) Math.floor(inter)]*Math.pow(inter-Math.floor(inter),3.0);
		}
		//System.out.println("Interpolate ... done");
		return result;
	}
}
