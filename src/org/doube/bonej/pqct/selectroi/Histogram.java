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

public class Histogram extends JPanel{
	ImageAndAnalysisDetails details;
	public double[] scaledImage;
	public double[] coefficients;
	public short[] histoImage;
	public short[] lbpImage;
	public double[] lbpMatch;
	public short[] histoSieve;
	public int histoWidth;
	public int step;
	public double[] softScaledImage;
	public int[] LBPmatchCoords;
	public double[] cortexROI;
	public int[] peaks;
	public int[] histoScaledImageSieve;
	public double[] calibObjectHUs;
	public double minimum;
	public double maximum;
	public Vector<Integer> iit;
	public Vector<Integer> jiit;
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
	public int[] histo;
	public int[] dummy;
	public int peak;

	public int calibHeight;
	public int calibWidth;
	
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
	public WritableRaster imageP;

	/*ROI Search with LBP and histogram combined*/
	public Histogram(ScaledImageData dataIn,ImageAndAnalysisDetails detailsIn, DrawImage showFigureIn,String imageSaveNameIn,double[] histoLBPTemplate,double[] histoGrayTemplate, int histolevels){
		
		details =detailsIn;
		scaledImage = (double[])dataIn.scaledImage.clone();
		
		//softScaledImage = (double[])dataIn.softScaledImage.clone();
		showFigure = showFigureIn;
		pixelSpacing = dataIn.pixelSpacing;
		imageSavePath = details.imageSavePath;
		imageSaveName = imageSaveNameIn;
		width =dataIn.width;
		height =dataIn.height;
		histoImage = new short[height*width];
		lbpImage= new short[height*width];
		
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
		if (maximum >2600){maximum = 2600;}
		//calibHeight = 40;	//Height of the window used for searching the phantoms
		//calibWidth = 510;	//Height of the window used for searching the phantoms 130 260
		calibHeight = 25;	//Height of the window used for searching the phantoms
		calibWidth = 25;	//Height of the window used for searching the phantoms 130 260
		step = 8;
		lbpMatch = new double[height/step*width/step];
		
		//Select calibration ROI using histogram feature vector
		//Normalize image histogram to 256 levels and get the histogram
		//System.out.println("Normalize");
		histoSieve = new short[width*height];		
		for (int j = 0; j<height;j++){ 
			for (int i = 0; i<width;i++){ 
				histoImage[i+j*width] = (short) (Math.floor(((scaledImage[i+j*width]-minimum)/(maximum-minimum+1))*((double)histolevels)));
				
				
				
				
				if (histoImage[i+j*width] < (short) (histolevels*0.33)){
					histoImage[i+j*width] = 0;
					histoSieve[i+j*width] = 0;
				}else{
					histoSieve[i+j*width] = 1;
				}
			}
		}

		//remove bone&muscle
		int poistoja =0;
		
		for (int j = 0; j<height;j++){ 
			for (int i = 0; i<width;i++){ 
				if (histoSieve[i+j*width] > 0 && histoImage[i+j*width] > (short) (histolevels*0.85)){
					fillConnected(histoImage, (short) (histolevels*0.33-1), histoSieve, i, j);		
					poistoja++;
				}
			}
		}
		
		
		
		erodeHistoImage(histoSieve,(short) 0);
		erodeHistoImage(histoSieve,(short) 0);
		/*
		erodeHistoImage(histoSieve,(short) 0);
		erodeHistoImage(histoSieve,(short) 0);
		erodeHistoImage(histoSieve,(short) 0);*/
		
		/*
		dilateHistoImage(histoSieve,(short) 0);
		dilateHistoImage(histoSieve,(short) 0);
		dilateHistoImage(histoSieve,(short) 0);
		
		dilateHistoImage(histoSieve,(short) 0);
		dilateHistoImage(histoSieve,(short) 0);
		*/
		for (int i = 0; i< width*height;++i){
			if (histoSieve[i] == 0){histoImage[i] = 0;}
		}
		/*Go through the whole image in 5 pixel steps. Highest proximity histogram is the phantoms and will be selected.. */
		//lbp_figure(histoImage, lbpImage);
		lbp_figureRI(histoImage, lbpImage);
		double bestLBPMatch = 0;
		LBPmatchCoords = new int[2];
		
		int[] histoLBP = new int[10];
		int[] histoGray = new int[histolevels];
		//System.out.println("Proximity");
		
		for (int j = 0; j<height-calibHeight;j=j+step){ 
			for (int i = 0; i<width-calibWidth;i=i+step){ 
				//Calculate the histogram of the rectangle
				histoLBP = new int[10];
				histoGray = new int[histolevels];
				for (int jj=0; jj<calibHeight;jj++){
					for (int ii=0; ii<calibWidth;ii++){
						++histoLBP[(int) lbpImage[(i+ii)+(j+jj)*width]];
						++histoGray[(int) histoImage[(i+ii)+(j+jj)*width]];
					}
				}
				
				//normalize the histogram and calculate proximity to histoTemplate
				double[] histoLBPNormalized = new double[10];
				double[] histoGrayNormalized = new double[histolevels];
				double LBPsum =0;
				double Graysum =0;
				double corrSum =0;
				
				for (int n = 0; n<10;n++){
					histoLBPNormalized[n] = ((double)(histoLBP[n]))/((double)(calibHeight*calibWidth));
					//Calculate LBP proximity
					LBPsum+=Math.min(histoLBPNormalized[n],histoLBPTemplate[n]);
				}
				
				//Calculate proximity after first determining lag with xcor...
				
				for (int n = 0; n<histolevels;n++){
					histoGrayNormalized[n] = ((double)(histoGray[n]))/((double)(calibWidth*calibHeight));
				}
				
				int maxLag = (int) (0.14*histolevels);
				double[] crossCorr = xcorr(histoGrayNormalized,histoGrayTemplate, histolevels, maxLag);
				int optimumLagIndex;
				optimumLagIndex = maxIndex(crossCorr,maxLag);
				
				
				corrSum = max(crossCorr,maxLag);
				
				
				//Calculate histogram proximity, take zeroth index to account for amount of background...
				if (optimumLagIndex != -maxLag){
					Graysum+=Math.min(histoGrayNormalized[0],histoGrayTemplate[0]);
				}
				for (int n = maxLag+optimumLagIndex; n<histolevels-maxLag+optimumLagIndex;n++){
					//Calculate histogram proximity from the matched histogram...
					Graysum+=Math.min(histoGrayNormalized[n],histoGrayTemplate[n]);
				}
				
				lbpMatch[(int) Math.floor(i/step+(j/step)*(width/step))] = LBPsum+Graysum+corrSum;
				
			}
		}
		
		/*Find the roi from lbpMatch with the highest areal proximity sum*/
		
		calibHeight = 60;
		calibWidth = 260;
		int[] testi = new int[2];
		double areaSum = 0;
		for (int j = 0; j<height/step-calibHeight/step;j++){
			for (int i = 0; i<width/step-calibWidth/step;i++){
				double aSum = 0;
				for (int jj=0; jj<calibHeight/step;jj++){
					for (int ii=0; ii<calibWidth/step;ii++){
						aSum +=lbpMatch[(i+ii)+(j+jj)*width/step];
					}
				}
				if (aSum > areaSum){
					areaSum = aSum;
					LBPmatchCoords[0] = i*step;	//X-coordinate of the upper left corner
					LBPmatchCoords[1] = j*step;	//Y-coordinate of the upper left corner
					testi[0] = i;	//X-coordinate of the upper left corner
					testi[1] = j;	//Y-coordinate of the upper left corner
				}
			}
		}
		
		/*
		for (int j = testi[1]; j<testi[1]+calibHeight/step;j++){
			for (int i = testi[0]; i<testi[0]+calibWidth/step;i++){
				lbpMatch[i+j*width/step] = 0;
			}
		}
		*/
		//System.out.println("Best match sum "+bestLBPMatch);
		/*Visualize result, remove later on...*/
		//LBPmatchCoords[0] = 130;
		//LBPmatchCoords[1] = 400;
		short lisa = (short) (0.14*histolevels);
		histoLBP = new int[10];
		histoGray = new int[histolevels];
		for (int j = 0; j<height;j++){ 
			for (int i = 0; i<width;i++){ 
				if ((i < LBPmatchCoords[0] || i > LBPmatchCoords[0]+calibWidth) || (j < LBPmatchCoords[1] || j > LBPmatchCoords[1]+calibHeight)){
					//lbpImage[i+j*width] = 0;
					if (histoImage[i+j*width] > (0.14*histolevels)){
						histoImage[i+j*width] = (short) (((int) histoImage[i+j*width])-lisa);
						//scaledImage[i+j*width] -=200;
					}
				}
				else{
					++histoLBP[(int) lbpImage[i+j*width]];
					++histoGray[(int) histoImage[i+j*width]];
					if (histoImage[i+j*width] < (0.86*histolevels)){
						histoImage[i+j*width] = (short)(((int) histoImage[i+j*width])+lisa);
						//scaledImage[i+j*width] +=500;
					}
				}
			
			}
		}
		
		/* //Writing histograms...
		try{//Get feature vector to be matched..
			FileWriter resultFile = new FileWriter("C:/Oma/Deakin/GENTS/LBPHistogrammit.xls",true);	//Open for appending
			BufferedWriter writer = new BufferedWriter(resultFile);
			for (int i = 0; i< 10; i++){
				writer.write(((double)(histoLBP[i]))/((double)(calibWidth*calibHeight))+"\t");
			}
			writer.write("\n");
			//Close the output stream
			writer.close();
		}catch (Exception err){//Catch exception if any
			System.err.println("Error: " + err.getMessage());
		}
		try{//Get feature vector to be matched..
			FileWriter resultFile = new FileWriter("C:/Oma/Deakin/GENTS/GrayHistogrammit.xls",true);	//Open for appending
			BufferedWriter writer = new BufferedWriter(resultFile);
			for (int i = 0; i< histolevels; i++){
				writer.write(((double)(histoGray[i]))/((double)(calibWidth*calibHeight))+"\t");
			}
			writer.write("\n");
			//Close the output stream
			writer.close();
		}catch (Exception err){//Catch exception if any
			System.err.println("Error: " + err.getMessage());
		}
		*/
		//ROI selected
		
		Vector<Double> bins = new Vector<Double>();
		histoWidth = 2000; //Used to be 1000
		int[] histo2 = new int[histoWidth];
		/*
		for (int i = 0; i < histoWidth; i++){
			bins.add((double) i/((double)histoWidth)*(maximum-minimum)+minimum);
		}
		*/
		peak =0;
		int indexOfPeak =0;
		
		histoScaledImageSieve = new int[width*height];
		// Calculate the histogram for the original image for obtaining the calibration equation
		for (int j = LBPmatchCoords[1]; j<LBPmatchCoords[1]+calibHeight;j++){ 
			for (int i = LBPmatchCoords[0]; i<LBPmatchCoords[0]+calibWidth;i++){ 
				if (histoSieve[i+j*width] > 0){
					histo2[(int) (Math.floor(((scaledImage[i+j*width]-minimum)/(maximum-minimum)*((double)histoWidth))))]+=1;
					histoScaledImageSieve[i+j*width] =(int) (Math.floor(((scaledImage[i+j*width]-minimum)/(maximum-minimum)*((double)histoWidth))));
					if (histo2[(int) (Math.floor(((scaledImage[i+j*width]-minimum)/(maximum-minimum)*((double)histoWidth))))] > peak && ((int) (Math.floor(((scaledImage[i+j*width]-minimum)/(maximum-minimum)*((double)histoWidth))))) > 0){
						peak = histo2[(int) (Math.floor(((scaledImage[i+j*width]-minimum)/(maximum-minimum)*((double)histoWidth))))];
						indexOfPeak =(int) (Math.floor(((scaledImage[i+j*width]-minimum)/(maximum-minimum)*((double)histoWidth))));
					}
				}
			}
		}
		/*Moving average for the histogram*/
		
		histo = new int[histoWidth];
		int movingAverageWidth = 8;
		for (int i =movingAverageWidth;i<histoWidth-movingAverageWidth;i++){
			double vali = 0;
			for (int j = -movingAverageWidth;j<=movingAverageWidth;++j){
				 vali+=histo2[i+j];
			}
			vali = vali/(((double)movingAverageWidth*2.0)+1.0);
			histo[i] =(int) Math.round(vali);
		}
		
		/*Work out the regions of interest of the individual calibration phantoms*/
		peaks = new int[12];

		try{
			peaks = findPeaks(histo,6);		//Returns 6 highest histogram peaks, should be the calibration phantoms...
			//System.out.println(peaks[0]+" "+peaks[1]+" "+peaks[2]+" "+peaks[3]+" "+peaks[4]+" "+peaks[5]+" "+peaks[6]+" "+peaks[7]+" "+peaks[8]+" "+peaks[9]+" "+peaks[10]+" "+peaks[11]);
		 }catch(Exception err){//Catch exception if any
			System.err.println("Piikit: " + err.getMessage());

		}

		
		/*Work out the calibration!!!*/
		
		calibObjectHUs = new double[6];

		for (int i =0;i<6;i++){
			byte[] tempSieve = new byte[width*height];
			for (int j=0;j<width*height;j++){
				if (histoScaledImageSieve[j] > peaks[2*i] && histoScaledImageSieve[j] < peaks[2*i+1]){
					tempSieve[j]=1;
				}else{
					tempSieve[j]=0;
				}
			}
			erodeBW(tempSieve);
			//erodeBW(tempSieve); /*Added another erode...*/
			int noOfValues =0;
			for (int j=0;j<width*height;j++){
				if(tempSieve[j] > 0){++noOfValues;}
			}	
			double[] values = new double[noOfValues];
			noOfValues =0;
			for (int j=0;j<width*height;j++){
				if(tempSieve[j] > 0){
					values[noOfValues]=scaledImage[j];
					++noOfValues;
				}
			}
			if (noOfValues > 0){
				calibObjectHUs[i] = median(values);
			}
		}
		/*Calculate the slope and the intercept...*/
		double[] phantomHUs = {0,50,100,150,250,350};
		
		coefficients = linearFit(calibObjectHUs,phantomHUs);
		
		//DRAW FIGURES
		if (details.imOn ==true){
			try{
				/*Try to draw several images to the same one...*/
				 BufferedImage img = new BufferedImage (1000, 1000+300,BufferedImage.TYPE_INT_RGB);	/*The "canvas"*/
			
				 //Save ROI selection image
				BufferedImage bi = getMyImage(histo,peak,peaks); // retrieve histo image
				BufferedImage bi2 = getMyImage(scaledImage,width,height,minimum, maximum,LBPmatchCoords,calibWidth, calibHeight, histoScaledImageSieve, peaks); // retrieve image
				img.createGraphics().drawImage(bi2, 0, 0, null);
				img.createGraphics().drawImage(bi, 0, 1000, null);
				File imageOutputFile = new File(details.imageSavePath+"/"+ imageSaveName+".png");
				ImageIO.write(img, "png", imageOutputFile);
			}catch (Exception err){System.err.println("Couldn't write image out: "+err.getMessage());}
		}
		
	}
	
	/*Linear Fit*/
	public double[] linearFit(double[] knownX, double[] knownY){
		double[][] xx = new double[knownX.length][2];
		
		for (int i = 0; i< knownX.length;i++){
			xx[i][0] =1;
			xx[i][1] =knownX[i];
		}
		Matrix Y = new Matrix(knownY,knownY.length);
		Matrix X = new Matrix(xx);
		double[] result = new double[2];
		try{
			Matrix coefficients = X.solve(Y);
			result[1] = coefficients.get(1,0);
			result[0] = coefficients.get(0,0);
		}catch(Exception err){
			System.err.println("Matrix coefficients: " + err.getMessage());
			return result;
		}
		
		return result;
	}
	
	void erodeBW(byte[] data){
		//Erode algorithm
		//Modified from the best dilate by one solution taken from http://ostermiller.org/dilate_and_erode.html
		byte background = 0;
		for (int i=0; i<height; i++){
			for (int j=0; j<width; j++){
				if (data[i*width+j] > background){
					if ((i>0 && data[(i-1)*width+j]<=background && data[(i-1)*width+j] > -1) ||
						(j>0 && data[(i)*width+j-1]<=background) && data[(i)*width+j-1] > -1 ||
						(i+1<height && data[(i+1)*width+j]<=background && data[(i+1)*width+j] > -1)||
						(j+1<width && data[(i)*width+j+1]<=background && data[(i)*width+j+1] > -1))
						{data[i*width+j] = -1;}	//Erode the pixel if any of the neighborhood pixels is background
				}
			}
		}
		for (int i=0; i<width*height; i++){
			if (data[i] < 0){
				data[i] = 0;
			}
		}
	}
	/*From http://www.leepoint.net/notes-java/data/arrays/30arraysExamples.html*/
	//================================================== median
	//   Precondition: Array must be sorted ADDED the sort...
	public static double median(double[] m) {
		Arrays.sort(m);
		int middle = m.length/2;  // subscript of middle element
		if (m.length%2 == 1) {
			// Odd number of elements -- return the middle one.
			return m[middle];
		} else {
		   // Even number -- return average of middle two
		   // Must cast the numbers to double before dividing.
		   return (m[middle-1] + m[middle]) / 2.0;
		}
	}//end method median

	
	
	public int maxInd(int[] peaks){
		int[] sorted = (int[])peaks.clone();
		Arrays.sort(sorted);
		int iiii = 0;
		
		while (peaks[iiii] != sorted[sorted.length-1]){
			++iiii;
		}
		
		return iiii;	
	}
	public int max(int[] peaks){
		int[] sorted = (int[])peaks.clone();
		Arrays.sort(sorted);
		return sorted[sorted.length-1];	
	}
	int[] findPeaks(int[] histo,int howMany){
		int[] peaks = new int[howMany*2];
		int[] temp = (int[]) histo.clone();
		temp[0] = 0; //Remove background...
		
		int peak=0;
		int maximum = 0;
		for (int i = 0;i<howMany;i++){
			peak = maxInd(temp);
			maximum = max(temp);
			if (peak > 0 && maximum > 0){
				int j = peak;
				int previous = maximum;
				//while (previous >= temp[j] && temp[j] > maximum/2){
				while (temp[j] > maximum/2){
					previous = temp[j];
					--j;
					if (j == 0){return peaks;}
				}
				peaks[i*2] = j;
				j = peak;
				previous = maximum;
				//while (previous >= temp[j] && temp[j] > maximum/2){
				while (temp[j] > maximum/2){
					previous = temp[j];
					++j;
					if (j == histo.length){return peaks;}
				}
				peaks[i*2+1] = j;
				for (int k = peaks[i*2]-5; k <=peaks[i*2+1]+5;++k){
					temp[k] = 0;
				}
			}
		}
		Arrays.sort(peaks);	//Return peaks in ascending order...
		return peaks;
	}
	
	void fillConnected(short[] hist, short thresh, short[] siev, int i, int j){
		Vector<Short> initialI = new Vector<Short>();
		Vector<Short> initialJ = new Vector<Short>();
		initialI.add((short)i);
		initialJ.add((short)j);
		while (initialI.size()>0){
			i =initialI.lastElement();
			j =initialJ.lastElement();
			initialI.remove(initialI.size()-1);
			initialJ.remove(initialJ.size()-1);
			
			if (hist[i+j*width] >= thresh && siev[i+j*width]==1){
				siev[i+j*width] = 0;
				//hist[i+j*width] = 0;
			}
			
			if (i>0 && hist[i-1+j*width] >= thresh && siev[i-1+j*width]==1) {
				initialI.add((short)(i-1));
				initialJ.add((short)j);
			}
		
			if (i < width && hist[i+1+j*width] >= thresh && siev[i+1+j*width]==1) {
				initialI.add((short)(i+1));
				initialJ.add((short)j);
			}
		
			if (j>0 && hist[i+(j-1)*width] >= thresh && siev[i+(j-1)*width]==1) {
				initialI.add((short)i);
				initialJ.add((short)(j-1));
			}
		
			if (j < height && hist[i+(j+1)*width] >= thresh && siev[i+(j+1)*width]==1) {
				initialI.add((short)i);
				initialJ.add((short)(j+1));
			}
		}
	}


	public int maxIndex(double[] crossCorr, int length){
		double[] sorted = (double[])crossCorr.clone();
		Arrays.sort(sorted);
		int i = 0;
		while (crossCorr[i] != sorted[length*2]){
			++i;
		}
		i -=length;
		return i;	
	}
	
	public double max(double[] crossCorr, int length){
		double[] sorted = (double[])crossCorr.clone();
		Arrays.sort(sorted);
		return sorted[length*2];	
	}
	/*Calculate cross-correlation, equations taken from http://paulbourke.net/miscellaneous/correlate/*/
	public double[] xcorr(double[] series1,double[] series2, int length, int maxDelay){
		double[] xcor = new double[maxDelay*2+1];
		double ms1 =0;
		double ms2 =0;
		/*calculate means*/
		for (int i = 0; i<length;i++){
			ms1+=series1[i];
			ms2+=series2[i];
		}
		ms1 /=((double)length);
		ms2 /=((double)length);
		double mx;
		double my;
		double summxmy;
		double summxSq;
		double summySq;
		double summxmySq;
		for (int i =-maxDelay;i<=maxDelay;i++){//ignore beginning and end of the signal...
			summxmy=0;
			summxSq=0;
			summySq=0;
			for (int j = maxDelay; j< length-maxDelay; j++){
				mx = series1[j]-ms1;
				my = series2[j+i]-ms1;
				summxmy+=mx*my;
				summxSq+=mx*mx;
				summySq+=my*my;
			}
			xcor[i+maxDelay]=summxmy/Math.sqrt(summxSq*summySq);
		}
		return xcor;
	}
	
	/*Rotation invariant LBP...
	static unsigned char mapping[256] = {0, 1, 1, 2, 1, 9, 2, 3, 1, 9, 9, 9, 2, 9, 3, 4, 1, 9, 9, 9, 9, 9, 9, 9, 2, 9, 9, 9, 3, 9, 4, 5, 1, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 2, 9, 9, 9, 9, 9, 9, 9, 3, 9, 9, 9, 4, 9, 5, 6, 1, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 2, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 3, 9, 9, 9, 9, 9, 9, 9, 4, 9, 9, 9, 5, 9, 6, 7, 1, 2, 9, 3, 9, 9, 9, 4, 9, 9, 9, 9, 9, 9, 9, 5, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 6, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 7, 2, 3, 9, 4, 9, 9, 9, 5, 9, 9, 9, 9, 9, 9, 9, 6, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 7, 3, 4, 9, 5, 9, 9, 9, 6, 9, 9, 9, 9, 9, 9, 9, 7, 4, 5, 9, 6, 9, 9, 9, 7, 5, 6, 9, 7, 6, 7, 7, 8};*/
	/*Returns LBP-figure for 3x3 neighborhood*/
	void lbp_figureRI(short[] img, short[] result)
	{
		short mapping[] = {0, 1, 1, 2, 1, 9, 2, 3, 1, 9, 9, 9, 2, 9, 3, 4, 1, 9, 9, 9, 9, 9, 9, 9, 2, 9, 9, 9, 3, 9, 4, 5, 1, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 2, 9, 9, 9, 9, 9, 9, 9, 3, 9, 9, 9, 4, 9, 5, 6, 1, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 2, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 3, 9, 9, 9, 9, 9, 9, 9, 4, 9, 9, 9, 5, 9, 6, 7, 1, 2, 9, 3, 9, 9, 9, 4, 9, 9, 9, 9, 9, 9, 9, 5, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 6, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 7, 2, 3, 9, 4, 9, 9, 9, 5, 9, 9, 9, 9, 9, 9, 9, 6, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 7, 3, 4, 9, 5, 9, 9, 9, 6, 9, 9, 9, 9, 9, 9, 9, 7, 4, 5, 9, 6, 9, 9, 9, 7, 5, 6, 9, 7, 6, 7, 7, 8};
		short value;
		//int[][] searchOrder = {{1,0,-1,-1,-1,0,1,1},{1,1,1,0,-1,-1,-1,0}};
		int[][] searchOrder = {{-1,-1,-1,0,1,1,1,0},{-1,0,1,1,1,0,-1,-1}};
		//int[][] searchOrder = {{-2,-2,-2,0,2,2,2,0},{-2,0,2,2,2,0,-2,-2}};
		short kertoimet[] = {1,2,4,8,16,32,64,128};
		for (int r=2;r<height-2;r++){
			for (int c=2;c<width-2;c++){
				value = 0;
				for (int z=0;z < 8;z++){ 
					if (img[(r+searchOrder[0][z])*width+c+searchOrder[1][z]] >= img[r*width+c]+20) {
						value = (short) (((short)value)+((short) kertoimet[z]));
					}
				}
				result[r*width+c] = mapping[value]; /* Increase histogram bin value */
			}
		}
	}
	
	
	/*Returns LBP-figure for 3x3 neighborhood*/
	void lbp_figure(short[] img, short[] result)
	{
		short value;
		int[][] searchOrder = {{1,0,-1,-1,-1,0,1,1},{1,1,1,0,-1,-1,-1,0}};
		short kertoimet[] = {1,2,4,8,16,32,64,128};
		for (int r=1;r<height-1;r++){
			for (int c=1;c<width-1;c++){
				value = 0;
				for (int z=0;z < 8;z++){ 
					if (img[(r+searchOrder[0][z])*width+c+searchOrder[1][z]] >= img[r*width+c]+20) {
						value = (short) (((short)value)+((short) kertoimet[z]));
					}
				}
				result[r*width+c] = value; /* Increase histogram bin value */
			}
		}
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
	
	
	void dilateHistoImage(short[] data,short background){
		//Dilate algorithm
		short temp = -1;
		//Modified from the best dilate by one solution taken from http://ostermiller.org/dilate_and_erode.html
		for (int i=0; i<height; i++){
			for (int j=0; j<width; j++){
				if (data[i*width+j] > background){
					if (i>0 && data[(i-1)*width+j]==background) {data[(i-1)*width+j] = temp;}
					if (j>0 && data[(i)*width+j-1]==background) {data[(i)*width+j-1] = temp;}
					if (i+1<height && data[(i+1)*width+j]==background) {data[(i+1)*width+j] = temp;}
					if (j+1<width && data[(i)*width+j+1]==background) {data[(i)*width+j+1] = temp;}
				}
			}
		}
		for (int i=0; i<width*height; i++){
			if (data[i] ==temp){
				data[i] = 1;
			}
		}
	}
	
	void erodeHistoImage(short[] data,short background){
		//Erode algorithm
		//Modified from the best dilate by one solution taken from http://ostermiller.org/dilate_and_erode.html
		for (int i=0; i<height; i++){
			for (int j=0; j<width; j++){
				if (data[i*width+j] > background){
					if ((i>0 && data[(i-1)*width+j]<=background && data[(i-1)*width+j] > -1) ||
						(j>0 && data[(i)*width+j-1]<=background) && data[(i)*width+j-1] > -1 ||
						(i+1<height && data[(i+1)*width+j]<=background && data[(i+1)*width+j] > -1)||
						(j+1<width && data[(i)*width+j+1]<=background && data[(i)*width+j+1] > -1))
						{data[i*width+j] = -1;}	//Erode the pixel if any of the neighborhood pixels is background
				}
			}
		}
		for (int i=0; i<width*height; i++){
			if (data[i] < 0){
				data[i] = 0;
			}
		}
	}
	
	
	void erodeHistoImage(byte[] data,byte background){
		//Erode algorithm
		//Modified from the best dilate by one solution taken from http://ostermiller.org/dilate_and_erode.html
		for (int i=0; i<height; i++){
			for (int j=0; j<width; j++){
				if (data[i*width+j] > background){
					if ((i>0 && data[(i-1)*width+j]<=background && data[(i-1)*width+j] > -1) ||
						(j>0 && data[(i)*width+j-1]<=background) && data[(i)*width+j-1] > -1 ||
						(i+1<height && data[(i+1)*width+j]<=background && data[(i+1)*width+j] > -1)||
						(j+1<width && data[(i)*width+j+1]<=background && data[(i)*width+j+1] > -1))
						{data[i*width+j] = -1;}	//Erode the pixel if any of the neighborhood pixels is background
				}
			}
		}
		for (int i=0; i<width*height; i++){
			if (data[i] < 0){
				data[i] = 0;
			}
		}
	}
	

	public BufferedImage getMyImage(int[] image,int[] histo,int peak, int[] peaks) {
			 int pixel;
			 for (int x = 1; x < histoWidth;x++) {
				for (int y = 299;y>=299-(int)((((double) histo[x])/((double)peak))*299.0);--y){
					image[x+y*histoWidth]= 255<<24 | 0 <<16| 0 <<8| 0; 
					if (x > peaks[0] && x < peaks[1]){image[x+y*histoWidth]= 255<<24 | 255 <<16| 0 <<8| 0;} 
					if (x > peaks[2] && x < peaks[3]){image[x+y*histoWidth]= 255<<24 | 0 <<16| 255 <<8| 0;} 
					if (x > peaks[4] && x < peaks[5]){image[x+y*histoWidth]= 255<<24 | 0 <<16| 0 <<8| 255;} 
					if (x > peaks[6] && x < peaks[7]){image[x+y*histoWidth]= 255<<24 | 255 <<16| 255 <<8| 0;} 
					if (x > peaks[8] && x < peaks[9]){image[x+y*histoWidth]= 255<<24 | 0 <<16| 255 <<8| 255;} 
					if (x > peaks[10] && x < peaks[11]){image[x+y*histoWidth]= 255<<24 | 255 <<16| 0 <<8| 255;} 
				}
			 }

			 Image imageToDraw = createImage(new MemoryImageSource(histoWidth,300,image,0,histoWidth));
			 imageToDraw= imageToDraw.getScaledInstance(1000, 300, Image.SCALE_SMOOTH);
			 
			 BufferedImage bufferedImage = (BufferedImage) showFigure.createImage(imageToDraw.getWidth(null), imageToDraw.getHeight(null));
			 
			 Graphics2D gbuf = bufferedImage.createGraphics();
			 
			 gbuf.drawImage(imageToDraw, 0, 0,null);
			 
			 return bufferedImage;
		}
		public BufferedImage getMyImage(int[] histo,int peak, int[] peaks) {
			 int pixel;
			 int[] image = new int[histoWidth*300];
			 for (int x = 1; x < histoWidth;x++) {
				for (int y = 299;y>=299-(int)((((double) histo[x])/((double)peak))*299.0);--y){
					image[x+y*histoWidth]= 255<<24 | 0 <<16| 0 <<8| 0; 
					if (x > peaks[0] && x < peaks[1]){image[x+y*histoWidth]= 255<<24 | 255 <<16| 0 <<8| 0;} 
					if (x > peaks[2] && x < peaks[3]){image[x+y*histoWidth]= 255<<24 | 0 <<16| 255 <<8| 0;} 
					if (x > peaks[4] && x < peaks[5]){image[x+y*histoWidth]= 255<<24 | 0 <<16| 0 <<8| 255;} 
					if (x > peaks[6] && x < peaks[7]){image[x+y*histoWidth]= 255<<24 | 255 <<16| 255 <<8| 0;} 
					if (x > peaks[8] && x < peaks[9]){image[x+y*histoWidth]= 255<<24 | 0 <<16| 255 <<8| 255;} 
					if (x > peaks[10] && x < peaks[11]){image[x+y*histoWidth]= 255<<24 | 255 <<16| 0 <<8| 255;} 
				}
			 }

			 Image imageToDraw = createImage(new MemoryImageSource(histoWidth,300,image,0,histoWidth));
			 imageToDraw= imageToDraw.getScaledInstance(1000, 300, Image.SCALE_SMOOTH);
			 
			 BufferedImage bufferedImage = (BufferedImage) showFigure.createImage(imageToDraw.getWidth(null), imageToDraw.getHeight(null));
			 
			 Graphics2D gbuf = bufferedImage.createGraphics();
			 
			 gbuf.drawImage(imageToDraw, 0, 0,null);
			 
			 return bufferedImage;
		}
		
	public BufferedImage getMyImage(int[] image,int[] histo,int peak) {
			 int pixel;
			 for (int x = 1; x < 1000;x++) {
				for (int y = 299;y>=299-(int)((((double) histo[x])/((double)peak))*299.0);--y){
					image[x+y*1000]= 255<<24 | 0 <<16| 255 <<8| 0; 
					
				}
				image[x]= 255<<24 | 255 <<16| 0 <<8| 0; 
			 }

			 Image imageToDraw = createImage(new MemoryImageSource(1000,300,image,0,1000));

			 
			 BufferedImage bufferedImage = (BufferedImage) showFigure.createImage(imageToDraw.getWidth(null), imageToDraw.getHeight(null));
			 
			 Graphics2D gbuf = bufferedImage.createGraphics();
			 
			 gbuf.drawImage(imageToDraw, 0, 0,null);
			 
			 return bufferedImage;
		}
	
    
		public BufferedImage getMyImage(double[] imageIn,int width, int height,double min, double max,int[] coords, int calibWidth, int calibHeight,int[] histoScaledImageSieve,int[] peaks) {
			int[] image = new int[width*height];
		int pixel;
		for (int x = 0; x < width*height;x++) {
         pixel = (int) (((((double) (imageIn[x] -min))/((double)(max-min)))*255.0));
			image[x]= 255<<24 | pixel <<16| pixel <<8| pixel; 
		}
		for (int j = coords[1]; j<coords[1]+calibHeight;j++){ 
			for (int i = coords[0]; i<coords[0]+calibWidth;i++){ 
				pixel = (int) (((((double) (imageIn[i+j*width] -min))/((double)(max-min)))*255.0));
				image[i+j*width] = 255<<24 | 128 <<16| pixel <<8| 0; 

			}
		}
		/*Mark calibration phantoms*/
		for (int x = 0; x < width*height;x++){
			
			if (histoScaledImageSieve[x] > peaks[0] && histoScaledImageSieve[x] < peaks[1]){image[x]= 255<<24 | 255 <<16| 0 <<8| 0;} 
			if (histoScaledImageSieve[x] > peaks[2] && histoScaledImageSieve[x] < peaks[3]){image[x]= 255<<24 | 0 <<16| 255 <<8| 0;} 
			if (histoScaledImageSieve[x] > peaks[4] && histoScaledImageSieve[x] < peaks[5]){image[x]= 255<<24 | 0 <<16| 0 <<8| 255;} 
			if (histoScaledImageSieve[x] > peaks[6] && histoScaledImageSieve[x] < peaks[7]){image[x]= 255<<24 | 255 <<16| 255 <<8| 0;} 
			if (histoScaledImageSieve[x] > peaks[8] && histoScaledImageSieve[x] < peaks[9]){image[x]= 255<<24 | 0 <<16| 255<<8| 255;} 
			if (histoScaledImageSieve[x] > peaks[10] && histoScaledImageSieve[x] < peaks[11]){image[x]= 255<<24 | 255 <<16| 0 <<8| 255;} 
		}

			 Image imageToDraw = createImage(new MemoryImageSource(width,height,image,0,width));
			 imageToDraw= imageToDraw.getScaledInstance(1000, -1, Image.SCALE_SMOOTH);
			 
			 BufferedImage bufferedImage = (BufferedImage) showFigure.createImage(imageToDraw.getWidth(null), imageToDraw.getHeight(null));
			 
			 Graphics2D gbuf = bufferedImage.createGraphics();
			 
			 gbuf.drawImage(imageToDraw, 0, 0,null);
			 
			 return bufferedImage;
		}

		public BufferedImage getMyImage(double[] imageIn,Vector<Integer> xx,Vector<Integer> yy, byte[] marrowSieve,byte[] marrowKernelSieve, byte[] stSieve, byte[] muscleSieve, byte[] legSieve) {
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

			 Image imageToDraw = createImage(new MemoryImageSource(width,height,image,0,width));
			 imageToDraw= imageToDraw.getScaledInstance(1000, -1, Image.SCALE_SMOOTH);
			 
			 BufferedImage bufferedImage = (BufferedImage) showFigure.createImage(imageToDraw.getWidth(null), imageToDraw.getHeight(null));
			 
			 Graphics2D gbuf = bufferedImage.createGraphics();
			 
			 gbuf.drawImage(imageToDraw, 0, 0,null);
			 
			 return bufferedImage;
		}


		
		public BufferedImage getMyImage(double[] imageIn,double[] marrowCenter,Vector<Integer> pind, double[] R, double[] R2, double[] Theta2, 
				byte[] marrowSieve,byte[] marrowKernelSieve, byte[] stSieve, byte[] muscleSieve, byte[] legSieve,
				int width, int height, double minimum, double maximum) {
			 int[] image = new int[width*height];
			 int pixel;
			 for (int x = 0; x < width*height;x++) {
				pixel = (int) (((((double) (imageIn[x] -minimum))/((double)(maximum-minimum)))*255.0)); //Korjaa tama...
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
