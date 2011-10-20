package org.doube.bonej.pqct;

import ij.*;
import ij.text.*;
import ij.process.*;
import ij.gui.*;
import java.util.*;	//Vector
import ij.plugin.filter.*;
import org.doube.bonej.pqct.analysis.*;		//Analysis stuff..
import org.doube.bonej.pqct.selectroi.*;		//ROI selection..
import org.doube.bonej.pqct.io.*;	//image data 

import java.awt.image.*; //Creating the result BufferedImage...
import ij.plugin.filter.Info;
import ij.io.*;
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

	ImageJ density distribution analysis plugin
    Copyright (C) 2011 Timo Rantalainen
*/


public class Distribution_Analysis implements PlugInFilter {
	ImagePlus imp;

	int sectorWidth;
	boolean cOn;
	boolean dOn;
	String resultString;
	String imageInfo;
	double resolution;
	double fatThreshold;
	double areaThreshold;
	double BMDThreshold;
	double scalingFactor;
	double constant;	
	boolean flipDistribution;
	boolean manualRotation;
	public int setup(String arg, ImagePlus imp) {
		this.imp = imp;
		//return DOES_32;
		return DOES_ALL;
	}

	/*
	//For Debugging
	TextWindow checkWindow = new TextWindow(new String("DICOM info"),new String(""),800,400);
	checkWindow.append((String) imp.getProperty("Info"));
	*/
	
	public void run(ImageProcessor ip) {
		imageInfo = new Info().getImageInfo(imp,imp.getChannelProcessor());		
		sectorWidth = 10;
		cOn = true;
		dOn = true;
		resolution = 0;
		if (getInfoProperty(imageInfo,"Pixel Spacing")!= null){
			String temp = getInfoProperty(imageInfo,"Pixel Spacing");
			if (temp.indexOf("\\")!=-1){
				temp = temp.substring(0,temp.indexOf("\\"));
			}
			resolution = Double.valueOf(temp);
		}
		//Get parameters for scaling the image and for thresholding
		GenericDialog dialog = new GenericDialog("Analysis parameters");
		dialog.addNumericField("Fat threshold", 40.0, 4, 8, null);
		dialog.addNumericField("Area threshold", 550.0, 4, 8, null);
		dialog.addNumericField("BMD threshold", 690.0, 4, 8, null);
		if (getInfoProperty(imageInfo,"Rescale Slope")!= null){//Suggest HU scaling for dicom Files
			dialog.addNumericField("Scaling_coefficient (slope)", Double.valueOf(getInfoProperty(imageInfo,"Rescale Slope")), 4, 8, null);
			dialog.addNumericField("Scaling_constant (intercept)",Double.valueOf(getInfoProperty(imageInfo,"Rescale Intercept")), 4, 8, null);
					
		}else{
			dialog.addNumericField("Scaling_coefficient (slope)", 1.724, 4, 8, null);
			dialog.addNumericField("Scaling_constant (intercept)", -322.0, 4, 8, null);
		}
		if (resolution != 0){
			dialog.addNumericField("In-plane_pixel_size [mm]", resolution, 4, 8, null);
		} else {
			dialog.addNumericField("In-plane_pixel_size [mm]", 1.0, 4, 8, null);
		}
		//Get ROI selection
		String[] choiceLabels = {"Bigger","Smaller","Left","Right","Top","Bottom","Central","Peripheral"};
		dialog.addChoice("Roi_selection", choiceLabels, "Bigger"); 
		String[] rotationLabels = {"According_to_Imax/Imin","Furthest_point"};
		dialog.addChoice("Rotation_selection", rotationLabels, "According_to_Imax/Imin");
		dialog.addCheckbox("Analyse_cortical_results",true);
		dialog.addCheckbox("Analyse_density_distribution",true);
		dialog.addCheckbox("Allow_cleaving",false);
		dialog.addCheckbox("Cleave_retain_smaller",false);
		dialog.addCheckbox("Suppress_result_image",false);
		dialog.addCheckbox("Limit_ROI_search_to_manually_selected",false);
		dialog.addCheckbox("Set_distribution_results_rotation_manually",false);
		dialog.addNumericField("Manual_rotation_[+-_180_deg]", 0.0, 4, 8, null);
		dialog.addCheckbox("Flip_distribution_results",false);
		dialog.addCheckbox("Save_visual_result_image_on_disk",false);
		dialog.addStringField("Image_save_path",Prefs.getDefaultDirectory(),40);
		dialog.showDialog();
		
		if (dialog.wasOKed()){ //Stop in case of cancel..
			fatThreshold				= dialog.getNextNumber();
			areaThreshold				= dialog.getNextNumber();
			BMDThreshold				= dialog.getNextNumber();
			scalingFactor				= dialog.getNextNumber();
			constant					= dialog.getNextNumber();
			resolution					= dialog.getNextNumber();
			String roiChoice			= dialog.getNextChoice();
			String rotationChoice		= dialog.getNextChoice();
			cOn							= dialog.getNextBoolean();
			dOn							= dialog.getNextBoolean();
			boolean allowCleaving		= dialog.getNextBoolean();
			boolean cleaveReturnSmaller = dialog.getNextBoolean();
			boolean suppressImages		= dialog.getNextBoolean();
			boolean manualRoi			= dialog.getNextBoolean();
			manualRotation				= dialog.getNextBoolean();
			double manualAlfa			= dialog.getNextNumber()*Math.PI/180.0;
			flipDistribution			= dialog.getNextBoolean();
			boolean saveImageOnDisk		= dialog.getNextBoolean();
			String imageSavePath 		= dialog.getNextString();
			ScaledImageData scaledImageData;
			
			String imageName;
			if (getInfoProperty(imageInfo,"File Name")!= null){
				imageName = getInfoProperty(imageInfo,"File Name");
			}else{
				if(imp.getImageStackSize() == 1){
					imageName = getInfoProperty(imageInfo,"Title");
				}else{
					imageName = imageInfo.substring(0,imageInfo.indexOf("\n"));
				}
			}
			
			if (imp.getBitDepth() ==16){ //For unsigned short Dicom, which appears to be the default ImageJ DICOM...
				short[] tempPointer = (short[]) imp.getProcessor().getPixels();
				int[] unsignedShort = new int[tempPointer.length];
				for (int i=0;i<tempPointer.length;++i){unsignedShort[i] = 0x0000FFFF & (int) (tempPointer[i]);}
				scaledImageData = new ScaledImageData(unsignedShort, imp.getWidth(), imp.getHeight(),resolution, scalingFactor, constant,3);	//Scale and 3x3 median filter the data
			}else{		
				scaledImageData = new ScaledImageData((float[]) imp.getProcessor().getPixels(), imp.getWidth(), imp.getHeight(),resolution, scalingFactor, constant,3);	//Scale and 3x3 median filter the data
			}
			ImageAndAnalysisDetails imageAndAnalysisDetails = new ImageAndAnalysisDetails(scalingFactor, constant,fatThreshold, 
															areaThreshold,BMDThreshold,roiChoice,rotationChoice,choiceLabels,
															allowCleaving,cleaveReturnSmaller,manualRoi,manualRotation,manualAlfa,flipDistribution);
			SelectROI roi = new SelectROI(scaledImageData, imageAndAnalysisDetails,imp);
			TextWindow textWindow = (TextWindow) ij.WindowManager.getFrame("Distribution Analysis Results");
			if (textWindow == null){
				textWindow = new TextWindow(new String("Distribution Analysis Results"),new String(""),500,200);
				writeHeader(textWindow);
			}
			resultString = new String();
			printResults(textWindow);
			ImagePlus resultImage = null;
			if (cOn ){
				CorticalAnalysis cortAnalysis =new CorticalAnalysis(roi);
				printCorticalResults(textWindow,cortAnalysis);
				
				if(!dOn){
					BufferedImage bi = roi.getMyImage(roi.scaledImage,roi.sieve,roi.width,roi.height,roi.minimum,roi.maximum,dialog.getParent());
					resultImage = new ImagePlus("Visual results",bi);
				}
				
			}
			if (dOn){
				AnalyzeROI analyzeRoi = new AnalyzeROI(roi,imageAndAnalysisDetails);
				printDistributionResults(textWindow,analyzeRoi);
				BufferedImage bi = roi.getMyImage(roi.scaledImage,analyzeRoi.marrowCenter,analyzeRoi.pind,analyzeRoi.R,analyzeRoi.R2,analyzeRoi.Theta2,roi.width,roi.height,roi.minimum,roi.maximum,dialog.getParent()); // retrieve image
				resultImage = new ImagePlus("Visual results",bi);
			}
			
			if (!suppressImages && resultImage!= null){
				resultImage.show();
			}
			if (saveImageOnDisk && resultImage!= null){
				FileSaver fSaver = new FileSaver(resultImage);
				fSaver.saveAsPng(imageSavePath+"/"+imageName+".png"); 
			}
			textWindow.append(resultString);
		}
	}
	
	void writeHeader(TextWindow textWindow){
		String headerRow = new String();
		if (imp!=null){
			headerRow += "Filename\tSubject name\tSubject ID\tSubject birthdate\tMeasurement date\t";
		}
		headerRow += "Fat Threshold\tArea Threshold\tBMD Threshold\tScaling Coefficient\tScaling Constant\tResolution\t";
		if (cOn){
			headerRow +="CoD [mg/cm3]\tCoA [mm2]\tSSI [mm3]\tToD [mg/cm3]\tToA[mm2]\tBSId[g2/cm4]\t";
		}
		if (dOn){
			headerRow +="Alfa [deg]\tManual rotation [true/false]\tFlip [true/false]\t";
		
			for (int pp = 0;pp<((int) 360/sectorWidth);pp++){
				headerRow +="Endocortical radius "+pp*sectorWidth+" - "+((pp+1)*sectorWidth)+" [mm]\t";
			}
			for (int pp = 0;pp<((int) 360/sectorWidth);pp++){
				headerRow +="Pericortical radius "+pp*sectorWidth+" - "+((pp+1)*sectorWidth)+" [mm]\t";
			}
			//Cortex BMD values			
			for (int pp = 0;pp<((int) 360/sectorWidth);pp++){
				headerRow +="Endocortical BMD "+pp*sectorWidth+" - "+((pp+1)*sectorWidth)+" [mg/cm3]\t";
			}
			for (int pp = 0;pp<((int) 360/sectorWidth);pp++){
				headerRow +="Midcortical BMD "+pp*sectorWidth+" - "+((pp+1)*sectorWidth)+" [mg/cm3]\t";
			}
			for (int pp = 0;pp<((int) 360/sectorWidth);pp++){
				headerRow +="Pericortical BMD "+pp*sectorWidth+" - "+((pp+1)*sectorWidth)+" [mg/cm3]\t";
			}
		}
		textWindow.append(headerRow);
	}
	
	void printResults(TextWindow textWindow){
		if (imp != null){
			if (getInfoProperty(imageInfo,"File Name")!= null){
				resultString += getInfoProperty(imageInfo,"File Name")+"\t";
			}else{
				if(imp.getImageStackSize() == 1){
					resultString += getInfoProperty(imageInfo,"Title")+"\t";
				}else{
					resultString += imageInfo.substring(0,imageInfo.indexOf("\n"))+"\t";
				}
			}
			resultString += getInfoProperty(imageInfo,"Patient's Name")+"\t";
			resultString += getInfoProperty(imageInfo,"Patient ID")+"\t";
			resultString += getInfoProperty(imageInfo,"Patient's Birth Date")+"\t";
			resultString += getInfoProperty(imageInfo,"Acquisition Date")+"\t";
		}
		
		resultString += Double.toString(fatThreshold)+"\t";
		resultString += Double.toString(areaThreshold)+"\t";
		resultString += Double.toString(BMDThreshold)+"\t";
		resultString += Double.toString(scalingFactor)+"\t";
		resultString += Double.toString(constant)+"\t";	
		resultString += Double.toString(resolution)+"\t";
	}
	
		String getInfoProperty(String properties,String propertyToGet){
			String toTokenize = properties;
			StringTokenizer st = new StringTokenizer(toTokenize,"\n");
			String currentToken = null;
			while (st.hasMoreTokens() ) {
				currentToken = st.nextToken();
				if (currentToken.indexOf(propertyToGet) != -1){break;}
			}
			if (currentToken.indexOf(propertyToGet) != -1){
				StringTokenizer st2 = new StringTokenizer(currentToken,":");
				String token2 = null;
				while (st2.hasMoreTokens()){
					token2 = st2.nextToken();
				}
				return token2.trim();
			}
			return null;
	}
	
	void printCorticalResults(TextWindow textWindow,CorticalAnalysis cortAnalysis){
		resultString +=Double.toString(cortAnalysis.BMD)+"\t";
		resultString +=Double.toString(cortAnalysis.AREA)+"\t";
		resultString +=Double.toString(cortAnalysis.SSI)+"\t";
		resultString +=Double.toString(cortAnalysis.ToD)+"\t";
		resultString +=Double.toString(cortAnalysis.ToA)+"\t";
		resultString +=Double.toString(cortAnalysis.BSId)+"\t";
	}
	
	void printDistributionResults(TextWindow textWindow,AnalyzeROI analyzeRoi){
		resultString += Double.toString(analyzeRoi.alfa*180/Math.PI)+"\t";
		resultString += Boolean.toString(manualRotation)+"\t";
		resultString += Boolean.toString(flipDistribution)+"\t";
		
		for (int pp = 0;pp<((int) 360/sectorWidth);pp++){
			resultString += analyzeRoi.endocorticalRadii[pp]+"\t";
		}
		for (int pp = 0;pp<((int) 360/sectorWidth);pp++){
			resultString += analyzeRoi.pericorticalRadii[pp]+"\t";
		}
		//Cortex BMD values			
		for (int pp = 0;pp<((int) 360/sectorWidth);pp++){
			resultString += analyzeRoi.endoCorticalBMDs[pp]+"\t";
		}
		for (int pp = 0;pp<((int) 360/sectorWidth);pp++){
			resultString += analyzeRoi.midCorticalBMDs[pp]+"\t";
		}
		for (int pp = 0;pp<((int) 360/sectorWidth);pp++){
			resultString += analyzeRoi.periCorticalBMDs[pp]+"\t";
		}
	}
}
