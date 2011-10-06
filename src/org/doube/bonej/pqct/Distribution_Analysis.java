package org.doube.bonej.pqct;

import ij.*;

import ij.text.*;
import ij.process.*;
import ij.gui.*;
import java.awt.*;
import java.util.*;	//Vector
import ij.plugin.filter.*;
import org.doube.bonej.pqct.analysis.*;		//Analysis stuff..
import org.doube.bonej.pqct.selectroi.*;		//ROI selection..
import org.doube.bonej.pqct.io.*;	//image data 
import java.awt.image.*; //Creating the result BufferedImage...

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
		sectorWidth = 10;
		cOn = true;
		dOn = true;
		double resolution = 0;
		if (imp.getOriginalFileInfo().pixelWidth != 0){
			resolution = imp.getOriginalFileInfo().pixelWidth;
		}
	
		//Get parameters for scaling the image and for thresholding
		GenericDialog dialog = new GenericDialog(new String("Analysis parameters"));
		dialog.addNumericField(new String("Fat threshold"), 40.0, 4);
		dialog.addNumericField(new String("Area threshold"), 550.0, 4);
		dialog.addNumericField(new String("BMD threshold"), 690.0, 4);
		if (imp.getOriginalFileInfo().fileFormat == ij.io.FileInfo.DICOM ){//Suggest HU scaling for dicom Files
			double[] coeffs = imp.getCalibration().getCoefficients();
			dialog.addNumericField(new String("Scaling_coefficient (slope)"), coeffs[1], 4);
			dialog.addNumericField(new String("Scaling_constant (intercept)"),coeffs[0], 4);
					
		}else{
			dialog.addNumericField(new String("Scaling_coefficient (slope)"), 1.724, 4);
			dialog.addNumericField(new String("Scaling_constant (intercept)"), -322.0, 4);
		}
		if (resolution != 0){
			dialog.addNumericField(new String("In-plane_pixel_size [mm]"), resolution, 4);
		} else {
			dialog.addNumericField(new String("In-plane_pixel_size [mm]"), 0.8, 4);
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
		dialog.showDialog();
		
		if (dialog.wasOKed()){ //Stop in case of cancel..
			double fatThreshold			= dialog.getNextNumber();
			double areaThreshold		= dialog.getNextNumber();
			double BMDThreshold			= dialog.getNextNumber();
			double scalingFactor		= dialog.getNextNumber();
			double constant				= dialog.getNextNumber();
			resolution					= dialog.getNextNumber();
			String roiChoice			= dialog.getNextChoice();
			String rotationChoice		= dialog.getNextChoice();
			cOn							= dialog.getNextBoolean();
			dOn							= dialog.getNextBoolean();
			boolean allowCleaving		= dialog.getNextBoolean();
			boolean cleaveReturnSmaller = dialog.getNextBoolean();
			boolean suppressImages		= dialog.getNextBoolean();
			ScaledImageData scaledImageData;
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
															allowCleaving,cleaveReturnSmaller);
			SelectROI roi = new SelectROI(scaledImageData, imageAndAnalysisDetails,imp);
			TextWindow textWindow = (TextWindow) ij.WindowManager.getFrame("Distribution Analysis Results");
			if (textWindow == null){
				textWindow = new TextWindow(new String("Distribution Analysis Results"),new String(""),500,200);
				writeHeader(textWindow);
			}
			resultString = new String();
			printResults(textWindow);
			if (cOn ){
				CorticalAnalysis cortAnalysis =new CorticalAnalysis(roi);
				printCorticalResults(textWindow,cortAnalysis);
				
				if(!dOn && !suppressImages){
					BufferedImage bi = roi.getMyImage(roi.scaledImage,roi.sieve,roi.width,roi.height,roi.minimum,roi.maximum,dialog.getParent());
					ImagePlus resultImage = new ImagePlus("Visual results",bi);
					resultImage.show();
				}
				
			}
			if (dOn){
				AnalyzeROI analyzeRoi = new AnalyzeROI(roi,imageAndAnalysisDetails);
				printDistributionResults(textWindow,analyzeRoi);
				if (!suppressImages){
					BufferedImage bi = roi.getMyImage(roi.scaledImage,analyzeRoi.marrowCenter,analyzeRoi.pind,analyzeRoi.R,analyzeRoi.R2,analyzeRoi.Theta2,roi.width,roi.height,roi.minimum,roi.maximum,dialog.getParent()); // retrieve image
					ImagePlus resultImage = new ImagePlus("Visual results",bi);
					resultImage.show();
				}
			}
			textWindow.append(resultString);
			//Display the analysis results...
		}
	}
	
	void writeHeader(TextWindow textWindow){
		String headerRow = new String();
		headerRow = "Filename\tSubject name\tSubject ID\tSubject birthdate\tMeasurement date\t";
		if (cOn){
			headerRow +="CoD [mg/cm3]\tCoA [mm2]\tSSI [mm3]\tToD [mg/cm3]\tToA[mm2]\tBSId[g2/cm4]\t";
		}
		if (dOn){
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
			resultString += imp.getOriginalFileInfo().fileName+"\t";
			resultString += getInfoProperty((String) imp.getProperty("Info"),"Patient's Name")+"\t";
			resultString += getInfoProperty((String) imp.getProperty("Info"),"Patient ID")+"\t";
			resultString += getInfoProperty((String) imp.getProperty("Info"),"Patient's Birth Date")+"\t";
			resultString += getInfoProperty((String) imp.getProperty("Info"),"Acquisition Date")+"\t";
		}
	}
	
		String getInfoProperty(String properties,String propertyToGet){
			String toTokenize = (String) imp.getProperty("Info");
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
