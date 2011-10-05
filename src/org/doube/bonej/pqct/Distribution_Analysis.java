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
	public int setup(String arg, ImagePlus imp) {
		this.imp = imp;
		//return DOES_32;
		return DOES_ALL;
	}

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
		if (imp.getOriginalFileInfo().fileFormat == ij.io.FileInfo.DICOM ){/*Suggest HU scaling for dicom Files*/
			double[] coeffs = imp.getCalibration().getCoefficients();
			dialog.addNumericField(new String("Scaling coefficient (slope)"), coeffs[1], 4);
			dialog.addNumericField(new String("Scaling constant (intercept)"),coeffs[0], 4);
			/*
			//For Debugging
			TextWindow checkWindow = new TextWindow(new String("test"),new String(""),200,200);			
			for (int i = 0; i<coeffs.length;++i){
				checkWindow.append("Coefficients "+i+" "+coeffs[i]);
			}
			*/			
		}else{
			dialog.addNumericField(new String("Scaling coefficient (slope)"), 1.724, 4);
			dialog.addNumericField(new String("Scaling constant (intercept)"), -322.0, 4);
		}
		if (resolution != 0){
			dialog.addNumericField(new String("In-plane pixel size [mm]"), resolution, 4);
		} else {
			dialog.addNumericField(new String("In-plane pixel size [mm]"), 0.8, 4);
		}
		//Get ROI selection
		String[] choiceLabels = {"Bigger","Smaller","Left","Right","Top","Bottom","Central","Peripheral"};
		dialog.addChoice("Roi selection", choiceLabels, "Bigger"); 
		String[] rotationLabels = {"According to Imax/Imin","Furthest point"};
		dialog.addChoice("Rotation selection", rotationLabels, "According to Imax/Imin");
		dialog.addCheckbox("Analyse cortical results",true);
		dialog.addCheckbox("Analyse density distribution",true);
		dialog.addCheckbox("Allow cleaving",false);
		dialog.addMessage("N.B. If a ROI has been manually defined,"); 
		dialog.addMessage("the Automated selection will still be used."); 
		dialog.addMessage("Manual ROI only allows preventing"); 
		dialog.addMessage("parts of the image from being considered in"); 
		dialog.addMessage("automatically detecting the bone."); 
		dialog.showDialog();
		
		if (dialog.wasOKed()){ //Stop in case of cancel..
			double fatThreshold		= dialog.getNextNumber();
			double areaThreshold	= dialog.getNextNumber();
			double BMDThreshold		= dialog.getNextNumber();
			double scalingFactor	= dialog.getNextNumber();
			double constant			= dialog.getNextNumber();
			resolution				= dialog.getNextNumber();
			String roiChoice		= dialog.getNextChoice();
			String rotationChoice	= dialog.getNextChoice();
			cOn						= dialog.getNextBoolean();
			dOn						= dialog.getNextBoolean();
			boolean allowCleaving	= dialog.getNextBoolean();
			/*
			//For debugging
			TextWindow checkWindow = new TextWindow(new String("Checkboxes"),new String(""),500,200);
			checkWindow.append("Selection "+roiChoice);
			*/
			ScaledImageData scaledImageData;
			if (imp.getBitDepth() ==16){ //For unsigned short Dicom, which appears to be the default ImageJ DICOM...
				short[] tempPointer = (short[]) imp.getProcessor().getPixels();
				int[] unsignedShort = new int[tempPointer.length];
				for (int i=0;i<tempPointer.length;++i){unsignedShort[i] = 0x0000FFFF & (int) (tempPointer[i]);}
				scaledImageData = new ScaledImageData(unsignedShort, imp.getWidth(), imp.getHeight(),resolution, scalingFactor, constant,3);	//Scale and 3x3 median filter the data
			}else{		
				scaledImageData = new ScaledImageData((float[]) imp.getProcessor().getPixels(), imp.getWidth(), imp.getHeight(),resolution, scalingFactor, constant,3);	//Scale and 3x3 median filter the data
			}
			ImageAndAnalysisDetails imageAndAnalysisDetails = new ImageAndAnalysisDetails(scalingFactor, constant,fatThreshold, areaThreshold,BMDThreshold,roiChoice,rotationChoice,choiceLabels,allowCleaving);
			SelectROI roi = new SelectROI(scaledImageData, imageAndAnalysisDetails,imp);
			TextWindow textWindow = new TextWindow(new String("Results"),new String(""),500,200);
			printResults(textWindow);
			if (cOn){
				CorticalAnalysis cortAnalysis =new CorticalAnalysis(roi);
				printCorticalResults(textWindow,cortAnalysis);
				/*
				if(!dOn){
					BufferedImage bi = roi.getMyImage(roi.scaledImage,analyzeRoi.marrowCenter,analyzeRoi.pind,analyzeRoi.R,analyzeRoi.R2,analyzeRoi.Theta2,roi.width,roi.height,roi.minimum,roi.maximum,dialog.getParent()); // retrieve image
					ImagePlus resultImage = new ImagePlus("Visual results",bi);
					resultImage.show();
				}
				*/
			}
			if (dOn){
				AnalyzeROI analyzeRoi = new AnalyzeROI(roi,imageAndAnalysisDetails);
				printDistributionResults(textWindow,analyzeRoi);
				BufferedImage bi = roi.getMyImage(roi.scaledImage,analyzeRoi.marrowCenter,analyzeRoi.pind,analyzeRoi.R,analyzeRoi.R2,analyzeRoi.Theta2,roi.width,roi.height,roi.minimum,roi.maximum,dialog.getParent()); // retrieve image
				ImagePlus resultImage = new ImagePlus("Visual results",bi);
				resultImage.show();
			}
			String resultString = new String();

			
			//Display the analysis results...
			
		}
	}
	
	void printResults(TextWindow textWindow){
		if (imp != null){
			textWindow.append("Filename\t"+(String) imp.getProperty(new String("FileName")));
			textWindow.append("Subject name\t"+(String) imp.getProperty(new String("PatName")));
			textWindow.append("Subject ID\t"+(String) imp.getProperty(new String("PatID")));
			textWindow.append("Subject birthdate\t"+(Long) imp.getProperty(new String("PatBirth")));
			textWindow.append("Measurement date\t"+(Long) imp.getProperty(new String("MeasDate")));
		}
	}
	void printCorticalResults(TextWindow textWindow,CorticalAnalysis cortAnalysis){
		textWindow.append("CoD\t"+Double.toString(cortAnalysis.BMD)+"\tmg/cm3");
		textWindow.append("CoA\t"+Double.toString(cortAnalysis.AREA)+"\tmm2");
		textWindow.append("SSI\t"+Double.toString(cortAnalysis.SSI)+"\tmm3");
		textWindow.append("ToD\t"+Double.toString(cortAnalysis.ToD)+"\tmg/cm3");
		textWindow.append("ToA\t"+Double.toString(cortAnalysis.ToA)+"\tmm2");
		textWindow.append("BSId\t"+Double.toString(cortAnalysis.BSId)+"\tg2/cm4");
	}
	void printDistributionResults(TextWindow textWindow,AnalyzeROI analyzeRoi){

		for (int pp = 0;pp<((int) 360/sectorWidth);pp++){
			textWindow.append("Endocortical radius "+pp*sectorWidth+" - "+((pp+1)*sectorWidth)+"\t"+analyzeRoi.endocorticalRadii[pp]+"\tmm");
		}
		for (int pp = 0;pp<((int) 360/sectorWidth);pp++){
			textWindow.append("Pericortical radius "+pp*sectorWidth+" - "+((pp+1)*sectorWidth)+"\t"+analyzeRoi.pericorticalRadii[pp]+"\tmm");
		}
		/*	//Radii after peeling off one layer of pixels
		for (int pp = 0;pp<((int) 360/sectorWidth);pp++){
			textWindow.append(analyzeRoi.peeledEndocorticalRadii[pp]+" mg/cm3");
		}
		for (int pp = 0;pp<((int) 360/sectorWidth);pp++){
			textWindow.append(analyzeRoi.peeledPericorticalRadii[pp]+" mg/cm3");
		}
		*/
		//Cortex BMD values			
		for (int pp = 0;pp<((int) 360/sectorWidth);pp++){
			textWindow.append("Endocortical BMD "+pp*sectorWidth+" - "+((pp+1)*sectorWidth)+"\t"+analyzeRoi.endoCorticalBMDs[pp]+"\tmg/cm3");
		}
		for (int pp = 0;pp<((int) 360/sectorWidth);pp++){
			textWindow.append("Midcortical BMD "+pp*sectorWidth+" - "+((pp+1)*sectorWidth)+"\t"+analyzeRoi.midCorticalBMDs[pp]+"\tmg/cm3");
		}
		for (int pp = 0;pp<((int) 360/sectorWidth);pp++){
			textWindow.append("Pericortical BMD "+pp*sectorWidth+" - "+((pp+1)*sectorWidth)+"\t"+analyzeRoi.periCorticalBMDs[pp]+"\tmg/cm3");
		}
	}

}
