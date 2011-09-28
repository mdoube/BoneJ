package org.doube.bonej.pqct;

import ij.*;
import ij.text.*;
import ij.process.*;
import ij.gui.*;
import java.awt.*;
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
	double areaThreshold;
	double BMDThreshold;
	double scalingFactor;
	double constant;
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
		double resolution;
		if (imp.getProperty(new String("VoxelSize")) != null){
			resolution = (Double) imp.getProperty(new String("VoxelSize"));
		}else{
			resolution =0;
		}
		//Get parameters for scaling the image and for thresholding
		GenericDialog dialog = new GenericDialog(new String("Analysis parameters"));
		dialog.addNumericField(new String("Area threshold"), 550.0, 4);
		dialog.addNumericField(new String("BMD threshold"), 690.0, 4);
		dialog.addNumericField(new String("Scaling coefficient (slope)"), 1.724, 4);
		dialog.addNumericField(new String("Scaling constant (intercept)"), -322.0, 4);
		if (resolution != 0){
			dialog.addNumericField(new String("In-plane pixel size [mm]"), resolution, 4);
		} else {
			dialog.addNumericField(new String("In-plane pixel size [mm]"), 0.8, 4);
		}
			
		dialog.showDialog();
		areaThreshold	= dialog.getNextNumber();
		BMDThreshold	= dialog.getNextNumber();
		scalingFactor	= dialog.getNextNumber();
		constant		= dialog.getNextNumber();
		resolution		= dialog.getNextNumber();
		//System.out.println(areaThreshold  +" "+ BMDThreshold +" "+ scalingFactor +" "+ constant);
		//ImagePlus dataIn = new ImagePlus("File in", ip);		//Does not give image properties
		//System.out.println(imp.getProperties());
		ScaledImageData scaledImageData;
		if (imp.getBitDepth() ==16){ //For unsigned short Dicom, which appears to be the default ImageJ DICOM...
			short[] tempPointer = (short[]) imp.getProcessor().getPixels();
			int[] unsignedShort = new int[tempPointer.length];
			for (int i=0;i<tempPointer.length;++i){unsignedShort[i] = 0x0000FFFF & (int) (tempPointer[i]);}
			scaledImageData = new ScaledImageData(unsignedShort, imp.getWidth(), imp.getHeight(),resolution, scalingFactor, constant,3);	//Scale and 3x3 median filter the data
		}else{		
			scaledImageData = new ScaledImageData((float[]) imp.getProcessor().getPixels(), imp.getWidth(), imp.getHeight(),resolution, scalingFactor, constant,3);	//Scale and 3x3 median filter the data
		}
		ImageAndAnalysisDetails imageAndAnalysisDetails = new ImageAndAnalysisDetails(scalingFactor, constant, areaThreshold,BMDThreshold);
		SelectROI roi = new SelectROI(scaledImageData, imageAndAnalysisDetails,imp);
		CorticalAnalysis cortAnalysis =new CorticalAnalysis(roi);
		AnalyzeROI analyzeRoi = new AnalyzeROI(roi,imageAndAnalysisDetails);
		
		String resultString = new String();
		TextWindow textWindow = new TextWindow(new String("Results"),new String(""),500,200);
		printResults(textWindow,cortAnalysis,analyzeRoi);
		
		/*Display the analysis results...*/
		BufferedImage bi = roi.getMyImage(roi.scaledImage,analyzeRoi.marrowCenter,analyzeRoi.pind,analyzeRoi.R,analyzeRoi.R2,analyzeRoi.Theta2,roi.width,roi.height,roi.minimum,roi.maximum,dialog.getParent()); // retrieve image
		ImagePlus resultImage = new ImagePlus("Visual results",bi);
		resultImage.show();
	}
	
	void printResults(TextWindow textWindow,CorticalAnalysis cortAnalysis,AnalyzeROI analyzeRoi){
		if (imp.getProperty(new String("FileName")) != null){
			textWindow.append("Filename\t"+(String) imp.getProperty(new String("FileName")));
			textWindow.append("Subject name\t"+(String) imp.getProperty(new String("PatName")));
			textWindow.append("Subject ID\t"+(String) imp.getProperty(new String("PatID")));
			textWindow.append("Subject birthdate\t"+(Long) imp.getProperty(new String("PatBirth")));
			textWindow.append("Measurement date\t"+(Long) imp.getProperty(new String("MeasDate")));
		}
		if (cOn){
			textWindow.append("BMD\t"+Double.toString(cortAnalysis.BMD)+"\tmg/cm3");
			textWindow.append("AREA\t"+Double.toString(cortAnalysis.AREA)+"\tmm2");
			textWindow.append("SSI\t"+Double.toString(cortAnalysis.SSI)+"\tmm3");
		}
		if (dOn){
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
}
