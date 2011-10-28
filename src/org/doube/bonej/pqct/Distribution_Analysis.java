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

package org.doube.bonej.pqct;

import ij.*;
import ij.text.*;
import ij.process.*;
import ij.gui.*;
import ij.measure.*;						//Calibration
import java.util.*;							//Vector

import ij.plugin.PlugIn;
import org.doube.bonej.pqct.analysis.*;		//Analysis stuff..
import org.doube.bonej.pqct.selectroi.*;	//ROI selection..
import org.doube.bonej.pqct.io.*;			//image data 
import java.awt.*;							//Image, component for debugging...
import ij.plugin.filter.Info;
import ij.io.*;

public class Distribution_Analysis implements PlugIn {

	int sectorWidth;
	boolean cOn;
	boolean	mOn;
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
	boolean guessFlip;
	boolean stacked;
	boolean manualRotation;
	boolean allowCleaving;
	String roiChoice;
	String rotationChoice;
	
	
	public void run(String arg) {
		ImagePlus imp = WindowManager.getCurrentImage();
		if (imp == null)
			return;
		if (imp.getType() != ImagePlus.GRAY16){
			IJ.error("Distribution analysis expects 16-bit greyscale data");
			return;
		}
		imageInfo = new Info().getImageInfo(imp,imp.getChannelProcessor());
		/*Check image calibration*/
		Calibration cal = imp.getCalibration();
		double[] calibrationCoefficients = {0,1};
		if (getInfoProperty(imageInfo,"Stratec File") == null){ 
				if (cal.getCoefficients() != null)
			calibrationCoefficients = cal.getCoefficients();
		} else {
			calibrationCoefficients = new double[2];
			calibrationCoefficients[0] = -322.0;
			calibrationCoefficients[1] = 1.724;
		}
		sectorWidth = 10;
		cOn = true;
		dOn = true;
		mOn = true;
		resolution = cal.pixelWidth;
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
		dialog.addNumericField("Scaling_coefficient (slope)", calibrationCoefficients[1], 4, 8, null);
		dialog.addNumericField("Scaling_constant (intercept)",calibrationCoefficients[0], 4, 8, null);
		//Get ROI selection
		String[] choiceLabels = {"Bigger","Smaller","Left","Right","Top","Bottom","Central","Peripheral"};
		dialog.addChoice("Roi_selection", choiceLabels, "Bigger"); 
		String[] rotationLabels = {"According_to_Imax/Imin","Furthest_point","All_Bones_Imax/Imin"};
		dialog.addChoice("Rotation_selection", rotationLabels, "According_to_Imax/Imin");
		dialog.addCheckbox("Analyse_cortical_results",true);
		dialog.addCheckbox("Analyse_mass_distribution",true);
		dialog.addCheckbox("Analyse_density_distribution",true);
		dialog.addCheckbox("Allow_cleaving",false);
		dialog.addCheckbox("Suppress_result_image",false);
		dialog.addCheckbox("Limit_ROI_search_to_manually_selected",false);
		dialog.addCheckbox("Set_distribution_results_rotation_manually",false);
		dialog.addNumericField("Manual_rotation_[+-_180_deg]", 0.0, 4, 8, null);
		dialog.addCheckbox("Flip_distribution_results",false);
		dialog.addCheckbox("Guess_right",false);
		dialog.addCheckbox("Stacked_bones",true);
		dialog.addCheckbox("Save_visual_result_image_on_disk",false);
		dialog.addStringField("Image_save_path",Prefs.getDefaultDirectory(),40);
		dialog.showDialog();
		
		if (dialog.wasOKed()){ //Stop in case of cancel..
			fatThreshold				= dialog.getNextNumber();
			areaThreshold				= dialog.getNextNumber();
			BMDThreshold				= dialog.getNextNumber();
			scalingFactor				= dialog.getNextNumber();
			constant					= dialog.getNextNumber();
			roiChoice			= dialog.getNextChoice();
			rotationChoice		= dialog.getNextChoice();
			cOn							= dialog.getNextBoolean();
			mOn							= dialog.getNextBoolean();
			dOn							= dialog.getNextBoolean();
			allowCleaving		= dialog.getNextBoolean();
			boolean suppressImages		= dialog.getNextBoolean();
			boolean manualRoi			= dialog.getNextBoolean();
			manualRotation				= dialog.getNextBoolean();
			double manualAlfa			= dialog.getNextNumber()*Math.PI/180.0;
			flipDistribution			= dialog.getNextBoolean();
			guessFlip					= dialog.getNextBoolean();
			stacked						= dialog.getNextBoolean();
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

			short[] tempPointer = (short[]) imp.getProcessor().getPixels();			
			int[] unsignedShort = new int[tempPointer.length];			
			if (getInfoProperty(imageInfo,"Stratec File") == null){ //For unsigned short Dicom, which appears to be the default ImageJ DICOM...
				for (int i=0;i<tempPointer.length;++i){unsignedShort[i] = 0x0000FFFF & (int) (tempPointer[i]);}
			}else{
				float[] floatPointer = (float[]) imp.getProcessor().toFloat(1,null).getPixels();
				for (int i=0;i<tempPointer.length;++i){unsignedShort[i] = (int) (floatPointer[i] - Math.pow(2.0,15.0));}
			}
			scaledImageData = new ScaledImageData(unsignedShort, imp.getWidth(), imp.getHeight(),resolution, scalingFactor, constant,3);	//Scale and 3x3 median filter the data
			
			ImageAndAnalysisDetails imageAndAnalysisDetails = new ImageAndAnalysisDetails(scalingFactor, constant,fatThreshold, 
															areaThreshold,BMDThreshold,roiChoice,rotationChoice,choiceLabels,
															allowCleaving,manualRoi,manualRotation,manualAlfa,flipDistribution,
															guessFlip, stacked);
			SelectROI roi = new SelectROI(scaledImageData, imageAndAnalysisDetails,imp);
			DetermineAlfa determineAlfa = new DetermineAlfa(roi,imageAndAnalysisDetails);
			imageAndAnalysisDetails.flipDistribution = roi.details.flipDistribution;
			flipDistribution = imageAndAnalysisDetails.flipDistribution;
			TextPanel textPanel = IJ.getTextPanel();
			if (textPanel == null) {textPanel = new TextPanel();}
			if (textPanel.getLineCount() == 0){writeHeader(textPanel);}
			
			String results = "";
			results = printResults(results,determineAlfa, imp);
			ImagePlus resultImage = null;
			if (cOn ){
				CorticalAnalysis cortAnalysis =new CorticalAnalysis(roi);
				results = printCorticalResults(results,cortAnalysis);
				if(!dOn && !mOn){
					resultImage = getResultImage(roi.scaledImage,roi.width,roi.height);
				}
				
			}
			if (mOn){
				MassDistribution massDistribution =new MassDistribution(roi,imageAndAnalysisDetails,determineAlfa);
				results = printMassDistributionResults(results,massDistribution);
				if(!dOn){
					resultImage = getResultImage(roi.scaledImage,roi.width,roi.height,roi.minimum,roi.maximum,roi.sieve,determineAlfa.alfa/Math.PI*180.0);

				}
			}
			
			if (dOn){
				DistributionAnalysis DistributionAnalysis = new DistributionAnalysis(roi,imageAndAnalysisDetails,determineAlfa);
				results = printDistributionResults(results,DistributionAnalysis);
				resultImage = getResultImage(roi.scaledImage,roi.width,roi.height,roi.minimum,roi.maximum,roi.sieve,determineAlfa.alfa/Math.PI*180.0,DistributionAnalysis.marrowCenter,determineAlfa.pind, determineAlfa.pindColor,DistributionAnalysis.R,DistributionAnalysis.R2,DistributionAnalysis.Theta);
			}
			
			if (!suppressImages && resultImage!= null){
				resultImage.show();
			}
			if (saveImageOnDisk && resultImage!= null){
				FileSaver fSaver = new FileSaver(resultImage);
				fSaver.saveAsPng(imageSavePath+"/"+imageName+".png"); 
			}
			textPanel.appendLine(results);
			textPanel.updateDisplay();			
		}
	}

	/*Cortical analysis result image*/
	ImagePlus getResultImage(double[] values,int width,int height){
		ImagePlus tempImage = new ImagePlus("Visual results");
		tempImage.setProcessor(new FloatProcessor(width,height,values));
		new ImageConverter(tempImage).convertToRGB();
		tempImage.setProcessor(tempImage.getProcessor().resize(1000));
		return tempImage;
	}
	
	/*Mass distribution result image*/
	ImagePlus getResultImage(double[] values,int width,int height, double min, double max, byte[] sieve, double alfa){
		ImagePlus tempImage = new ImagePlus("Visual results");
		tempImage.setProcessor(new FloatProcessor(width,height,values));
		new ImageConverter(tempImage).convertToRGB();
		for (int y = 0; y < height;++y) {
			for (int x = 0; x < width;++x) {
				if (sieve[x+y*width] == 1){   //Tint roi area color with violet
					double scale = (values[x+y*tempImage.getWidth()]-min)/(max-min);
					tempImage.getProcessor().setColor(new Color((int) (127.0*scale),0,(int) (255.0*scale)));
					tempImage.getProcessor().drawPixel(x,y);
				}
			}
		}
		tempImage.getProcessor().setInterpolate(true);
		tempImage.getProcessor().rotate(alfa);
		tempImage.setProcessor(tempImage.getProcessor().resize(1000));
		return tempImage;
	}
	
	/*Density distribution result image*/
	ImagePlus getResultImage(double[] values,int width,int height, double min, double max, byte[] sieve, double alfa,
							double[] marrowCenter,Vector<Integer> pind,Vector<Integer> pindColor, double[] R, double[] R2, double[] Theta){
		ImagePlus tempImage = new ImagePlus("Visual results");
		tempImage.setProcessor(new FloatProcessor(width,height,values));
		new ImageConverter(tempImage).convertToRGB();
		for (int y = 0; y < height;++y) {
			for (int x = 0; x < width;++x) {
				if (sieve[x+y*width] == 1){   //Tint roi area color with violet
					double scale = (values[x+y*tempImage.getWidth()]-min)/(max-min);
					tempImage.getProcessor().setColor(new Color((int) (127.0*scale),0,(int) (255.0*scale)));
					tempImage.getProcessor().drawPixel(x,y);
				}
			}
		}
		 
		//Draw unrotated radii
		for(int i = 0; i< 360;i++) {//45;i++) {//
			int x = ((int) (marrowCenter[0]+R[i]*Math.cos(Theta[i])));
			int y = ((int) (marrowCenter[1]+R[i]*Math.sin(Theta[i])));
			double colorScale = ((double) pindColor.get(i))/359.0;
			tempImage.getProcessor().setColor(new Color((int) (255.0*colorScale),0,(int) (255.0*(1.0-colorScale))));
			tempImage.getProcessor().drawPixel(x,y);
			x = ((int) (marrowCenter[0]+R2[i]*Math.cos(Theta[i])));
			y = ((int) (marrowCenter[1]+R2[i]*Math.sin(Theta[i])));
			tempImage.getProcessor().setColor(new Color(0,(int) (255.0*colorScale),(int) (255.0*(1.0-colorScale))));
			tempImage.getProcessor().drawPixel(x,y);
		}
		/*plot marrowCenter*/
		for(int i = 0; i< 10;i++) {//45;i++) {//
			int x = ((int) (marrowCenter[0]+i));
			int y = ((int) (marrowCenter[1]));
			tempImage.getProcessor().setColor(new Color(0,255,255));
			tempImage.getProcessor().drawPixel(x,y);
			x = (int) (marrowCenter[0]);
			y = (int) (marrowCenter[1]+i);
			tempImage.getProcessor().setColor(new Color(255,0,255));
			tempImage.getProcessor().drawPixel(x,y);
			/*Plot rotated axes...*/
			x = ((int) ((double) marrowCenter[0]+((double) i)*Math.cos(-alfa/180*Math.PI)));
			y = ((int) ((double) marrowCenter[1]+((double) i)*Math.sin(-alfa/180*Math.PI)));
			tempImage.getProcessor().setColor(new Color(0,255,0));
			tempImage.getProcessor().drawPixel(x,y);
			x = ((int) ((double) marrowCenter[0]+((double) -i)*Math.sin(-alfa/180*Math.PI)));
			y = ((int) ((double) marrowCenter[1]+((double) i)*Math.cos(-alfa/180*Math.PI)));
			tempImage.getProcessor().setColor(new Color(0,0,255));
			tempImage.getProcessor().drawPixel(x,y);

		}
		tempImage.getProcessor().setInterpolate(true);
		tempImage.getProcessor().rotate(alfa);
		tempImage.setProcessor(tempImage.getProcessor().resize(1000));
		return tempImage;
	}
	
	void writeHeader(TextPanel textPanel){
		String[] propertyNames = {"File Name","Patient's Name","Patient ID","Patient's Birth Date","Acquisition Date","Pixel Spacing"};
		String[] parameterNames = {"Fat Threshold","Area Threshold","BMD Threshold","Scaling Coefficient","Scaling Constant"};
		String[] dHeadings = {"Alfa [deg]","Rotation correction [deg]","Manual Rotation","Flip Distribution","Guess right"
		,"Stacked bones","Allow Cleaving","Roi choice","Rotation choice"};
			
		String headings = "";
		for (int i = 0;i<propertyNames.length;++i){
			headings+=propertyNames[i]+"\t";
		}
		for (int i = 0;i<parameterNames.length;++i){
			headings+=parameterNames[i]+"\t";
		}
		for (int i = 0;i<dHeadings.length;++i){
				headings+=dHeadings[i]+"\t";
		}
		
		if(cOn){
			String[] coHeadings = {"CoD [mg/cm³]","CoA [mm²]","SSI [mm³]","ToD [mg/cm³]","ToA[mm²]","BSId[g²/cm4]"};
			for (int i = 0;i<coHeadings.length;++i){
				headings+=coHeadings[i]+"\t";
			}
		}
		if(mOn){
			for (int i = 0;i<((int) 360/sectorWidth);++i){
				headings+=i*sectorWidth+"° - "+((i+1)*sectorWidth)+"° mineral mass [mg]\t";
			}
		}
		
		if(dOn){
			for (int i = 0;i<((int) 360/sectorWidth);++i){
				headings+=i*sectorWidth+"° - "+((i+1)*sectorWidth)+"° endocortical radius [mm]\t";
			}
			for (int i = 0;i<((int) 360/sectorWidth);++i){
				headings+=i*sectorWidth+"° - "+((i+1)*sectorWidth)+"° pericortical radius [mm]\t";
			}
			//Cortex BMD values			
			for (int i = 0;i<((int) 360/sectorWidth);++i){
				headings+=i*sectorWidth+"° - "+((i+1)*sectorWidth)+"° endocortical vBMD [mg/cm³]\t";
			}
			for (int i = 0;i<((int) 360/sectorWidth);++i){
				headings+=i*sectorWidth+"° - "+((i+1)*sectorWidth)+"° midcortical vBMD [mg/cm³]\t";
			}
			for (int i = 0;i<((int) 360/sectorWidth);++i){
				headings+=i*sectorWidth+"° - "+((i+1)*sectorWidth)+"° pericortical vBMD [mg/cm³]\t";
			}

		}
		textPanel.setColumnHeadings(headings);
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

	String printResults(String results,DetermineAlfa determineAlfa, ImagePlus imp){
		String[] propertyNames = {"File Name","Patient's Name","Patient ID","Patient's Birth Date","Acquisition Date","Pixel Spacing"};
		String[] parameters = {Double.toString(fatThreshold),Double.toString(areaThreshold),Double.toString(BMDThreshold),Double.toString(scalingFactor),Double.toString(constant)};

		if (imp != null){
			if (getInfoProperty(imageInfo,"File Name")!= null){
				results+=getInfoProperty(imageInfo,"File Name")+"\t";
			}else{
				if(imp.getImageStackSize() == 1){
					results+=getInfoProperty(imageInfo,"Title")+"\t";
				}else{
					results+=imageInfo.substring(0,imageInfo.indexOf("\n"))+"\t";
				}
			}
			for (int i = 1;i<propertyNames.length;++i){
				results+=getInfoProperty(imageInfo,propertyNames[i])+"\t";
			}
		}
		
		for (int i = 0;i<parameters.length;++i){
			results+=parameters[i]+"\t";
		}
		
		results += Double.toString(determineAlfa.alfa*180/Math.PI)+"\t";
		results += Double.toString(determineAlfa.rotationCorrection)+"\t";
		results += Boolean.toString(manualRotation)+"\t";
		results += Boolean.toString(flipDistribution)+"\t";
		results += Boolean.toString(guessFlip)+"\t";
		results += Boolean.toString(stacked)+"\t";
		results += Boolean.toString(allowCleaving)+"\t";
		results += roiChoice+"\t";
		results += rotationChoice+"\t";
		return results;
	}
	
	String printCorticalResults(String results,CorticalAnalysis cortAnalysis){
		results+=cortAnalysis.BMD+"\t";
		results+=cortAnalysis.AREA+"\t";
		results+=cortAnalysis.SSI+"\t";
		results+=cortAnalysis.ToD+"\t";
		results+=cortAnalysis.ToA+"\t";
		results+=cortAnalysis.BSId+"\t";
		return results;
	}
	
	String printMassDistributionResults(String results,MassDistribution massDistribution){
		for (int pp = 0;pp<((int) 360/sectorWidth);pp++){
			results += massDistribution.BMCs[pp]+"\t";
		}
		return results;
	}		
	
	String printDistributionResults(String results,DistributionAnalysis DistributionAnalysis){
		for (int pp = 0;pp<((int) 360/sectorWidth);pp++){
			results += DistributionAnalysis.endocorticalRadii[pp]+"\t";
		}
		for (int pp = 0;pp<((int) 360/sectorWidth);pp++){
			results += DistributionAnalysis.pericorticalRadii[pp]+"\t";
		}
		//Cortex BMD values			
		for (int pp = 0;pp<((int) 360/sectorWidth);pp++){
			results += DistributionAnalysis.endoCorticalBMDs[pp]+"\t";
		}
		for (int pp = 0;pp<((int) 360/sectorWidth);pp++){
			results += DistributionAnalysis.midCorticalBMDs[pp]+"\t";
		}
		for (int pp = 0;pp<((int) 360/sectorWidth);pp++){
			results += DistributionAnalysis.periCorticalBMDs[pp]+"\t";
		}
		return results;
	}
}
