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
import java.io.*;							//File IO for reading .TYP files
import ij.plugin.PlugIn;
import org.doube.bonej.pqct.analysis.*;		//Analysis stuff..
import org.doube.bonej.pqct.selectroi.*;	//ROI selection..
import org.doube.bonej.pqct.io.*;			//image data 
import java.awt.*;							//Image, component for debugging...
import ij.plugin.filter.Info;
import ij.io.*;

public class Distribution_Analysis implements PlugIn {

	
	int sectorWidth;
	int divisions;
	int concentricSector;
	int concentricDivisions;
	boolean cOn;		//Basic analyses
	boolean	mOn;		//Mass distribution
	boolean	conOn;		//Concentric rings analysis
	boolean dOn;		//Distribution analysis
	boolean stOn;		//Soft tissue analysis
	String resultString;
	String imageInfo;
	boolean flipHorizontal;
	double resolution;
	//Thresholds
	double airThreshold;
	double fatThreshold;
	double muscleThreshold;
	double marrowThreshold;
	double softThreshold;	
	double rotationThreshold;
	double areaThreshold;
	double BMDThreshold;

	double scalingFactor;
	double constant;	
	boolean flipDistribution;
	boolean guessFlip;
	boolean guessLarger;
	boolean stacked;
	boolean guessStacked;
	boolean invertGuess;
	boolean manualRotation;
	boolean allowCleaving;
	boolean preventPeeling;
	String roiChoice;
	String roiChoiceSt;
	String rotationChoice;
	
	
	public void run(String arg) {
		ImagePlus imp = WindowManager.getCurrentImage();
		if (imp == null)
			return;
		if (imp.getType() != ImagePlus.GRAY16){
			IJ.error("Distribution analysis expects 16-bit greyscale data");
			return;
		}
		/*Set sector widths and division numbers*/
		sectorWidth 		= 10;
		divisions			= 3;
		concentricSector	= 10;
		concentricDivisions = 10;		
		
		imageInfo = new Info().getImageInfo(imp,imp.getChannelProcessor());
		/*Check image calibration*/
		Calibration cal = imp.getCalibration();
		double[] calibrationCoefficients = {0,1};
		if (getInfoProperty(imageInfo,"Stratec File") == null){ 
				if (cal.getCoefficients() != null)
			calibrationCoefficients = cal.getCoefficients();
		} else {
			calibrationCoefficients = new double[2];
			/*Read calibration from TYP file database*/
			String typFileName = getInfoProperty(imageInfo,"Device");
			try {
				InputStream ir = this.getClass().getClassLoader().getResourceAsStream("org/doube/bonej/pqct/typ/"+typFileName);
				byte[] typFileData = new byte[ir.available()];
				ir.read(typFileData);
				ir.close();
				String typFiledDataString = new String(typFileData,"ISO-8859-1");
				//break the typFileDataString into lines
				StringTokenizer st = new StringTokenizer(typFiledDataString, "\n");
				Vector<String> typFileLines = new Vector<String>();
				while(st.hasMoreTokens()){
					typFileLines.add(st.nextToken());
				}
				//Search for XSlope and XInter
				String[] searchFor = {"XInter","XSlope"};
				for (int i = 0;i<searchFor.length;++i){
					int index = 0;
					String temp = typFileLines.get(index);
					while (temp.indexOf(searchFor[i]) < 0 && index < typFileLines.size()){
						++index;
						temp = typFileLines.get(index);
					}
					if (temp.indexOf(searchFor[i]) >= 0){	//Found line
						StringTokenizer st2 = new StringTokenizer(temp, "=");
						Vector<String> typFileLineTokens = new Vector<String>();
						while(st2.hasMoreTokens()){
							typFileLineTokens.add(st2.nextToken().trim());
						}
						calibrationCoefficients[i] = Double.valueOf(typFileLineTokens.get(1));
					} else {
						calibrationCoefficients[i] = (double) i*1000.0;
					}
				}
				calibrationCoefficients[1] /= 1000.0;		//1.495
			} catch (Exception err){System.err.println("Error: "+err.getMessage());}
		}
		
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
		dialog.addCheckbox("Flip_horizontal",false);
		dialog.addCheckbox("No_filtering",false);
		dialog.addNumericField("Air_threshold", -155.3, 4, 8, null);	//Anything above this is fat or more dense
		dialog.addNumericField("Fat threshold", -23.295, 4, 8, null);		//Anything between this and air threshold is fat
		dialog.addNumericField("Muscle_threshold", 0.0, 4, 8, null);		//Anything above this is muscle or more dense
		dialog.addNumericField("Marrow_threshold", 54.35493, 4, 8, null);		//Anything above this is muscle or more dense		
		dialog.addNumericField("Soft_tissue_threshold", 155.2998, 4, 8, null);		//Anything  between this and muscle threshold is muscle
		dialog.addNumericField("Rotation_threshold", 200.0, 4, 8, null);
		dialog.addNumericField("Area threshold", 600.0, 4, 8, null); 	//550.0
		dialog.addNumericField("BMD threshold", 600.0, 4, 8, null);		//690.0
		/*
		dialog.addNumericField("Scaling_coefficient (slope)", calibrationCoefficients[1], 4, 8, null);
		dialog.addNumericField("Scaling_constant (intercept)",calibrationCoefficients[0], 4, 8, null);
		*/
		dialog.addNumericField("Scaling_coefficient (slope)", 0.743, 4, 8, null);
		dialog.addNumericField("Scaling_constant (intercept)",-751.873, 4, 8, null);
		//Get ROI selection
		String[] choiceLabels = {"Bigger","Smaller","Left","Right","Top","Bottom","Central","Peripheral","SecondLargest","TwoLargestLeft","TwoLargestRight"};
		dialog.addChoice("Roi_selection", choiceLabels, choiceLabels[2]); 
		dialog.addChoice("Soft_Tissue_Roi_selection", choiceLabels, choiceLabels[9]); 
		String[] rotationLabels = {"According_to_Imax/Imin","Furthest_point","All_Bones_Imax/Imin","Not_selected_to_right","Selected_to_right"};
		dialog.addChoice("Rotation_selection", rotationLabels, rotationLabels[1]); //"According_to_Imax/Imin"
		dialog.addCheckbox("Analyse_cortical_results",true);
		dialog.addCheckbox("Analyse_mass_distribution",true);
		dialog.addCheckbox("Analyse_concentric_density_distribution",false);
		dialog.addCheckbox("Analyse_density_distribution",true);	//true
		dialog.addCheckbox("Analyse_soft_tissues",true);	//true
		dialog.addCheckbox("Prevent_peeling_PVE_pixels",false);	//true
		dialog.addCheckbox("Allow_cleaving",true);					//false
		dialog.addCheckbox("Suppress_result_image",false);
		dialog.addCheckbox("Limit_ROI_search_to_manually_selected",false);
		dialog.addCheckbox("Set_distribution_results_rotation_manually",true);
		dialog.addNumericField("Manual_rotation_[+-_180_deg]", 0.0, 4, 8, null);
		dialog.addCheckbox("Flip_distribution_results",false);
		dialog.addCheckbox("Guess_right",false);
		dialog.addCheckbox("Guess_larger",false);
		dialog.addCheckbox("Stacked_bones",false);
		dialog.addCheckbox("Guess_stacked",false);
		dialog.addCheckbox("Invert_flip_guess",false);
		dialog.addCheckbox("Save_visual_result_image_on_disk",false);
		dialog.addStringField("Image_save_path",Prefs.getDefaultDirectory(),40);
		dialog.addHelp("http://bonej.org/densitydistribution");
		dialog.showDialog();
		
		if (dialog.wasOKed()){ //Stop in case of cancel..
			flipHorizontal				= dialog.getNextBoolean();
			boolean noFiltering			= dialog.getNextBoolean();
			airThreshold				= dialog.getNextNumber();
			fatThreshold				= dialog.getNextNumber();
			muscleThreshold				= dialog.getNextNumber();
			marrowThreshold				= dialog.getNextNumber();
			softThreshold				= dialog.getNextNumber();
			rotationThreshold			= dialog.getNextNumber();
			areaThreshold				= dialog.getNextNumber();
			BMDThreshold				= dialog.getNextNumber();
			scalingFactor				= dialog.getNextNumber();
			constant					= dialog.getNextNumber();
			roiChoice					= dialog.getNextChoice();
			roiChoiceSt					= dialog.getNextChoice();
			rotationChoice				= dialog.getNextChoice();
			cOn							= dialog.getNextBoolean();
			mOn							= dialog.getNextBoolean();
			conOn						= dialog.getNextBoolean();
			dOn							= dialog.getNextBoolean();
			stOn						= dialog.getNextBoolean();
			preventPeeling				= dialog.getNextBoolean();
			allowCleaving				= dialog.getNextBoolean();
			boolean suppressImages		= dialog.getNextBoolean();
			boolean manualRoi			= dialog.getNextBoolean();
			manualRotation				= dialog.getNextBoolean();
			double manualAlfa			= dialog.getNextNumber()*Math.PI/180.0;
			flipDistribution			= dialog.getNextBoolean();
			guessFlip					= dialog.getNextBoolean();
			guessLarger					= dialog.getNextBoolean();
			stacked						= dialog.getNextBoolean();
			guessStacked				= dialog.getNextBoolean();
			invertGuess					= dialog.getNextBoolean();
			boolean saveImageOnDisk		= dialog.getNextBoolean();
			String imageSavePath 		= dialog.getNextString();
			ScaledImageData scaledImageData;
			
			String imageName;
			if (getInfoProperty(imageInfo,"File Name")!= null){
				imageName = getInfoProperty(imageInfo,"File Name");
			}else{
				if(imp.getImageStackSize() == 1){
					imageName = imp.getTitle();
					imageInfo+="File Name:"+imageName+"\n";
				}else{
					imageName = imageInfo.substring(0,imageInfo.indexOf("\n"));
					imageInfo+="File Name:"+imageName+"\n";
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
			
			ImageAndAnalysisDetails imageAndAnalysisDetails = new ImageAndAnalysisDetails(flipHorizontal,noFiltering,scalingFactor, constant,
															airThreshold, fatThreshold, muscleThreshold,marrowThreshold, softThreshold,	rotationThreshold, areaThreshold, BMDThreshold,
															roiChoice,roiChoiceSt,rotationChoice,choiceLabels,rotationLabels,
															preventPeeling,allowCleaving,manualRoi,manualRotation,manualAlfa,flipDistribution,
															guessFlip,guessLarger, stacked,guessStacked,invertGuess,sectorWidth,divisions,concentricSector,concentricDivisions,stOn);
			scaledImageData = new ScaledImageData(unsignedShort, imp.getWidth(), imp.getHeight(),resolution, scalingFactor, constant,3,flipHorizontal,noFiltering);	//Scale and 3x3 median filter the data
			SelectROI roi = new SelectROI(scaledImageData, imageAndAnalysisDetails,imp,imageAndAnalysisDetails.boneThreshold,true);
			/*testing*/
			/*
			ImagePlus tempImage = new ImagePlus("Sieve");
			tempImage.setProcessor(new ByteProcessor(roi.width,roi.height));
			tempImage.getProcessor().setBackgroundValue(0.0);
			tempImage.getProcessor().setValue(255.0);

			for (int y = 0; y < roi.height;++y) {
				for (int x = 0; x < roi.width;++x) {
					if (roi.sieve[x+y*roi.width] == 1){   //Tint roi area color with violet
						tempImage.getProcessor().drawPixel(x,y);
					}
				}
			}
			tempImage.show();
			*/
			
			
			DetermineAlfa determineAlfa = new DetermineAlfa(roi,imageAndAnalysisDetails);
			
			/*
			ImagePlus tempImage2 = new ImagePlus("Determine");
			tempImage2.setProcessor(new ByteProcessor(roi.width,roi.height));
			tempImage2.getProcessor().setBackgroundValue(0.0);
			tempImage2.getProcessor().setValue(255.0);

			for (int y = 0; y < roi.height;++y) {
				for (int x = 0; x < roi.width;++x) {
					if (roi.sieve[x+y*roi.width] == 1){   //Tint roi area color with violet
						tempImage2.getProcessor().drawPixel(x,y);
					}
				}
			}
			tempImage2.show();
			*/
			
			imageAndAnalysisDetails.flipDistribution = roi.details.flipDistribution;
			flipDistribution = imageAndAnalysisDetails.flipDistribution;
			imageAndAnalysisDetails.stacked = roi.details.stacked;
			stacked = imageAndAnalysisDetails.stacked;
			TextPanel textPanel = IJ.getTextPanel();
			if (textPanel == null) {textPanel = new TextPanel();}
			if (textPanel.getLineCount() == 0){writeHeader(textPanel);}
			
			String results = "";
			results = printResults(results,determineAlfa, imp);
			ImagePlus resultImage = null;
			boolean makeImage = true;
			if(suppressImages && !saveImageOnDisk){
				makeImage = false;
			}else{
				resultImage = getRGBResultImage(roi.scaledImage,roi.width,roi.height);
			}

			if(stOn){
				SoftTissueAnalysis softTissueAnalysis = new SoftTissueAnalysis(roi);
				results = printSoftTissueResults(results,softTissueAnalysis);
				if(makeImage){
					resultImage = addSoftTissueSieve(resultImage,roi.softSieve);
				}
			}
			
			if (cOn ){
				CorticalAnalysis cortAnalysis =new CorticalAnalysis(roi);
				results = printCorticalResults(results,cortAnalysis);
				if(makeImage){
					resultImage = addBoneSieve(resultImage,roi.sieve);
				}
				
			}
			if (mOn){
				MassDistribution massDistribution =new MassDistribution(roi,imageAndAnalysisDetails,determineAlfa);
				results = printMassDistributionResults(results,massDistribution);
				/*
				if(!dOn && !conOn && makeImage){
					resultImage = getResultImage(roi.scaledImage,roi.width,roi.height,roi.minimum,roi.maximum,roi.sieve,determineAlfa.alfa/Math.PI*180.0);

				}
				*/
			}
			if (conOn){
				ConcentricRingAnalysis concentricRingAnalysis =new ConcentricRingAnalysis(roi,imageAndAnalysisDetails,determineAlfa);
				results = printConcentricRingResults(results,concentricRingAnalysis);
				/*
				if(!dOn && makeImage){
					resultImage = getResultImage(roi.scaledImage,roi.width,roi.height,roi.minimum,roi.maximum,roi.sieve,determineAlfa.alfa/Math.PI*180.0,concentricRingAnalysis.boneCenter,determineAlfa.pind, determineAlfa.pindColor,concentricRingAnalysis.Ru,concentricRingAnalysis.Theta);

				}
				*/
			}
			
			
			if (dOn){
				DistributionAnalysis DistributionAnalysis = new DistributionAnalysis(roi,imageAndAnalysisDetails,determineAlfa);
				results = printDistributionResults(results,DistributionAnalysis);
				if (makeImage){
					resultImage = addDensityDistribution(resultImage,determineAlfa.alfa/Math.PI*180.0,DistributionAnalysis.marrowCenter,determineAlfa.pind, determineAlfa.pindColor,DistributionAnalysis.R,DistributionAnalysis.R2,DistributionAnalysis.Theta);
				}
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

	/*Get image into which we'll start adding stuff*/
	ImagePlus getRGBResultImage(double[] values,int width,int height){
		ImagePlus tempImage = new ImagePlus("Visual results");
		tempImage.setProcessor(new FloatProcessor(width,height,values));
		new ImageConverter(tempImage).convertToRGB();
		return tempImage;
	}
	
	/*Add soft sieve*/
	ImagePlus addSoftTissueSieve(ImagePlus tempImage, byte[] sieve){
		for (int y = 0; y < tempImage.getHeight();++y) {
			for (int x = 0; x < tempImage.getWidth();++x) {
				if (sieve[x+y*tempImage.getWidth()] == 2){   //Tint fat area with yellow
					int value = tempImage.getProcessor().getPixel(x,y);
					int[] rgb = new int[3];
					for (int i = 0; i<3;++i){
						rgb[i] = (value >>(i*8))& 0XFF;
					}
					tempImage.getProcessor().setColor(new Color(rgb[2],rgb[1],0));
					tempImage.getProcessor().drawPixel(x,y);
				}
				if (sieve[x+y*tempImage.getWidth()] == 3){   //Tint muscle area with red
					int value = tempImage.getProcessor().getPixel(x,y);
					int[] rgb = new int[3];
					for (int i = 0; i<3;++i){
						rgb[i] = (value >>(i*8))& 0XFF;
					}
					tempImage.getProcessor().setColor(new Color(rgb[2],0,0));
					tempImage.getProcessor().drawPixel(x,y);
				}
			}
		}
		//tempImage.setProcessor(tempImage.getProcessor().resize(1000));
		return tempImage;
	}
	
	/*Add bone sieve*/
	ImagePlus addBoneSieve(ImagePlus tempImage, byte[] sieve){
		for (int y = 0; y < tempImage.getHeight();++y) {
			for (int x = 0; x < tempImage.getWidth();++x) {
				if (sieve[x+y*tempImage.getWidth()] == 1){   //Tint bone area with purple
					int value = tempImage.getProcessor().getPixel(x,y);
					int[] rgb = new int[3];
					for (int i = 0; i<3;++i){
						rgb[i] = (value >>(i*8))& 0XFF;
					}
					tempImage.getProcessor().setColor(new Color(rgb[2],0,rgb[1]));
					tempImage.getProcessor().drawPixel(x,y);
				}
			}
		}
		return tempImage;
	}
		
	/*addDenstiyDistribution*/
	ImagePlus addDensityDistribution(ImagePlus tempImage,double alfa,double[] marrowCenter,Vector<Integer> pind,Vector<Integer> pindColor,
										double[] R, double[] R2, double[] Theta){
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
		//tempImage.setProcessor(tempImage.getProcessor().resize(1000));
		return tempImage;
	}

	/*Concentric rings distribution result image*/
	ImagePlus getResultImage(double[] values,int width,int height, double min, double max, byte[] sieve, double alfa,
							double[] marrowCenter,Vector<Integer> pind,Vector<Integer> pindColor, double[] R, double[] Theta){
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
		for(int i = 0; i< Theta.length;i++) {//45;i++) {//
			int x = ((int) (marrowCenter[0]+R[i]*Math.cos(Theta[i])));
			int y = ((int) (marrowCenter[1]+R[i]*Math.sin(Theta[i])));
			double colorScale = ((double) pindColor.get(i))/359.0;
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
		//tempImage.setProcessor(tempImage.getProcessor().resize(1000));
		return tempImage;
	}
	
	void writeHeader(TextPanel textPanel){
		String[] propertyNames = {"File Name","Patient's Name","Patient ID","Patient's Birth Date","Acquisition Date","Pixel Spacing","Object Length"};
		String[] parameterNames = {"Air Threshold","Fat Threshold","Muscle Threshold","Marrow Threshold","Soft Threshold","Rotation Threshold","Area Threshold","BMD Threshold","Scaling Coefficient","Scaling Constant"};
		String[] dHeadings = {"Alpha [deg]","Rotation correction [deg]","Distance between bones[mm]","Manual Rotation","Flip Distribution","Guess right","Guess larger"
		,"Stacked bones","Invert guess","Allow Cleaving","Prevent PVE peeling","Roi choice","Rotation choice"};
			
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
		
		if(stOn){
			String[] coHeadings = {"MuD [mg/cm³]","MuA [mm²]","FatD [mg/cm³]","FatA [mm²]","LimbD [mg/cm³]","LimbA [mm²]"};
			for (int i = 0;i<coHeadings.length;++i){
				headings+=coHeadings[i]+"\t";
			}
		}
		
		if(cOn){
			String[] coHeadings = {"MedMassD [g/cm³]","MeD [mg/cm³]","MeA [mm²]","CoD [mg/cm³]","CoA [mm²]","SSI [mm³]","ToD [mg/cm³]","ToA[mm²]","BSId[g²/cm4]"};
			for (int i = 0;i<coHeadings.length;++i){
				headings+=coHeadings[i]+"\t";
			}
		}
		if(mOn){
			for (int i = 0;i<((int) 360/sectorWidth);++i){
				headings+=i*sectorWidth+"° - "+((i+1)*sectorWidth)+"° mineral mass [mg]\t";
			}
		}
		
		if(conOn){
			for (int i = 0;i<((int) 360/concentricSector);++i){
				headings+=i*concentricSector+"° - "+((i+1)*concentricSector)+"° concentric analysis pericortical radius [mm]\t";
			}
			for (int j = 0;j<concentricDivisions;++j){
				for (int i = 0;i<((int) 360/concentricSector);++i){
					headings+="Division "+(j+1)+" sector "+i*concentricSector+"° - "+((i+1)*concentricSector)+"° vBMD [mg/cm³]\t";
				}
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
		String[] propertyNames = {"File Name","Patient's Name","Patient ID","Patient's Birth Date","Acquisition Date","Pixel Spacing","ObjLen"};
		String[] parameters = {Double.toString(airThreshold)
								,Double.toString(fatThreshold),Double.toString(muscleThreshold)
								,Double.toString(marrowThreshold)
								,Double.toString(softThreshold),Double.toString(rotationThreshold)
								,Double.toString(areaThreshold),Double.toString(BMDThreshold)
								,Double.toString(scalingFactor),Double.toString(constant)};

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
		results += Double.toString(determineAlfa.distanceBetweenBones)+"\t";
		results += Boolean.toString(manualRotation)+"\t";
		results += Boolean.toString(flipDistribution)+"\t";
		results += Boolean.toString(guessFlip)+"\t";
		results += Boolean.toString(guessLarger)+"\t";
		results += Boolean.toString(stacked)+"\t";
		results += Boolean.toString(invertGuess)+"\t";
		results += Boolean.toString(allowCleaving)+"\t";
		results += Boolean.toString(preventPeeling)+"\t";
		results += roiChoice+"\t";
		results += rotationChoice+"\t";
		return results;
	}
	
	String printSoftTissueResults(String results,SoftTissueAnalysis softTissueAnalysis){
		results+=softTissueAnalysis.MuD+"\t";
		results+=softTissueAnalysis.MuA+"\t";
		results+=softTissueAnalysis.FatD+"\t";
		results+=softTissueAnalysis.FatA+"\t";
		results+=softTissueAnalysis.LimbD+"\t";
		results+=softTissueAnalysis.LimbA+"\t";
		return results;
	}
	
	String printCorticalResults(String results,CorticalAnalysis cortAnalysis){
		results+=cortAnalysis.medMassD+"\t";
		results+=cortAnalysis.MeD+"\t";
		results+=cortAnalysis.MeA+"\t";
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
	
	
	String printConcentricRingResults(String results,ConcentricRingAnalysis concentricRingAnalysis){
		for (int i = 0;i<((int) 360/concentricSector);++i){
			results += concentricRingAnalysis.pericorticalRadii[i]+"\t";
		}
		for (int j = 0;j<concentricDivisions;++j){
			for (int i = 0;i<((int) 360/concentricSector);++i){
				results += concentricRingAnalysis.BMDs.get(j)[i]+"\t";
			}
		}
		return results;
	}
	

	
	String printDistributionResults(String results,DistributionAnalysis DistributionAnalysis){
		for (int pp = 0;pp<((int) 360/sectorWidth);++pp){
			results += DistributionAnalysis.endocorticalRadii[pp]+"\t";
		}
		for (int pp = 0;pp<((int) 360/sectorWidth);++pp){
			results += DistributionAnalysis.pericorticalRadii[pp]+"\t";
		}
		//Cortex BMD values			
		for (int pp = 0;pp<((int) 360/sectorWidth);++pp){
			results += DistributionAnalysis.endoCorticalBMDs[pp]+"\t";
		}
		for (int pp = 0;pp<((int) 360/sectorWidth);++pp){
			results += DistributionAnalysis.midCorticalBMDs[pp]+"\t";
		}
		for (int pp = 0;pp<((int) 360/sectorWidth);++pp){
			results += DistributionAnalysis.periCorticalBMDs[pp]+"\t";
		}
		return results;
	}
}
