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
	boolean alphaOn;	//Rotation angle
	String resultString;
	String imageInfo;
	boolean flipHorizontal;
	boolean flipVertical;
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
	boolean sleeveOn;
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
			if (cal.getCoefficients() != null){
				calibrationCoefficients = cal.getCoefficients();
			}
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
		String[] topLabels = new String[4];
		boolean[] defaultTopValues = new boolean[4];
		topLabels[0] = "Flip_horizontal";
		defaultTopValues[0] = false;
		topLabels[1] = "Flip_vertical";
		defaultTopValues[1] = false;
		topLabels[2] = "No_filtering";
		defaultTopValues[2] = false;
		topLabels[3] = "Measurement_tube";
		defaultTopValues[3] = false;
		dialog.addCheckboxGroup(1, 4, topLabels, defaultTopValues);
		
		dialog.addNumericField("Air_threshold", -40, 4, 8, null);	//Anything above this is fat or more dense
		dialog.addNumericField("Fat threshold", 40, 4, 8, null);		//Anything between this and air threshold is fat
		dialog.addNumericField("Muscle_threshold", 40, 4, 8, null);		//Anything above this is muscle or more dense
		dialog.addNumericField("Marrow_threshold", 80, 4, 8, null);		//Anything above this is muscle or more dense		
		dialog.addNumericField("Soft_tissue_threshold", 200.0, 4, 8, null);		//Anything  between this and muscle threshold is muscle
		dialog.addNumericField("Rotation_threshold", 200.0, 4, 8, null);
		dialog.addNumericField("Area threshold", 550.0, 4, 8, null); 	//550.0
		dialog.addNumericField("BMD threshold", 690.0, 4, 8, null);		//690.0
		
		
		dialog.addNumericField("Scaling_coefficient (slope)", calibrationCoefficients[1], 4, 8, null);
		dialog.addNumericField("Scaling_constant (intercept)",calibrationCoefficients[0], 4, 8, null);
		
		/*
		//Debugging
		dialog.addNumericField("Scaling_coefficient (slope)", 0.821, 4, 8, null);
		dialog.addNumericField("Scaling_constant (intercept)",-856.036, 4, 8, null);
		*/
		//Get ROI selection
		String[] choiceLabels = {"Bigger","Smaller","Left","Right","Top","Bottom","Central","Peripheral","SecondLargest","TwoLargestLeft","TwoLargestRight"};
		dialog.addChoice("Roi_selection", choiceLabels, choiceLabels[0]); 
		dialog.addChoice("Soft_Tissue_Roi_selection", choiceLabels, choiceLabels[0]); 
		String[] rotationLabels = {"According_to_Imax/Imin","Furthest_point","All_Bones_Imax/Imin","Not_selected_to_right","Selected_to_right"};
		dialog.addChoice("Rotation_selection", rotationLabels, rotationLabels[0]); //"According_to_Imax/Imin"
		
		String[] middleLabels = new String[10];
		boolean[] middleDefaults = new boolean[10];
		middleLabels[0] = "Analyse_cortical_results";
		middleDefaults[0] = false;
		middleLabels[1] = "Analyse_mass_distribution";
		middleDefaults[1] = false;
		middleLabels[2] = "Analyse_concentric_density_distribution";
		middleDefaults[2] = false;
		middleLabels[3] = "Analyse_density_distribution";
		middleDefaults[3] = false;
		middleLabels[4] = "Analyse_soft_tissues";
		middleDefaults[4] = false;
		middleLabels[5] = "Prevent_peeling_PVE_pixels";
		middleDefaults[5] = false;
		middleLabels[6] = "Allow_cleaving";
		middleDefaults[6] = false;
		middleLabels[7] = "Suppress_result_image";
		middleDefaults[7] = false;
		middleLabels[8] = "Limit_ROI_search_to_manually_selected";
		middleDefaults[8] = false;
		middleLabels[9] = "Set_distribution_results_rotation_manually";
		middleDefaults[9] = false;
		dialog.addCheckboxGroup(4, 3, middleLabels, middleDefaults);
		
		dialog.addNumericField("Manual_rotation_[+-_180_deg]", 0.0, 4, 8, null);
		
		String[] bottomLabels = new String[7];
		boolean[] bottomDefaults = new boolean[7];
		bottomLabels[0] = "Guess_right";
		bottomDefaults[0] = false;
		bottomLabels[1] = "Guess_larger";
		bottomDefaults[1] = false;
		bottomLabels[2] = "Stacked_bones";
		bottomDefaults[2] = false;
		bottomLabels[3] = "Guess_stacked";
		bottomDefaults[3] = false;
		bottomLabels[4] = "Invert_flip_guess";
		bottomDefaults[4] = false;
		bottomLabels[5] = "Flip_distribution_results";
		bottomDefaults[5] = false;
		bottomLabels[6] = "Save_visual_result_image_on_disk";
		bottomDefaults[6] = false;
		dialog.addCheckboxGroup(2, 5, bottomLabels, bottomDefaults);
		
		dialog.addStringField("Image_save_path",Prefs.getDefaultDirectory(),40);
		dialog.addHelp("http://bonej.org/densitydistribution");
		dialog.showDialog();
		
		if (dialog.wasOKed()){ //Stop in case of cancel..
			flipHorizontal				= dialog.getNextBoolean();
			flipVertical				= dialog.getNextBoolean();
			boolean noFiltering			= dialog.getNextBoolean();
			sleeveOn					= dialog.getNextBoolean();
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
			
			/*Look study special, images acquired prior to 2007 need to be flipped horizontally
			String checkDate = getInfoProperty(imageInfo,"Acquisition Date");
			if (checkDate.indexOf("2005") >-1 || checkDate.indexOf("2006") >-1){
				flipHorizontal = true;
			}
			*/
			short[] tempPointer = (short[]) imp.getProcessor().getPixels();			
			int[] signedShort = new int[tempPointer.length];			
		
			if (imp.getOriginalFileInfo().fileType == ij.io.FileInfo.GRAY16_SIGNED || cal.isSigned16Bit()){
				float[] floatPointer = (float[]) imp.getProcessor().toFloat(1,null).getPixels();
				for (int i=0;i<tempPointer.length;++i){signedShort[i] = (int) (floatPointer[i] - Math.pow(2.0,15.0));}
			} else {
				/*
				Apply the original calibration of the image prior to applying the calibration got from the user
				-> enables using ImageJ for figuring out the calibration without too much fuss.
				*/
				double[] origCalCoeffs = imp.getOriginalFileInfo().coefficients;
				if (origCalCoeffs == null){origCalCoeffs = cal.getCoefficients();}
				float[] floatPointer = (float[]) imp.getProcessor().toFloat(1,null).getPixels();
				for (int i=0;i<tempPointer.length;++i){signedShort[i] = (int) (floatPointer[i]*origCalCoeffs[1]+origCalCoeffs[0]);}
			}
			
			ImageAndAnalysisDetails imageAndAnalysisDetails = new ImageAndAnalysisDetails(flipHorizontal,flipVertical,noFiltering,sleeveOn
															,scalingFactor, constant,
															airThreshold, fatThreshold, muscleThreshold,marrowThreshold, softThreshold,	rotationThreshold, areaThreshold, BMDThreshold,
															roiChoice,roiChoiceSt,rotationChoice,choiceLabels,rotationLabels,
															preventPeeling,allowCleaving,manualRoi,manualRotation,manualAlfa,flipDistribution,
															guessFlip,guessLarger, stacked,guessStacked,invertGuess,sectorWidth,divisions,concentricSector,concentricDivisions,stOn);
			scaledImageData = new ScaledImageData(signedShort, imp.getWidth(), imp.getHeight(),resolution, scalingFactor, constant,3,flipHorizontal,flipVertical,noFiltering);	//Scale and 3x3 median filter the data
			RoiSelector roi = null;
			if(cOn || mOn || conOn || dOn){
				roi = new SelectROI(scaledImageData, imageAndAnalysisDetails,imp,imageAndAnalysisDetails.boneThreshold,true);
			}
			RoiSelector softRoi = null;
			if(stOn){
				softRoi = new SelectSoftROI(scaledImageData, imageAndAnalysisDetails,imp,imageAndAnalysisDetails.boneThreshold,true);
				if (roi == null){
					roi = softRoi;
				}
			}
			
			if (roi != null){	/*An analysis was conducted*/
				alphaOn = false;
				DetermineAlfa determineAlfa = null;
				if (cOn || mOn || conOn || dOn){
					determineAlfa = new DetermineAlfa((SelectROI) roi,imageAndAnalysisDetails);
					alphaOn = true;
				}
				
				imageAndAnalysisDetails.flipDistribution = roi.details.flipDistribution;
				flipDistribution = imageAndAnalysisDetails.flipDistribution;
				imageAndAnalysisDetails.stacked = roi.details.stacked;
				stacked = imageAndAnalysisDetails.stacked;
				TextPanel textPanel = IJ.getTextPanel();
				if (textPanel == null) {textPanel = new TextPanel();}
				if (textPanel.getLineCount() == 0){writeHeader(textPanel);}
				
				String results = "";
				results = printResults(results, imp);
				if (determineAlfa != null){
					results = printAlfa(results,determineAlfa);
				}
				
				ImagePlus resultImage = null;
				boolean makeImage = true;
				if(suppressImages && !saveImageOnDisk && roi != null){
					makeImage = false;
				}else{
					resultImage = getRGBResultImage(roi.scaledImage,roi.width,roi.height);
				}

				if(stOn){
					SoftTissueAnalysis softTissueAnalysis = new SoftTissueAnalysis((SelectSoftROI) softRoi);
					results = printSoftTissueResults(results,softTissueAnalysis);
					if(makeImage && resultImage != null){
						resultImage = addSoftTissueSieve(resultImage,softRoi.softSieve);
					}
				}
				
				if (cOn){
					CorticalAnalysis cortAnalysis =new CorticalAnalysis((SelectROI) roi);
					results = printCorticalResults(results,cortAnalysis);
					if(makeImage && resultImage != null){
						resultImage = addBoneSieve(resultImage,roi.sieve,roi.scaledImage,roi.details.marrowThreshold);
					}
					
				}
				if (mOn){
					MassDistribution massDistribution =new MassDistribution((SelectROI) roi,imageAndAnalysisDetails,determineAlfa);
					results = printMassDistributionResults(results,massDistribution);
				}
				if (conOn){
					ConcentricRingAnalysis concentricRingAnalysis =new ConcentricRingAnalysis((SelectROI) roi,imageAndAnalysisDetails,determineAlfa);
					results = printConcentricRingResults(results,concentricRingAnalysis);
					if(!dOn && makeImage && resultImage != null){
						resultImage = addPeriRadii(resultImage,concentricRingAnalysis.boneCenter, determineAlfa.pindColor,concentricRingAnalysis.Ru,concentricRingAnalysis.Theta);
						resultImage = addMarrowCenter(resultImage,determineAlfa.alfa/Math.PI*180.0,concentricRingAnalysis.boneCenter);
					}
				}
				
				
				if (dOn){
					DistributionAnalysis DistributionAnalysis = new DistributionAnalysis((SelectROI) roi,imageAndAnalysisDetails,determineAlfa);
					results = printDistributionResults(results,DistributionAnalysis);
					if (makeImage && resultImage != null){
						resultImage = addRadii(resultImage,determineAlfa.alfa/Math.PI*180.0,DistributionAnalysis.marrowCenter, determineAlfa.pindColor,DistributionAnalysis.R,DistributionAnalysis.R2,DistributionAnalysis.Theta);
						resultImage = addMarrowCenter(resultImage,determineAlfa.alfa/Math.PI*180.0,DistributionAnalysis.marrowCenter);
					}
				}
				
				if ((dOn || conOn) && makeImage && resultImage != null){
					resultImage = addRotate(resultImage,determineAlfa.alfa/Math.PI*180.0);
				}
				
				if (!suppressImages && resultImage!= null){
					resultImage = addScale(resultImage,roi.pixelSpacing);	//Add scale after rotating
					resultImage.show();
				}
				if (saveImageOnDisk && resultImage!= null){
					resultImage = addScale(resultImage,roi.pixelSpacing);	//Add scale after rotating
					FileSaver fSaver = new FileSaver(resultImage);
					fSaver.saveAsPng(imageSavePath+"/"+imageName+".png"); 
				}
				textPanel.appendLine(results);
				textPanel.updateDisplay();
			}else{
				IJ.log("No analysis was selected.");
			}
		}
	}

	/*Get image into which we'll start adding stuff*/
	ImagePlus getRGBResultImage(double[] values,int width,int height){
		ImagePlus tempImage = new ImagePlus("Visual results");
		tempImage.setProcessor(new FloatProcessor(width,height,values));
		new ImageConverter(tempImage).convertToRGB();
		return tempImage;
	}
	
	ImagePlus addScale(ImagePlus tempImage, double pixelSpacing){
		Calibration cal = new Calibration();
		cal.setUnit("mm");
		cal.pixelWidth = cal.pixelHeight = pixelSpacing;
		tempImage.setCalibration(cal);
		tempImage.getProcessor().setColor(new Color(255,0,0));
		tempImage.getProcessor().drawLine(5, 5, (int)(5.0+10.0/pixelSpacing), 5);
		tempImage.getProcessor().drawString("1 cm", 5, 20);
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
				if (sieve[x+y*tempImage.getWidth()] == 4){   //Tint intra fat area with green
					int value = tempImage.getProcessor().getPixel(x,y);
					int[] rgb = new int[3];
					for (int i = 0; i<3;++i){
						rgb[i] = (value >>(i*8))& 0XFF;
					}
					tempImage.getProcessor().setColor(new Color(0,rgb[1],0));
					tempImage.getProcessor().drawPixel(x,y);
				}
				if (sieve[x+y*tempImage.getWidth()] == 5){   //Tint subcut fat area with purple
					int value = tempImage.getProcessor().getPixel(x,y);
					int[] rgb = new int[3];
					for (int i = 0; i<3;++i){
						rgb[i] = (value >>(i*8))& 0XFF;
					}
					tempImage.getProcessor().setColor(new Color(rgb[2],0,rgb[0]));
					tempImage.getProcessor().drawPixel(x,y);
				}
			}
		}
		//tempImage.setProcessor(tempImage.getProcessor().resize(1000));
		return tempImage;
	}
	
	/*Add bone sieve*/
	ImagePlus addBoneSieve(ImagePlus tempImage, byte[] sieve,double[] scaledImage, double marrowThreshold){
		for (int y = 0; y < tempImage.getHeight();++y) {
			for (int x = 0; x < tempImage.getWidth();++x) {
				if (sieve[x+y*tempImage.getWidth()] == 1){   //Tint bone area with purple
					int value = tempImage.getProcessor().getPixel(x,y);
					int[] rgb = new int[3];
					for (int i = 0; i<3;++i){
						rgb[i] = (value >>(i*8))& 0XFF;
					}
					tempImage.getProcessor().setColor(new Color(rgb[2],0,rgb[0]));
					tempImage.getProcessor().drawPixel(x,y);
				}
				if (sieve[x+y*tempImage.getWidth()] == 1 && scaledImage[x+y*tempImage.getWidth()] <=marrowThreshold){   //Tint marrow area with green
					int value = tempImage.getProcessor().getPixel(x,y);
					int[] rgb = new int[3];
					for (int i = 0; i<3;++i){
						rgb[i] = (value >>(i*8))& 0XFF;
					}
					if (rgb[0] < 255-50){
						rgb[0]+=50;
					}
					tempImage.getProcessor().setColor(new Color(0,0,rgb[0]));
					tempImage.getProcessor().drawPixel(x,y);
				}
			}
		}
		return tempImage;
	}
		
	/*addDenstiyDistribution*/
	ImagePlus addRadii(ImagePlus tempImage,double alfa,double[] marrowCenter,Vector<Integer> pindColor,
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
		return tempImage;
	}
	ImagePlus addMarrowCenter(ImagePlus tempImage,double alfa,double[] marrowCenter){
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
		return tempImage;
	}
	
	ImagePlus addRotate(ImagePlus tempImage,double alfa){
		tempImage.getProcessor().setBackgroundValue(0.0);
		IJ.run(tempImage, "Arbitrarily...", "angle=" + alfa + " grid=1 interpolation=Bilinear enlarge");  
		return tempImage;
	}
	
	/*Concentric rings distribution result image*/
	ImagePlus addPeriRadii(ImagePlus tempImage,double[] marrowCenter,Vector<Integer> pindColor, double[] R, double[] Theta){
		//Draw unrotated radii
		for(int i = 0; i< Theta.length;i++) {//45;i++) {//
			int x = ((int) (marrowCenter[0]+R[i]*Math.cos(Theta[i])));
			int y = ((int) (marrowCenter[1]+R[i]*Math.sin(Theta[i])));
			double colorScale = ((double) pindColor.get(i))/359.0;
			tempImage.getProcessor().setColor(new Color(0,(int) (255.0*colorScale),(int) (255.0*(1.0-colorScale))));
			tempImage.getProcessor().drawPixel(x,y);
		}
		return tempImage;
	}
	
	void writeHeader(TextPanel textPanel){
		String[] propertyNames = {"File Name","Patient's Name","Patient ID","Patient's Birth Date","Acquisition Date","Pixel Spacing","Object Length"};
		String[] parameterNames = {"Air Threshold","Fat Threshold","Muscle Threshold","Marrow Threshold","Soft Threshold","Rotation Threshold","Area Threshold","BMD Threshold","Scaling Coefficient","Scaling Constant"};
		String[] dHeadings = {"Manual Rotation","Flip Distribution","Guess right","Guess larger"
		,"Stacked bones","Invert guess","Allow Cleaving","Prevent PVE peeling","Roi choice","Rotation choice","Flip Horizontal","Flip Vertical"};
		
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
		if(alphaOn){
			String[] rHeadings = {"Alpha [deg]","Rotation correction [deg]","Distance between bones[mm]"};	
			for (int i = 0;i<rHeadings.length;++i){
				headings+=rHeadings[i]+"\t";
			}
		}
		
		if(stOn){
			String[] coHeadings = {"MuD [mg/cm³]","MuA [cm²]","LeanMuD [mg/cm³]","LeanMuA [cm²]","IntraFatD [mg/cm³]","IntraFatA [cm²]","FatD [mg/cm³]","FatA [cm²]","SubCutFatD [mg/cm³]","SubCutFatA [cm²]","LimbD [mg/cm³]","LimbA [cm²]","Density weighted fat percentage [%]"};
			for (int i = 0;i<coHeadings.length;++i){
				headings+=coHeadings[i]+"\t";
			}
		}
		
		if(cOn){
			String[] coHeadings = {"MaMassD [g/cm³]","StratecMaMassD [g/cm³]","MaD [mg/cm³]","MaA [mm²]","CoD [mg/cm³]","CoA [mm²]","SSI [mm³]","ToD [mg/cm³]","ToA[mm²]","MeA [mm²]","BSId[g²/cm4]"};
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
			headings+="Peeled mean vBMD [mg/cm³]\t";
			//Radial distribution
			for (int i =0; i < (int) divisions; ++i){
				headings+= "Radial division "+i+" vBMD [mg/cm³]\t";
			}
			//Polar distribution
			for (int i = 0;i<((int) 360/sectorWidth);++i){
				headings+= "Polar sector "+i+" vBMD [mg/cm³]\t";
			}
			
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

	String printResults(String results, ImagePlus imp){
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
		results += Boolean.toString(flipHorizontal)+"\t";
		results += Boolean.toString(flipVertical)+"\t";
		return results;
	}
	
	String printAlfa(String results,DetermineAlfa determineAlfa){
		results += Double.toString(determineAlfa.alfa*180/Math.PI)+"\t";
		results += Double.toString(determineAlfa.rotationCorrection)+"\t";
		results += Double.toString(determineAlfa.distanceBetweenBones)+"\t";
		return results;
	}
	String printSoftTissueResults(String results,SoftTissueAnalysis softTissueAnalysis){
		results+=softTissueAnalysis.TotalMuD+"\t";
		results+=softTissueAnalysis.TotalMuA+"\t";
		results+=softTissueAnalysis.MuD+"\t";
		results+=softTissueAnalysis.MuA+"\t";
		results+=softTissueAnalysis.IntraMuFatD+"\t";
		results+=softTissueAnalysis.IntraMuFatA+"\t";
		results+=softTissueAnalysis.FatD+"\t";
		results+=softTissueAnalysis.FatA+"\t";
		results+=softTissueAnalysis.SubCutFatD+"\t";
		results+=softTissueAnalysis.SubCutFatA+"\t";
		results+=softTissueAnalysis.LimbD+"\t";
		results+=softTissueAnalysis.LimbA+"\t";
		results+=softTissueAnalysis.FatPercentage+"\t";
		return results;
	}
	
	String printCorticalResults(String results,CorticalAnalysis cortAnalysis){
		results+=cortAnalysis.MaMassD+"\t";
		results+=cortAnalysis.StratecMaMassD+"\t";
		results+=cortAnalysis.MaD+"\t";
		results+=cortAnalysis.MaA+"\t";
		results+=cortAnalysis.BMD+"\t";
		results+=cortAnalysis.AREA+"\t";
		results+=cortAnalysis.SSI+"\t";
		results+=cortAnalysis.ToD+"\t";
		results+=cortAnalysis.ToA+"\t";
		results+=cortAnalysis.MeA+"\t";
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
		results+= DistributionAnalysis.peeledBMD+"\t";
		//Radial distribution
		for (int i =0; i < (int) divisions; ++i){
			results+= DistributionAnalysis.radialDistribution[i]+"\t";
		}
		//Polar distribution
		for (int i = 0;i<((int) 360/sectorWidth);++i){
			results+= DistributionAnalysis.polarDistribution[i]+"\t";
		}
		
		
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
