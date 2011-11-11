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

package org.doube.bonej.pqct.io;
public class ImageAndAnalysisDetails{
	public double scalingFactor;
	public double constant;
	public double marrowThreshold;
	public double airThreshold;
	public double fatThreshold;
	public double muscleThreshold;
	public double areaThreshold;	//For cortical AREA analyses (CoA, SSI, I) + peeling distal pixels
	public double BMDthreshold;		//For cortical BMD analyses
	public double softThreshold;	//Thresholding soft tissues + marrow from bone
	public double boneThreshold;	//Thresholding bone from the rest and cortical AREA analyses (CoA, SSI, I)
	public int filterSize;
	public int softFilterSize;
	public int sectorWidth;
	public int divisions;
	public int concentricSector;
	public int concentricDivisions;
	public String imageSavePath;
	public String roiChoice;
	public String rotationChoice;
	public String[] choiceLabels;
	public boolean allowCleaving;
	public boolean manualRoi;
	public boolean manualRotation;
	public double manualAlfa;
	public boolean flipDistribution;
	public boolean guessFlip;
	public boolean guessLarger;
	public boolean stacked;
	public boolean invertGuess;
	
	//ImageJ plugin constructor
	public ImageAndAnalysisDetails(double scalingFactor, double constant,double fatThreshold,double areaThreshold,
									double BMDthreshold, String roiChoice,String rotationChoice,String[] choiceLabels,
									boolean allowCleaving, boolean manualRoi,
									boolean manualRotation, double manualAlfa, boolean flipDistribution, 
									boolean guessFlip,boolean guessLarger,boolean stacked, boolean invertGuess,
									int sectorWidth,int divisions,int concentricSector,int concentricDivisions){
		this.scalingFactor			= scalingFactor;
		this.constant 				= constant;
		this.airThreshold			= -100;
		this.fatThreshold	 		= fatThreshold;
		this.muscleThreshold 		= 200;
		this.marrowThreshold 		= 300;
		this.areaThreshold 			= areaThreshold;	//For cortical AREA analyses (CoA, SSI, I) + peeling distal pixels
		this.BMDthreshold 			= BMDthreshold;		//For cortical BMD analyses
		this.softThreshold 			= 300;	//Thresholding soft tissues + marrow from bone
		this.boneThreshold 			= areaThreshold;
		this.filterSize				= 3;
		this.softFilterSize			= 7;
		this.sectorWidth			=sectorWidth;
		this.divisions				= divisions;
		this.concentricSector		= concentricSector;
		this.concentricDivisions	= concentricDivisions;
		this.roiChoice				= roiChoice;
		this.rotationChoice			= rotationChoice;
		this.choiceLabels			= choiceLabels;
		this.imageSavePath 			= new String("");
		this.allowCleaving			= allowCleaving;
		this.manualRoi				= manualRoi;
		this.manualRotation			= manualRotation;
		this.manualAlfa				= manualAlfa;
		this.flipDistribution		= flipDistribution;
		this.guessFlip				= guessFlip;
		this.guessLarger			= guessLarger;
		this.stacked				= stacked;
		this.invertGuess			= invertGuess;
		this.sectorWidth			=sectorWidth;
		this.divisions				= divisions;
		this.concentricSector		= concentricSector;
		this.concentricDivisions	= concentricDivisions;
	}
}
