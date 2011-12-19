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
	public boolean flipHorizontal;
	public boolean noFiltering;
	public boolean sleeveOn;
	public double scalingFactor;
	public double constant;
	
	public double airThreshold;		//Fat lower threshold
	public double fatThreshold;		//Fat higher threshold
	public double muscleThreshold;	//Muscle lower threshold
	public double marrowThreshold;	//Marrow higher threshold
	public double softThreshold;	//Soft tissues higher threshold
	public double areaThreshold;	//For cortical AREA analyses (CoA, SSI, I) + peeling distal pixels
	public double rotationThreshold;
	public double BMDthreshold;		//For cortical BMD analyses	
	public double boneThreshold;	//Thresholding bone from the rest and cortical AREA analyses (CoA, SSI, I)
	public int filterSize;
	public int softFilterSize;
	public int sectorWidth;
	public int divisions;
	public int concentricSector;
	public int concentricDivisions;
	public String imageSavePath;
	public String roiChoice;
	public String roiChoiceSt;
	public String rotationChoice;
	public String[] choiceLabels;
	public String[] rotationLabels;
	public boolean preventPeeling;
	public boolean allowCleaving;
	public boolean manualRoi;
	public boolean manualRotation;
	public double manualAlfa;
	public boolean flipDistribution;
	public boolean guessFlip;
	public boolean guessLarger;
	public boolean stacked;
	public boolean guessStacked;
	public boolean invertGuess;
	public boolean stOn;
	
	//ImageJ plugin constructor
	public ImageAndAnalysisDetails(boolean flipHorizontal,boolean noFiltering,boolean sleeveOn
									,double scalingFactor, double constant,
									double airThreshold,double fatThreshold, double muscleThreshold, double marrowThreshold, double softThreshold,double rotationThreshold,double areaThreshold, double BMDthreshold, 
									String roiChoice,String roiChoiceSt,String rotationChoice,String[] choiceLabels,
									String[] rotationLabels,boolean preventPeeling, boolean allowCleaving, boolean manualRoi,
									boolean manualRotation, double manualAlfa, boolean flipDistribution, 
									boolean guessFlip,boolean guessLarger,boolean stacked,boolean guessStacked, boolean invertGuess,
									int sectorWidth,int divisions,int concentricSector,int concentricDivisions, boolean stOn){
		this.flipHorizontal			= flipHorizontal;
		this.noFiltering			= noFiltering;
		this.sleeveOn				= sleeveOn;
		this.scalingFactor			= scalingFactor;
		this.constant 				= constant;
		this.airThreshold			= airThreshold;
		this.rotationThreshold		= rotationThreshold;
		this.fatThreshold	 		= fatThreshold;
		this.muscleThreshold 		= muscleThreshold;
		this.marrowThreshold 		= marrowThreshold;
		this.areaThreshold 			= areaThreshold;	//For cortical AREA analyses (CoA, SSI, I) + peeling distal pixels
		this.BMDthreshold 			= BMDthreshold;		//For cortical BMD analyses
		this.softThreshold 			= softThreshold;	//Thresholding soft tissues + marrow from bone
		this.boneThreshold 			= areaThreshold;
		this.filterSize				= 3;
		this.softFilterSize			= 7;
		this.sectorWidth			= sectorWidth;
		this.divisions				= divisions;
		this.concentricSector		= concentricSector;
		this.concentricDivisions	= concentricDivisions;
		this.roiChoice				= roiChoice;
		this.roiChoiceSt			= roiChoiceSt;
		this.rotationChoice			= rotationChoice;
		this.choiceLabels			= choiceLabels;
		this.rotationLabels			= rotationLabels;
		this.imageSavePath 			= new String("");
		this.preventPeeling			= preventPeeling;
		this.allowCleaving			= allowCleaving;
		this.manualRoi				= manualRoi;
		this.manualRotation			= manualRotation;
		this.manualAlfa				= manualAlfa;
		this.flipDistribution		= flipDistribution;
		this.guessFlip				= guessFlip;
		this.guessLarger			= guessLarger;
		this.stacked				= stacked;
		this.guessStacked			= guessStacked;
		this.invertGuess			= invertGuess;
		this.sectorWidth			=sectorWidth;
		this.divisions				= divisions;
		this.concentricSector		= concentricSector;
		this.concentricDivisions	= concentricDivisions;
		this.stOn					= stOn;
	}
}
