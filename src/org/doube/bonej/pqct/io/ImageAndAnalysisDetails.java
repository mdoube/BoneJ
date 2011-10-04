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
	public String imageSavePath;
	public String roiChoice;
	public String rotationChoice;
	public String[] choiceLabels;
	//ImageJ plugin constructor
	public ImageAndAnalysisDetails(double scalingFactorIn, double constantIn,double areaThresholdIn,double BMDthresholdIn, String roiChoiceIn,String rotationChoice,String[] choiceLabels){
		scalingFactor	= scalingFactorIn;
		constant 		= constantIn;
		airThreshold	= -100;
		fatThreshold 	= 40;
		muscleThreshold = 200;
		marrowThreshold = 300;
		areaThreshold 	= areaThresholdIn;	//For cortical AREA analyses (CoA, SSI, I) + peeling distal pixels
		BMDthreshold 	= BMDthresholdIn;		//For cortical BMD analyses
		softThreshold 	= 300;	//Thresholding soft tissues + marrow from bone
		boneThreshold 	= areaThresholdIn;
		filterSize		= 3;
		softFilterSize	= 7;
		sectorWidth 	= 10;
		roiChoice		= roiChoiceIn;
		this.rotationChoice = rotationChoice;
		this.choiceLabels = choiceLabels;
		imageSavePath 	= new String("");
	}
}