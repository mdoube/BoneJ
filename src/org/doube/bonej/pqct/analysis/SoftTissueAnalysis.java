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

package Analysis;
import SelectRoi.*;	//ROI selection..
public class SoftTissueAnalysis{
	public double stBMD;
	public double stAREA;
	public double subBMD;
	public double subAREA;
	public double muscleBMD;
	public double muscleAREA;
	public double legAREA;
	
	public double fatPercentage;
	public double fatPercentageWeighted;
	public double weightedLegArea;
	public double weightedFatArea;
	public SoftTissueAnalysis(SelectROI roi){
	
		//Analyze soft tissue BMDs and Areas here...
		stBMD=0;
		stAREA=0;
		subBMD=0;
		subAREA=0;
		muscleBMD=0;
		muscleAREA=0;
		legAREA=0;
		weightedLegArea =0;
		weightedFatArea =0;
		for (int i =0;i<roi.width*roi.height;i++){
			if (roi.stSieve[i]==1 && roi.softScaledImage[i] >= roi.airThreshold && roi.softScaledImage[i] < roi.fatThreshold){
				subBMD+=roi.softScaledImage[i];
				subAREA+=1;
				weightedFatArea += roi.softScaledImage[i]+1000.0;
				stBMD+=roi.softScaledImage[i];
				stAREA+=1;
			}
			if (roi.muscleSieve[i]==1 && roi.softScaledImage[i] >= roi.fatThreshold && roi.softScaledImage[i] < roi.muscleThreshold){
				muscleBMD+=roi.softScaledImage[i];
				muscleAREA+=1;
				stBMD+=roi.softScaledImage[i];
				stAREA+=1;
			}
			if (roi.legSieve[i] == 1){
				legAREA+=1;
				weightedLegArea += roi.softScaledImage[i]+1000.0;
				
			}
		}
		subBMD/=subAREA;
		muscleBMD/=muscleAREA;
		stBMD/=stAREA;
		subAREA*=roi.pixelSpacing*roi.pixelSpacing;
		muscleAREA*=roi.pixelSpacing*roi.pixelSpacing;
		stAREA*=roi.pixelSpacing*roi.pixelSpacing;
		legAREA*=roi.pixelSpacing*roi.pixelSpacing;
		fatPercentage = 100*(subAREA/legAREA);
		fatPercentageWeighted= 100*(weightedFatArea/weightedLegArea);
	}
}