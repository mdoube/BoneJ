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


package org.doube.bonej.pqct.analysis;
import java.util.*;	//Vector, Collections
import org.doube.bonej.pqct.selectroi.*;	//ROI selection..
//import SelectROI.*;
public class MedullaryAnalysis{
	public double medBMD;
	public double medAREA;
	public double kernelBMD;
	public double kernelAREA;
	public double ringBMD;
	public double ringAREA;
	public double[] radialBMD;
	public Vector<Double> marrowDensities;
	public MedullaryAnalysis(SelectROI roi)
	{
		marrowDensities = roi.marrowDensities;
		
		medBMD = 0;	
		medAREA = 0;
		kernelBMD = 0;	
		kernelAREA = 0;
		ringBMD = 0;	
		ringAREA = 0;
      
		for (int i =0;i<roi.width*roi.height;i++){
			if (roi.marrowSieve[i] > 0){
				medAREA++;
				medBMD+=roi.scaledImage[i];
				if (roi.marrowKernelSieve[i]==1){
					kernelAREA++;
					kernelBMD+=roi.scaledImage[i];
				}else{
					ringAREA++;
					ringBMD+=roi.scaledImage[i];
				}
			}
			 
		}
      //System.out.println("Result after for "+medBMD+" "+medAREA+"\n");

		medBMD/=medAREA;
		medAREA*=roi.pixelSpacing*roi.pixelSpacing;
		kernelBMD/=kernelAREA;
		kernelAREA*=roi.pixelSpacing*roi.pixelSpacing;
		ringBMD/=ringAREA;
		ringAREA*=roi.pixelSpacing*roi.pixelSpacing;
		radialBMD = roi.interpolateMarrow(roi.marrowDensities);
      //System.out.println("Result "+medBMD+" "+medAREA+"\n");
	}
}