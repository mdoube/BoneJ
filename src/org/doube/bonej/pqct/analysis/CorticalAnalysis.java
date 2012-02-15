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
import java.lang.Math; //atan2
import org.doube.bonej.pqct.selectroi.*;	//ROI selection..
//import SelectROI.*;
public class CorticalAnalysis{
	public double BMD;
	public double AREA;
	public double MeA;	//Medullary area = ToA-CoA
	public double MaA;	//Marrow area
	public double MaD;	//Marrow density
	public double MaMassD;
	public double StratecMaMassD;
	public double ToA;
	public double ToD;
	public double maxRadiusY;
	public double[] cortexCenter;
	public double SSI;
	public double BSId;
	public double IPo;
	public double IMax;
	public double IMin;
	public double dwIPo;
	public double dwIMax; 
	public double dwIMin;
	public double alfa;
	
	public CorticalAnalysis(SelectROI roi)
	{
		ToA =0;
		ToD = 0;
		MaA =0;
		MaD = 0;
		for (int i =0;i<roi.width*roi.height;i++){
			if (roi.sieve[i] >0){
				ToA +=1;
				ToD +=roi.scaledImage[i];
				if (roi.scaledImage[i] < roi.details.marrowThreshold){ //Marrow analysis
					MaA +=1;
					MaD +=roi.scaledImage[i];
				}
			}
		}
		ToD/=ToA;
		ToA*=roi.pixelSpacing*roi.pixelSpacing;
		MaD/=MaA;
		MaA*=roi.pixelSpacing*roi.pixelSpacing;
		/*Mass density is calculated by converting the BMD to Hounsfield Units, and scaling the HUs to comparable HUs between machines
		**HUs are then scaled to mass density as described in Schneider et al. Phys. Med. Biol. 45 (2000) 459�478.
		*/
		double H;		//Machine comparable Hounsfield Unit
		double mu;		//Med BMD as attenuation coefficient
		double muH2O;	//Water as attenuation coefficient
		muH2O = (0.0-roi.details.constant)/roi.details.scalingFactor;
		mu = (MaD-roi.details.constant)/roi.details.scalingFactor;
		H = mu/muH2O-1.0;	//Equation 6 in Schneider et al. 2000 *1000 omitted
		MaMassD = 1.018+0.893*H;					//Equation 21 in Schneider et al. 2000 *10^-3 omitted
		
		/*Stratec pQCT is calibrated so that fat is 0 vBMD and water is 50 vBMD, Sievanen J Bone Miner Res. 1998 May;13(5):871-82.*/
		muH2O = (50.0-roi.details.constant)/roi.details.scalingFactor;
		mu = (MaD-roi.details.constant)/roi.details.scalingFactor;
		H = mu/muH2O-1.0;	//Equation 6 in Schneider et al. 2000 *1000 omitted
		StratecMaMassD = 1.018+0.893*H;					//Equation 21 in Schneider et al. 2000 *10^-3 omitted
		BSId = ToD*ToD*ToA/100000000.0; //To make it look nicer, we'll use a unit of g^2/cm^4
		BMD = 0;
		AREA = 0;
		cortexCenter = new double[2];
		for (int j = 0;j< roi.cortexRoiI.size();j++){
			BMD+=roi.cortexROI[roi.cortexRoiI.get(j)+roi.cortexRoiJ.get(j)*roi.width];
		}
		BMD/=(double)roi.cortexRoiI.size();
		//Calculate cortical area from 550 threshold...
		for (int j = 0;j< roi.cortexAreaRoiI.size();j++){
			cortexCenter[0]+=(double)roi.cortexAreaRoiI.get(j);
			cortexCenter[1]+=(double)roi.cortexAreaRoiJ.get(j);
		}
		AREA=(double)roi.cortexAreaRoiI.size()*roi.pixelSpacing*roi.pixelSpacing;
		MeA = ToA - AREA;
		cortexCenter[0] /=(double)roi.cortexAreaRoiI.size();
		cortexCenter[1] /=(double)roi.cortexAreaRoiJ.size();
		maxRadiusY = 0; //y for cortical pixels. used for BSI calculations, i.e. density weighted section modulus
		for (int i = 0; i< roi.cortexAreaRoiI.size();i++){
			if (Math.sqrt(((double)roi.cortexAreaRoiI.get(i)-cortexCenter[0])*((double)roi.cortexAreaRoiI.get(i)-cortexCenter[0])
				+((double)roi.cortexAreaRoiJ.get(i)-cortexCenter[1])*((double)roi.cortexAreaRoiJ.get(i)-cortexCenter[1])) > maxRadiusY){
				maxRadiusY = Math.sqrt(((double)roi.cortexAreaRoiI.get(i)-cortexCenter[0])*((double)roi.cortexAreaRoiI.get(i)-cortexCenter[0])
				+((double)roi.cortexAreaRoiJ.get(i)-cortexCenter[1])*((double)roi.cortexAreaRoiJ.get(i)-cortexCenter[1]));
			}
		}
		//Calculate CSMIs and rotation angle to align maximal and minimal bending axes with X and Y axes
		double moment = 0;
		double xmax = 0;
		double ymax = 0;
		double dwmo = 0;
		double dwxmax = 0;
		double dwymax = 0;
		SSI = 0;
		//Calculating cross-sectional moment of inertia in the original image orientation
		for (int i = 0;i<roi.cortexAreaRoiI.size();i++){
			xmax = xmax+((roi.cortexAreaRoiI.get(i)-cortexCenter[0])*roi.pixelSpacing)*((roi.cortexAreaRoiI.get(i)-cortexCenter[0])*roi.pixelSpacing)*roi.pixelSpacing*roi.pixelSpacing;
			ymax = ymax+((roi.cortexAreaRoiJ.get(i)-cortexCenter[1])*roi.pixelSpacing)*((roi.cortexAreaRoiJ.get(i)-cortexCenter[1])*roi.pixelSpacing)*roi.pixelSpacing*roi.pixelSpacing;
			moment = moment+((roi.cortexAreaRoiI.get(i)-cortexCenter[0])*roi.pixelSpacing)*((roi.cortexAreaRoiJ.get(i)-cortexCenter[1])*roi.pixelSpacing)*roi.pixelSpacing*roi.pixelSpacing;
			dwxmax = dwxmax+((roi.cortexAreaRoiI.get(i)-cortexCenter[0])*roi.pixelSpacing/10)*((roi.cortexAreaRoiI.get(i)-cortexCenter[0])*roi.pixelSpacing/10)*(roi.pixelSpacing/10)*(roi.pixelSpacing/10)*(roi.scaledImage[roi.cortexAreaRoiI.get(i)+roi.cortexAreaRoiJ.get(i)*roi.width]);
			dwymax = dwymax+((roi.cortexAreaRoiJ.get(i)-cortexCenter[1])*roi.pixelSpacing/10)*((roi.cortexAreaRoiJ.get(i)-cortexCenter[1])*roi.pixelSpacing/10)*(roi.pixelSpacing/10)*(roi.pixelSpacing/10)*(roi.scaledImage[roi.cortexAreaRoiI.get(i)+roi.cortexAreaRoiJ.get(i)*roi.width]);
			dwmo = dwmo+((roi.cortexAreaRoiI.get(i)-cortexCenter[0])*roi.pixelSpacing/10)*((roi.cortexAreaRoiJ.get(i)-cortexCenter[1])*roi.pixelSpacing/10)*(roi.pixelSpacing/10)*(roi.pixelSpacing/10)*(roi.scaledImage[roi.cortexAreaRoiI.get(i)+roi.cortexAreaRoiJ.get(i)*roi.width]);
			SSI = SSI+(( ((roi.cortexAreaRoiI.get(i)-cortexCenter[0])*roi.pixelSpacing)*((roi.cortexAreaRoiI.get(i)-cortexCenter[0])*roi.pixelSpacing) + ((roi.cortexAreaRoiJ.get(i)-cortexCenter[1])*roi.pixelSpacing)*((roi.cortexAreaRoiJ.get(i)-cortexCenter[1])*roi.pixelSpacing))*roi.pixelSpacing*roi.pixelSpacing*(roi.scaledImage[roi.cortexAreaRoiI.get(i)+roi.cortexAreaRoiJ.get(i)*roi.width]/1200))/(maxRadiusY*roi.pixelSpacing);
		}


		double vali1,vali2,vali3,vali4;		
		IPo = xmax+ymax;
		dwIPo = dwxmax+dwymax;
		//Ipolar caclulated
		//Calculation of Imax and Imin
		//Calculate rotation required to align rotation axes
		alfa = Math.atan(2*moment/(ymax-xmax))/2;
		//Calculate the maximal and minimial cross-sectional moments of inertia
		vali1 = (ymax+xmax)/2+(ymax-xmax)/2*Math.cos(2*(-alfa))-moment*Math.sin(2*(-alfa));
		vali2 =(ymax+xmax)/2-(ymax-xmax)/2*Math.cos(2*(-alfa))+moment*Math.sin(2*(-alfa));
		if (vali1 > vali2){
			IMax = vali1;
			IMin = vali2;
		} else {
			IMax = vali2;
			IMin = vali1;
		}
		//Calculate the maximal and minimal density weighted cross-sectional moments of inertia
		vali3 = (dwymax+dwxmax)/2+(dwymax-dwxmax)/2*Math.cos(2*(-alfa))-dwmo*Math.sin(2*(-alfa));
		vali4 =(dwymax+dwxmax)/2-(dwymax-dwxmax)/2*Math.cos(2*(-alfa))+dwmo*Math.sin(2*(-alfa));
		if (vali3 > vali4){
			dwIMax = vali3;
			dwIMin = vali4;
		} else {
			dwIMax = vali4;
			dwIMin = vali3;
		}			
		
	}

}
