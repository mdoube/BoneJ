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

	Results writer for ImageJ density distribution analysis plugin
    Copyright (C) 2011 Timo Rantalainen
*/

package org.doube.bonej.pqct.utils;
import org.doube.bonej.pqct.io.*;
import org.doube.bonej.pqct.analysis.*;
import org.doube.bonej.pqct.*;
import ij.text.TextPanel;
import ij.ImagePlus;

public class ResultsWriter{
	public String imageInfo;
	public boolean alphaOn;
	public ResultsWriter(String imageInfo, boolean alphaOn){
		this.imageInfo = imageInfo;
		this.alphaOn = alphaOn;
	}
	
	public void writeHeader(TextPanel textPanel,ImageAndAnalysisDetails imageAndAnalysisDetails){
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
		
		if(imageAndAnalysisDetails.stOn){
			String[] coHeadings = {"MuD [mg/cm�]","MuA [cm�]","LeanMuD [mg/cm�]","LeanMuA [cm�]","IntraFatD [mg/cm�]","IntraFatA [cm�]","FatD [mg/cm�]","FatA [cm�]","SubCutFatD [mg/cm�]","SubCutFatA [cm�]","LimbD [mg/cm�]","LimbA [cm�]","Density weighted fat percentage [%]"};
			for (int i = 0;i<coHeadings.length;++i){
				headings+=coHeadings[i]+"\t";
			}
		}
		
		if(imageAndAnalysisDetails.cOn){
			String[] coHeadings = {"MaMassD [g/cm�]","StratecMaMassD [g/cm�]","MaD [mg/cm�]","MaA [mm�]","CoD [mg/cm�]","CoA [mm�]","SSI [mm�]","ToD [mg/cm�]","ToA[mm�]","MeA [mm�]","BSId[g�/cm4]"};
			for (int i = 0;i<coHeadings.length;++i){
				headings+=coHeadings[i]+"\t";
			}
		}
		if(imageAndAnalysisDetails.mOn){
			for (int i = 0;i<((int) 360/imageAndAnalysisDetails.sectorWidth);++i){
				headings+=i*imageAndAnalysisDetails.sectorWidth+"� - "+((i+1)*imageAndAnalysisDetails.sectorWidth)+"� mineral mass [mg]\t";
			}
		}
		
		if(imageAndAnalysisDetails.conOn){
			for (int i = 0;i<((int) 360/imageAndAnalysisDetails.concentricSector);++i){
				headings+=i*imageAndAnalysisDetails.concentricSector+"� - "+((i+1)*imageAndAnalysisDetails.concentricSector)+"� concentric analysis pericortical radius [mm]\t";
			}
			for (int j = 0;j<imageAndAnalysisDetails.concentricDivisions;++j){
				for (int i = 0;i<((int) 360/imageAndAnalysisDetails.concentricSector);++i){
					headings+="Division "+(j+1)+" sector "+i*imageAndAnalysisDetails.concentricSector+"� - "+((i+1)*imageAndAnalysisDetails.concentricSector)+"� vBMD [mg/cm�]\t";
				}
			}
		}
		
		if(imageAndAnalysisDetails.dOn){
			headings+="Peeled mean vBMD [mg/cm�]\t";
			//Radial distribution
			for (int i =0; i < (int) imageAndAnalysisDetails.divisions; ++i){
				headings+= "Radial division "+i+" vBMD [mg/cm�]\t";
			}
			//Polar distribution
			for (int i = 0;i<((int) 360/imageAndAnalysisDetails.sectorWidth);++i){
				headings+= "Polar sector "+i+" vBMD [mg/cm�]\t";
			}
			
			for (int i = 0;i<((int) 360/imageAndAnalysisDetails.sectorWidth);++i){
				headings+=i*imageAndAnalysisDetails.sectorWidth+"� - "+((i+1)*imageAndAnalysisDetails.sectorWidth)+"� endocortical radius [mm]\t";
			}
			for (int i = 0;i<((int) 360/imageAndAnalysisDetails.sectorWidth);++i){
				headings+=i*imageAndAnalysisDetails.sectorWidth+"� - "+((i+1)*imageAndAnalysisDetails.sectorWidth)+"� pericortical radius [mm]\t";
			}
			//Cortex BMD values			
			for (int i = 0;i<((int) 360/imageAndAnalysisDetails.sectorWidth);++i){
				headings+=i*imageAndAnalysisDetails.sectorWidth+"� - "+((i+1)*imageAndAnalysisDetails.sectorWidth)+"� endocortical vBMD [mg/cm�]\t";
			}
			for (int i = 0;i<((int) 360/imageAndAnalysisDetails.sectorWidth);++i){
				headings+=i*imageAndAnalysisDetails.sectorWidth+"� - "+((i+1)*imageAndAnalysisDetails.sectorWidth)+"� midcortical vBMD [mg/cm�]\t";
			}
			for (int i = 0;i<((int) 360/imageAndAnalysisDetails.sectorWidth);++i){
				headings+=i*imageAndAnalysisDetails.sectorWidth+"� - "+((i+1)*imageAndAnalysisDetails.sectorWidth)+"� pericortical vBMD [mg/cm�]\t";
			}

		}
		textPanel.setColumnHeadings(headings);
	}
	
	public String printResults(String results,ImageAndAnalysisDetails imageAndAnalysisDetails, ImagePlus imp){
		String[] propertyNames = {"File Name","Patient's Name","Patient ID","Patient's Birth Date","Acquisition Date","Pixel Spacing","ObjLen"};
		String[] parameters = {Double.toString(imageAndAnalysisDetails.airThreshold)
								,Double.toString(imageAndAnalysisDetails.fatThreshold),Double.toString(imageAndAnalysisDetails.muscleThreshold)
								,Double.toString(imageAndAnalysisDetails.marrowThreshold)
								,Double.toString(imageAndAnalysisDetails.softThreshold),Double.toString(imageAndAnalysisDetails.rotationThreshold)
								,Double.toString(imageAndAnalysisDetails.areaThreshold),Double.toString(imageAndAnalysisDetails.BMDthreshold)
								,Double.toString(imageAndAnalysisDetails.scalingFactor),Double.toString(imageAndAnalysisDetails.constant)};

		if (imp != null){
			if (Distribution_Analysis.getInfoProperty(imageInfo,"File Name")!= null){
				results+=Distribution_Analysis.getInfoProperty(imageInfo,"File Name")+"\t";
			}else{
				if(imp.getImageStackSize() == 1){
					results+=Distribution_Analysis.getInfoProperty(imageInfo,"Title")+"\t";
				}else{
					results+=imageInfo.substring(0,imageInfo.indexOf("\n"))+"\t";
				}
			}
			for (int i = 1;i<propertyNames.length;++i){
				results+=Distribution_Analysis.getInfoProperty(imageInfo,propertyNames[i])+"\t";
			}
		}
		
		for (int i = 0;i<parameters.length;++i){
			results+=parameters[i]+"\t";
		}

		results += Boolean.toString(imageAndAnalysisDetails.manualRotation)+"\t";
		results += Boolean.toString(imageAndAnalysisDetails.flipDistribution)+"\t";
		results += Boolean.toString(imageAndAnalysisDetails.guessFlip)+"\t";
		results += Boolean.toString(imageAndAnalysisDetails.guessLarger)+"\t";
		results += Boolean.toString(imageAndAnalysisDetails.stacked)+"\t";
		results += Boolean.toString(imageAndAnalysisDetails.invertGuess)+"\t";
		results += Boolean.toString(imageAndAnalysisDetails.allowCleaving)+"\t";
		results += Boolean.toString(imageAndAnalysisDetails.preventPeeling)+"\t";
		results += imageAndAnalysisDetails.roiChoice+"\t";
		results += imageAndAnalysisDetails.rotationChoice+"\t";
		results += Boolean.toString(imageAndAnalysisDetails.flipHorizontal)+"\t";
		results += Boolean.toString(imageAndAnalysisDetails.flipVertical)+"\t";
		return results;
	}
	
	public String printAlfa(String results,DetermineAlfa determineAlfa){
		results += Double.toString(determineAlfa.alfa*180/Math.PI)+"\t";
		results += Double.toString(determineAlfa.rotationCorrection)+"\t";
		results += Double.toString(determineAlfa.distanceBetweenBones)+"\t";
		return results;
	}
	
	public String printSoftTissueResults(String results,SoftTissueAnalysis softTissueAnalysis){
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
	
	public String printCorticalResults(String results,CorticalAnalysis cortAnalysis){
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
		
	public String printMassDistributionResults(String results,MassDistribution massDistribution,ImageAndAnalysisDetails imageAndAnalysisDetails){
		for (int pp = 0;pp<((int) 360/imageAndAnalysisDetails.sectorWidth);pp++){
			results += massDistribution.BMCs[pp]+"\t";
		}
		return results;
	}		
	
	
	public String printConcentricRingResults(String results,ConcentricRingAnalysis concentricRingAnalysis,ImageAndAnalysisDetails imageAndAnalysisDetails){
		for (int i = 0;i<((int) 360/imageAndAnalysisDetails.concentricSector);++i){
			results += concentricRingAnalysis.pericorticalRadii[i]+"\t";
		}
		for (int j = 0;j<imageAndAnalysisDetails.concentricDivisions;++j){
			for (int i = 0;i<((int) 360/imageAndAnalysisDetails.concentricSector);++i){
				results += concentricRingAnalysis.BMDs.get(j)[i]+"\t";
			}
		}
		return results;
	}
	

	
	public String printDistributionResults(String results,DistributionAnalysis DistributionAnalysis,ImageAndAnalysisDetails imageAndAnalysisDetails){
		results+= DistributionAnalysis.peeledBMD+"\t";
		//Radial distribution
		for (int i =0; i < (int) imageAndAnalysisDetails.divisions; ++i){
			results+= DistributionAnalysis.radialDistribution[i]+"\t";
		}
		//Polar distribution
		for (int i = 0;i<((int) 360/imageAndAnalysisDetails.sectorWidth);++i){
			results+= DistributionAnalysis.polarDistribution[i]+"\t";
		}
		
		
		for (int pp = 0;pp<((int) 360/imageAndAnalysisDetails.sectorWidth);++pp){
			results += DistributionAnalysis.endocorticalRadii[pp]+"\t";
		}
		for (int pp = 0;pp<((int) 360/imageAndAnalysisDetails.sectorWidth);++pp){
			results += DistributionAnalysis.pericorticalRadii[pp]+"\t";
		}
		//Cortex BMD values			
		for (int pp = 0;pp<((int) 360/imageAndAnalysisDetails.sectorWidth);++pp){
			results += DistributionAnalysis.endoCorticalBMDs[pp]+"\t";
		}
		for (int pp = 0;pp<((int) 360/imageAndAnalysisDetails.sectorWidth);++pp){
			results += DistributionAnalysis.midCorticalBMDs[pp]+"\t";
		}
		for (int pp = 0;pp<((int) 360/imageAndAnalysisDetails.sectorWidth);++pp){
			results += DistributionAnalysis.periCorticalBMDs[pp]+"\t";
		}
		return results;
	}
}