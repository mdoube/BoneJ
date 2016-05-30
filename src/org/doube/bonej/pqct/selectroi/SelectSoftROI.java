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

package org.doube.bonej.pqct.selectroi;
import java.util.*;	//Vector, Collections
import java.lang.Math; //atan2
import java.awt.*;			//Polygon, Rectangle
import org.doube.bonej.pqct.io.*;	//image data
import ij.*;		//ImagePlus
import ij.gui.*;	//ImagePlus ROI
import ij.text.*; 	//Debugging ...
import ij.process.*;	//Debugging
import org.doube.bonej.pqct.selectroi.liveWireEngine.*;	//LiveWire
import java.util.concurrent.ExecutionException;

@SuppressWarnings(value ={"serial","unchecked"}) //Unchecked for obtaining Vector<Object> as a returnvalue

public class SelectSoftROI extends RoiSelector{
	int radialDivisions = 720;
	public byte[] eroded = null;
	//ImageJ constructor
	public SelectSoftROI(ScaledImageData dataIn,ImageAndAnalysisDetails detailsIn, ImagePlus imp,double boneThreshold,boolean setRoi) throws ExecutionException{
		super(dataIn,detailsIn, imp,boneThreshold,setRoi);
		//Soft tissue analysis
		softSieve = null;
		byte[] softResult = null;
		if (details.stOn){

			/*Get rid of measurement tube used at the UKK institute*/
			byte[] sleeve = null;
			if (details.sleeveOn){
				sleeve = removeSleeve(softScaledImage,sleeve,25.0);
				int removed=0;
				for (int ii =0;ii<width*height;ii++){
					if(sleeve[ii]==1){
						softScaledImage[ii]=minimum;
						++removed;
					}
				}
			}


			/**Ignore data outside manually selected ROI, if manualRoi has been selected*/
			Roi ijROI = imp.getRoi();
			if (ijROI != null && details.manualRoi){	/*Set pixels outside the manually selected ROI to zero*/
			double[] tempScaledImage = Arrays.copyOf(softScaledImage,softScaledImage.length);
				/*Check whether pixel is within ROI, mark with bone threshold*/
				for (int j = 0;j< height;j++){
					for (int i = 0; i < width;i++){
						if (ijROI.contains(i,j)){
						}else{
							softScaledImage[i+j*width] = minimum;
						}
					}
				}
				/*Check whether a polygon can be acquired and include polygon points too*/
				Polygon polygon = ijROI.getPolygon();
				if (polygon != null){
					for (int j = 0;j< polygon.npoints;j++){
						softScaledImage[polygon.xpoints[j]+polygon.ypoints[j]*width] = tempScaledImage[polygon.xpoints[j]+polygon.ypoints[j]*width];
					}
				}
			}



			try{
				Vector<Object> masks = getSieve(softScaledImage,airThreshold,details.roiChoiceSt,details.guessStacked,details.stacked,false,true);
				softSieve						= (byte[]) masks.get(0);
				softResult					 	= (byte[]) masks.get(1);
				Vector<DetectedEdge> stEdges	= (Vector<DetectedEdge>) masks.get(2);

				/*Erode three layers of pixels from the fat sieve to get rid of higher density layer (i.e. skin)
				on top of fat to enable finding muscle border
				*/
				byte[] muscleSieve = Arrays.copyOf(softSieve,softSieve.length);
				double[] muscleImage = Arrays.copyOf(softScaledImage,softScaledImage.length);
				byte[] subCutaneousFat = null;
				//Remove skin by eroding three layers of pixels
				for (int i = 0;i< 3;++i){
					muscleSieve = erode(muscleSieve);
					//Changed 2016/01/08
					/*
					if (i == 0){
						//Add subcut fat sieve... Skin has already been removed by eroding one layer of pixels-> remove muscle later on
						subCutaneousFat = (byte[]) muscleSieve.clone();
					}
					*/
				}
				subCutaneousFat = Arrays.copyOf(muscleSieve,muscleSieve.length);
				/*Remove everything other than the selected limb from the image*/
				for (int i = 0; i<muscleSieve.length;++i){
					if (muscleSieve[i] < 1){
						muscleImage[i] = minimum;
					}
				}
				/*Look for muscle outline*/
				Vector<Object> muscleMasks = getSieve(muscleImage,details.muscleThreshold,"Bigger",details.guessStacked,details.stacked,false,false);
				//muscleSieve		= (byte[]) muscleMasks.get(0);
				Vector<DetectedEdge> muscleEdges = (Vector<DetectedEdge>) muscleMasks.get(2);
				Collections.sort(muscleEdges,Collections.reverseOrder());
				int tempMuscleArea=0;
				muscleSieve = new byte[softSieve.length];
				int areaToAdd=0;
				/*Include areas that contribute more than 1.0% on top of what is already included*/
				while (areaToAdd< muscleEdges.size() && tempMuscleArea*0.01 < muscleEdges.get(areaToAdd).area){
					byte[] tempMuscleSieve = fillSieve(muscleEdges.get(areaToAdd).iit, muscleEdges.get(areaToAdd).jiit,width,height,muscleImage,details.muscleThreshold);
					for (int i = 0; i<tempMuscleSieve.length;++i){
						if (tempMuscleSieve[i] > 0){muscleSieve[i] = tempMuscleSieve[i];}
					}
					tempMuscleArea+=muscleEdges.get(areaToAdd).area;
					areaToAdd++;
				}

				//Visualise muscleSieve
				/*
				ImagePlus tempImage = NewImage.createByteImage("MuscleSieve",width,height,1, NewImage.FILL_BLACK);
				byte[] rPixels = (byte[])tempImage.getProcessor().getPixels();
				for (int i = 0;i<muscleSieve.length;++i){
					for (int c = 0;c<width;++c){
						rPixels[i] = muscleSieve[i];
					}
				}
				tempImage.setDisplayRange(0,10);
				tempImage.show();
				*/
				
				/**Re-segment soft-tissues using livewire based on the muscleSieve
					1) bring rays back from image edges to centre of soft-tissue mask 1 deg apart
					2) use livewire on the 360 edge pixels
					3) rotate livewire init pixels around a few times to get the segmentation to go through subcut/intramuscular fat
				
				*/
				Vector<Object> masks2 = getSieve(softScaledImage,softThreshold,details.roiChoiceSt,details.guessStacked,details.stacked,false,false);
				byte[] boneResult	= (byte[]) masks2.get(1);
				ArrayList<Double> ii = new ArrayList<Double>();
				ArrayList<Double> jj = new ArrayList<Double>();
				double[] softCentre = new double[2];
				
				for (int i = 0;i<width;++i){
					for (int j = 0; j<height;++j){
						if (boneResult[i+j*width] > 0){
						//if (muscleSieve[i+j*width] > 0){
							ii.add((double) i);
							jj.add((double) j);						
							softCentre[0]+=(double) i;
							softCentre[1]+=(double) j;
						}
					}
				}
				
				softCentre[0]/=(double)ii.size();
				softCentre[1]/=(double)jj.size();
				
				
				double maxR = Math.sqrt(Math.pow(max(Math.abs(softCentre[0]-width),softCentre[0]),2d)+Math.pow(max(Math.abs(softCentre[1]-height),softCentre[1]),2d));
				double[] rs = new double[radialDivisions];
				double r,t;
				double[] theta = new double[radialDivisions];
				int[][] edgeCoords = new int[radialDivisions][2];
				//Get the extremes of muscle area with 1 deg increments in polar coordinates
				for (int i = 0;i<radialDivisions;++i){
					r = maxR;
					t = ((double)i)/((double)radialDivisions)*2d*Math.PI;
					while ( Math.round(r*Math.cos(t)+softCentre[0]) < 0 || Math.round(r*Math.sin(t)+softCentre[1]) < 0 ||
							Math.round(r*Math.cos(t)+softCentre[0]) >= width || Math.round(r*Math.sin(t)+softCentre[1]) >= height){
						r-=0.1;
					}
					while ( Math.round(r*Math.cos(t)+softCentre[0]) >= 0 && Math.round(r*Math.sin(t)+softCentre[1]) >= 0 &&
							Math.round(r*Math.cos(t)+softCentre[0]) < width && Math.round(r*Math.sin(t)+softCentre[1]) < height &&
							muscleSieve[(int) (Math.round(r*Math.cos(t)+softCentre[0])+Math.round(r*Math.sin(t)+softCentre[1])*width)] < 1){
						r-=0.1;
					}
					rs[i] = r;
					theta[i] = t;
					edgeCoords[i][0]=(int) (Math.round(rs[i]*Math.cos(theta[i])+softCentre[0]));
					edgeCoords[i][1]=(int) (Math.round(rs[i]*Math.sin(theta[i])+softCentre[1]));
				}

				double[][] pixels = new double[width][height];
				for (int rr = 0;rr<height;++rr){
					for (int c = 0;c<width;++c){
						pixels[c][rr] = (double) muscleSieve[c+rr*width];
					}
				}
				//Arraylists for edge, and livewire seed coordinates
				ArrayList<Integer> edgeii = new ArrayList<Integer>();
				ArrayList<Integer> edgejj = new ArrayList<Integer>();
				ArrayList<Integer> seedii = new ArrayList<Integer>();
				ArrayList<Integer> seedjj = new ArrayList<Integer>();
				
				double lwSteps = 6;
				
				//Create list of seed coordinates
				for (int i = 0;i<radialDivisions;++i){
					edgeii.add(edgeCoords[i][0]);
					edgejj.add(edgeCoords[i][1]);
				}
				//addTrace(tempImage,edgeii,edgejj);
				Vector<Object> tempVO = getLassoEdge(edgeii, edgejj, softCentre, muscleSieve);
				edgeii = (ArrayList<Integer>) tempVO.get(0);
				edgejj = (ArrayList<Integer>) tempVO.get(1);
				
				
				/*
				//Loop livewire 6 times over here
				for (int l = 0; l< ((int) lwSteps);++l){//details.edgeDivisions;++l){
					//Set new seeds
					Vector<Object> tempVO = getSeedPoints(edgeii, edgejj, softCentre,details.edgeDivisions, lwSteps, (double) l);
					seedii = (ArrayList<Integer>) tempVO.get(0);
					seedjj = (ArrayList<Integer>) tempVO.get(1);
				
					edgeii.clear();
					edgejj.clear();
					LiveWireCosts lwc = new LiveWireCosts(pixels);
					int[][] fromSeedToCursor;
					int prevLength = 0;
					for (int i = 0;i<seedii.size()-1; ++i){
						lwc.setSeed(seedii.get(i),seedjj.get(i));
						try{Thread.sleep(10);}catch(Exception e){}
						while ((fromSeedToCursor = lwc.returnPath(seedii.get(i+1),seedjj.get(i+1))) == null){
							try{Thread.sleep(1);}catch(Exception e){}
						}
						for (int iii = 0;iii< fromSeedToCursor.length;++iii){
							edgeii.add((int) fromSeedToCursor[iii][0]);
							edgejj.add((int) fromSeedToCursor[iii][1]);
						}
						if (l == (((int) lwSteps)-1)){
							//IJ.log("Edge Length "+(edgeii.size()-prevLength)+" x "+seedii.get(i)+" y "+seedjj.get(i)+" to x "+seedii.get(i+1)+" y "+seedjj.get(i+1));
						}
						prevLength = edgeii.size();
						
					}
					//Connect the last bit
					lwc.setSeed(seedii.get(seedii.size()-1),seedjj.get(seedjj.size()-1));
					try{Thread.sleep(10);}catch(Exception e){}
					while ((fromSeedToCursor = lwc.returnPath(seedii.get(0),seedjj.get(0))) == null){
						try{Thread.sleep(1);}catch(Exception e){}
					}
					
					
					for (int i = 0;i< fromSeedToCursor.length;++i){
						edgeii.add((int) fromSeedToCursor[i][0]);
						edgejj.add((int) fromSeedToCursor[i][1]);
					}
					if (l == (((int) lwSteps)-1)){
						//IJ.log("Remainder Edge Length "+(edgeii.size()-prevLength)+" x "+seedii.get(seedii.size()-1)+" y "+seedjj.get(seedii.size()-1)+" to x "+seedii.get(0)+" y "+seedjj.get(0));
					}
					
					
				}
				*/
				/*
				IJ.log("Seed points");
				for (int i = 0; i<seedii.size();++i){
					IJ.log(seedii.get(i)+" "+seedjj.get(i));
				}
				
				IJ.log("Edge points");
				for (int i = 0; i<edgeii.size();++i){
					IJ.log(edgeii.get(i)+" "+edgejj.get(i));
				}
				
				
				addTrace(tempImage,edgeii,edgejj);
				*/
				//Fill in muscle mask with inter-muscular fat 
				muscleSieve = getByteMask(width,height,edgeii,edgejj);
				
				/*
				//Visualise the segmentation result
				ImagePlus muscleImage2 = NewImage.createByteImage("muscleImage",width,height,1, NewImage.FILL_BLACK);
				byte[] rPixels3 = (byte[])muscleImage2.getProcessor().getPixels();
				for (int i = 0;i<muscleSieve.length;++i){
					for (int c = 0;c<width;++c){
						rPixels3[i] = muscleSieve[i];
					}
				}
				muscleImage2.setDisplayRange(0,1);
				muscleImage2.show();
				*/
				
				muscleSieve = dilateMuscleMask(muscleSieve,softScaledImage,width,height,muscleThreshold); //Dilate the sieve to include all muscle pixels
				/*Re-segmenting done*/
				
				//Wipe muscle area +3 layer of pixels away from subcut.
				byte[] tempMuscleSieve = Arrays.copyOf(muscleSieve,muscleSieve.length);
				eroded = new byte[softSieve.length];
				//dilate(tempMuscleSieve,(byte)1,(byte)0,(byte)2);
				//dilate(tempMuscleSieve,(byte)1,(byte)0,(byte)2);
				//dilate(tempMuscleSieve,(byte)1,(byte)0,(byte)2);
				for (int i = 0;i<tempMuscleSieve.length;++i){
					if (tempMuscleSieve[i] == 1){subCutaneousFat[i] = 0;}
				}
				
				for (int i = 0;i<softSieve.length;++i){
					if (softSieve[i] ==1 && softScaledImage[i] >= airThreshold && softScaledImage[i] < fatThreshold){
						softSieve[i] =2;	//Fat
					}
					if (muscleSieve[i] ==1 && boneResult[i] ==0 && softScaledImage[i] >= muscleThreshold && softScaledImage[i] < softThreshold){
						softSieve[i] = 3;	//Muscle
					}
					if (muscleSieve[i] ==1 && boneResult[i] ==0 && softScaledImage[i] >= airThreshold && softScaledImage[i] < muscleThreshold){
						softSieve[i] = 4;	//Intra/Intermuscular fat
					}
					if (subCutaneousFat[i] ==1 ){
						softSieve[i] = 5;	//Subcut fat
					}
					if (boneResult[i] ==1 ){
						softSieve[i] = 6;	//Bone & marrow
					}
					if (boneResult[i] ==1 && softScaledImage[i] < fatThreshold){
						softSieve[i] = 7;	//Marrow fat
					}
					if (softSieve[i] > 0 && subCutaneousFat[i] ==0 && tempMuscleSieve[i] ==0){
						eroded[i] = 1;	//Skin eroded pixels
					}
				}
				
				/*
				//Visualise the segmentation result
				ImagePlus softImage = NewImage.createByteImage("SoftSieve",width,height,1, NewImage.FILL_BLACK);
				byte[] rPixels2 = (byte[])softImage.getProcessor().getPixels();
				for (int i = 0;i<softSieve.length;++i){
					for (int c = 0;c<width;++c){
						rPixels2[i] = softSieve[i];
					}
				}
				softImage.setDisplayRange(0,6);
				softImage.show();
				*/
			}catch (ExecutionException err){
				throw err;
			}
		}
	}
	
	/*Debugging*/	
	public void addTrace(ImagePlus tempImage,ArrayList<Integer> edgeii,ArrayList<Integer> edgejj){
		int[] tempii = new int[edgeii.size()];
		int[] tempjj = new int[edgeii.size()];
		for (int i = 0; i<edgeii.size();++i){
			tempii[i] = edgeii.get(i);
			tempjj[i] = edgejj.get(i);
		}
	
		Polygon polygon = new Polygon(tempii,tempjj,tempii.length);
		/*Create the ROI*/
		PolygonRoi roi = new PolygonRoi(polygon,Roi.POLYGON);
		/*Set roi color to differentiate ROIs from each other*/
		/**Get the overlay or create one, if needed*/
		Overlay over;
		if (tempImage.getOverlay() == null){
			over = new Overlay();
			tempImage.setOverlay(over);
		}else{
			over = tempImage.getOverlay();
		}
		
		int colorInd = over.size();
		float[] colors = new float[]{
										0.5f+0.5f*((float) Math.sin(2d*Math.PI*((double)(colorInd-5))/10.0)),	/*R*/
										0.5f+0.5f*((float) Math.cos(2d*Math.PI*((double)colorInd)/10.0)), 	/*G*/
										0.5f+0.5f*((float) Math.sin(2d*Math.PI*((double)colorInd)/10.0))	/*B*/
									};
		roi.setStrokeColor(new Color(colors[0],colors[1],colors[2]));
		/*Add the roi to an overlay, and set the overlay active*/
		tempImage.setRoi(roi,true);
		over.add(roi);
		tempImage.updateAndRepaintWindow();
	}
	
	/*
		Get the line coordinates from origin to target
		Digital Differential Analyzer (DDA) algorithm for line
		http://www.tutorialspoint.com/computer_graphics/line_generation_algorithm.htm
	*/
	private ArrayList<Coordinate> getLine(Coordinate origin, Coordinate target){
		Coordinate difference = target.subtract(origin);
		double steps = difference.maxVal();
		double[] increments = new double[]{difference.ii/ steps,difference.jj/steps};
		ArrayList<Coordinate> coordinates = new ArrayList<Coordinate>();
		coordinates.add(origin);
		for (int i = 0;i<(int) steps;++i){
			coordinates.add(new Coordinate(coordinates.get(coordinates.size()-1).ii+increments[0],coordinates.get(coordinates.size()-1).jj+increments[1]));
		}
		//round
		for (int i = 0;i<coordinates.size();++i){
			coordinates.get(i).ii = Math.round(coordinates.get(i).ii);
			coordinates.get(i).jj = Math.round(coordinates.get(i).jj);
		}
		return coordinates;
	}
	
	public Vector<Object> getLassoEdge(ArrayList<Integer> edgeii, ArrayList<Integer> edgejj, double[] softCentre, byte[] image){
		//IJ.log("Lasso");
		int sectorToConsider = (int) (details.edgeDivisions/360d*((double)radialDivisions));
		/*
		//Visualise the segmentation result
		ImagePlus muscleImage2 = NewImage.createByteImage("muscleImage",width,height,1, NewImage.FILL_BLACK);
		byte[] rPixels3 = (byte[])muscleImage2.getProcessor().getPixels();
		for (int i = 0;i<image.length;++i){
				rPixels3[i] = image[i];
		}
		muscleImage2.setDisplayRange(0,10);
		muscleImage2.show();
		//addTrace(muscleImage2,edgeii,edgejj);
		*/
		
		
		//Pop edge into DetectedRadialEdgeTheta, sort by incrementing radius
		//Start the algorithm from the most distant point from the centre of area
		Vector<DetectedRadialEdgeTheta> radialEdge = new Vector<DetectedRadialEdgeTheta>();
		double theta;
		double r;
		double ii,jj;
		for (int i =0; i<edgeii.size();++i){
			ii = edgeii.get(i)-softCentre[0];
			jj = edgejj.get(i)-softCentre[1];
			radialEdge.add(new DetectedRadialEdgeTheta(edgeii.get(i),edgejj.get(i),Math.atan2(ii,jj),Math.sqrt(Math.pow(ii,2d)+Math.pow(jj,2d)),i));
		}
		//IJ.log("Calculated radii");
		//Get maximal radius, and it's location
		Collections.sort(radialEdge);
		int maxIndex = radialEdge.get(radialEdge.size()-1).index;
		//IJ.log("maxIndex "+maxIndex+" r "+(radialEdge.get(radialEdge.size()-1).radius/2));
		int[] indices = new int[radialDivisions];
		int count = 0;
		for (int i = maxIndex;i<radialDivisions;++i){
			indices[count] = i;
			++count;
		}
		for (int i = 0;i< maxIndex;++i){
			indices[count] = i;
			++count;
		}
		//IJ.log("sorted indices");
		int currentI = 0;
		ArrayList<Integer> boundaryIndices = new ArrayList<Integer>();
		ArrayList<Integer> currentIndices = new ArrayList<Integer>();
		//Loop through the boundary pixels. Look for the furthest point which can be seen
		//From the current point by checking if a line between the points has any roi
		//points in it
		boundaryIndices.add(0);
		while (currentI < (radialDivisions-1)){
			//IJ.log("Start of loop "+currentI);
		  //Go through all point pairs from the furthest, select the first, which
		  //can be seen. Check the next 179 points
		  currentIndices.clear();
		  for (int i = currentI; i<min(currentI+(sectorToConsider-1),radialDivisions); ++i){
			currentIndices.add(i);
		  }
		  
		  int sightIndice = 1;
		  //IJ.log("Checking "+currentI+" currentIndices.size() "+currentIndices.size());
		  for (int i = currentIndices.size()-1; i>0;--i){
			ArrayList<Coordinate> tempCoordinates = getLine(new Coordinate(edgeii.get(indices[currentI]),edgejj.get(indices[currentI])),
new Coordinate(edgeii.get(indices[currentIndices.get(i)]),edgejj.get(indices[currentIndices.get(i)])));
			//IJ.log("Line x "+tempCoordinates.get(t).ii+" y "+tempCoordinates.get(t).jj+" val "+image[(int) (tempCoordinates.get(t).ii+tempCoordinates.get(t).jj*width)]);
			/*
			if (currentI >330){
				for (int t = 0;t<tempCoordinates.size();++t){
					IJ.log("Line x "+tempCoordinates.get(t).ii+" y "+tempCoordinates.get(t).jj+" val "+image[(int) (tempCoordinates.get(t).ii+tempCoordinates.get(t).jj*width)]);
				}
			}
			*/
			int lineOfSight = checkPath(image,tempCoordinates);	//1 = path unblocked, 0 = path blocked
			if (lineOfSight == 1){
				//Line found
				sightIndice = i;
				break;
			}
		  }
		  //IJ.log("currentI "+currentI+" sightIndice "+sightIndice);
		  currentI = currentI+sightIndice;
		  if (currentI > (radialDivisions-1)){
			//Break the loop here
			//boundaryIndices.add(0);	//Connect the ends
			//IJ.log("Break loop "+boundaryIndices.size());
			break;
		  }
		  boundaryIndices.add(currentI);
		  
		}
		
		//Debugging
		/*
		for (int i =0;i<boundaryIndices.size();++i){
			IJ.log("Boundary coordinates i "+i+" x "+(edgeii.get(indices[boundaryIndices.get(i)])/2)+" y "+(edgejj.get(indices[boundaryIndices.get(i)])/2));
		}
		*/	
				
		
		//Create the boundary
		ArrayList<Integer> returnii = new ArrayList<Integer>();
		ArrayList<Integer> returnjj = new ArrayList<Integer>();
		for (int i =1;i<boundaryIndices.size();++i){
			ArrayList<Coordinate> tempCoordinates = getLine(new Coordinate(edgeii.get(indices[boundaryIndices.get(i-1)]),edgejj.get(indices[boundaryIndices.get(i-1)])),new Coordinate(edgeii.get(indices[boundaryIndices.get(i)]),edgejj.get(indices[boundaryIndices.get(i)])));
			for (int j =1;j<tempCoordinates.size();++j){
				returnii.add((int) tempCoordinates.get(j).ii);
				returnjj.add((int) tempCoordinates.get(j).jj);
			}
		}

		//Add the final missing bit
		ArrayList<Coordinate> tempCoordinates = getLine(new Coordinate(edgeii.get(indices[boundaryIndices.get(boundaryIndices.size()-1)]),edgejj.get(indices[boundaryIndices.get(boundaryIndices.size()-1)])),new Coordinate(edgeii.get(indices[boundaryIndices.get(0)]),edgejj.get(indices[boundaryIndices.get(0)])));
			for (int j =1;j<tempCoordinates.size();++j){
				returnii.add((int) tempCoordinates.get(j).ii);
				returnjj.add((int) tempCoordinates.get(j).jj);
			}
		
		//addTrace(muscleImage2,returnii,returnjj);
		
		Vector<Object> returnV = new Vector<Object>();
		returnV.add(returnii);
		returnV.add(returnjj);
		return returnV;
	}
	
	//Helper function to check whether a path is blocked
	public byte checkPath(byte[] image, ArrayList<Coordinate> pathCoordinates){
		int blocked = 0;
		for (int t = pathCoordinates.size()-2;t>0;--t){
		  if (image[(int) (pathCoordinates.get(t).ii+pathCoordinates.get(t).jj*width)] == 1){
			++blocked;
			if (blocked > 2){
				return (byte)  0; //cannot see the point
			}
		  }
		}
		return (byte) 1;	//Can get through the path
	}
	
	//Helper function to get seed points for livewire
	public Vector<Object> getSeedPoints(ArrayList<Integer> edgeii, ArrayList<Integer> edgejj, double[] softCentre,double divisions, double steps, double l){
		ArrayList<Integer> seedii = new ArrayList<Integer>();
		ArrayList<Integer> seedjj = new ArrayList<Integer>();
		
		//Pop edge into DetectedRadialEdge, sort by incrementing radius
		Vector<DetectedRadialEdge> radialEdge = new Vector<DetectedRadialEdge>();
		double theta;
		double r;
		double ii,jj;
		for (int i =0; i<edgeii.size();++i){
			ii = edgeii.get(i)-softCentre[0];
			jj = edgejj.get(i)-softCentre[1];
			radialEdge.add(new DetectedRadialEdge(edgeii.get(i),edgejj.get(i),Math.atan2(ii,jj),Math.sqrt(Math.pow(ii,2d)+Math.pow(jj,2d))));
		}
		Collections.sort(radialEdge);
		
		double tempR,tempR2;
		int maxInd;
		//Get appropriate seeds, select the furthest point from the centre
		for (int i = (int) (l*(((double)radialEdge.size())/divisions)/steps);i<(int) (radialEdge.size()-(((double)radialEdge.size())/divisions));i+=((int)(((double)radialEdge.size())/divisions))){
			
			//Look for the furthest point in this bracket
			maxInd = i;
			tempR = radialEdge.get(i).radius;
			for (int j = i+1;j<i+((int) (((double)edgeii.size())/divisions));++j){
				tempR2 = radialEdge.get(j).radius;
				if (tempR2 > tempR){
					maxInd =j;
					tempR = tempR2;
				}
			}
			seedii.add(radialEdge.get(maxInd).ii);
			seedjj.add(radialEdge.get(maxInd).jj);
		}
		Vector<Object> returnVal = new Vector<Object>();
		returnVal.add(seedii);
		returnVal.add(seedjj);
		return returnVal;		
	}
	
	
	int min(int a,int b){
		return a <= b ? a:b;
	}
	
	double max(double a,double b){
		return a >= b ? a:b;
	}
	
	
	public byte[] dilateMuscleMask(byte[] mask,double[] softScaledImage,int width, int height,double threshold){
		ArrayList<Integer> initialI = new ArrayList<Integer>();
		ArrayList<Integer> initialJ= new ArrayList<Integer>();
		int i,j;
		for (i = 0;i<width;++i){
			for (j = 0;j<height;++j){
				if (mask[i+j*width]==1){
					initialI.add(i);
					initialJ.add(j);
				}
			}
		}
		
		while (initialI.size() >0 && initialI.get(initialI.size()-1) > 0 &&  initialI.get(initialI.size()-1) < width-1
			&& initialJ.get(initialJ.size()-1) > 0 && initialJ.get(initialJ.size()-1) < height-1){
			i =initialI.get(initialI.size()-1);
			j = initialJ.get(initialJ.size()-1);
			initialI.remove( initialI.size()-1);
			initialJ.remove( initialJ.size()-1);

			if (mask[i+j*width] == 0 && softScaledImage[i+j*width] >=threshold){
				mask[i+j*width] = 1;
			}

			if (mask[i-1+j*width] == 0 && softScaledImage[i+j*width] >=threshold) {
				initialI.add(i-1);
				initialJ.add(j);
			}

			if (mask[i+1+j*width] == 0 && softScaledImage[i+j*width] >=threshold) {
				initialI.add(i+1);
				initialJ.add(j);
			}
			
			if (mask[i+(j-1)*width] == 0 && softScaledImage[i+j*width] >=threshold) {
				initialI.add(i);
				initialJ.add(j-1);
			}
			
			if (mask[i+(j+1)*width] == 0 && softScaledImage[i+j*width] >=threshold) {
				initialI.add(i);
				initialJ.add(j+1);
			}

		}
		return mask;
	}
	
	byte[] getByteMask(int width,int height,ArrayList<Integer> edgeii,ArrayList<Integer> edgejj){
		byte[] mask=new byte[width*height];
		for (int i = 0; i<edgeii.size();++i){
			mask[edgeii.get(i)+edgejj.get(i)*width] = (byte) 1;
		}
		int[] fillInitCoords = findMaskFillInit(mask,width,height,edgeii,edgejj);
		//IJ.log("Init x "+fillInitCoords[0]+" y "+fillInitCoords[1]);
		if (fillInitCoords != null){
			return fillMask(fillInitCoords[0],fillInitCoords[1],mask,width,height);
		}else{
			return mask;
		}
	}
	
	byte[] fillMask(int i, int j, byte[] mask, int width,int height){
		ArrayList<Integer> initialI = new ArrayList<Integer>();
		ArrayList<Integer> initialJ= new ArrayList<Integer>();
		initialI.add(i);
		initialJ.add(j);
		while (initialI.size() >0 && initialI.get(initialI.size()-1) > 0 &&  initialI.get(initialI.size()-1) < width-1
			&& initialJ.get(initialJ.size()-1) > 0 && initialJ.get(initialJ.size()-1) < height-1){
			i =initialI.get(initialI.size()-1);
			j = initialJ.get(initialJ.size()-1);
			initialI.remove( initialI.size()-1);
			initialJ.remove( initialJ.size()-1);

			if (mask[i+j*width] == 0 ){
				mask[i+j*width] = 1;
			}

			if (mask[i-1+j*width] == 0) {
				initialI.add(i-1);
				initialJ.add(j);
			}

			if (mask[i+1+j*width] == 0) {
				initialI.add(i+1);
				initialJ.add(j);
			}
			
			if (mask[i+(j-1)*width] == 0) {
				initialI.add(i);
				initialJ.add(j-1);
			}
			
			if (mask[i+(j+1)*width] == 0) {
				initialI.add(i);
				initialJ.add(j+1);
			}

		}
		return mask;
	}
	
	int[] findMaskFillInit(byte[] mask,int width,int height,ArrayList<Integer> edgeii,ArrayList<Integer> edgejj){
		byte[] tempMask = Arrays.copyOf(mask,mask.length);
		tempMask = fillBorder(tempMask,width,height);
		int[] returnCoordinates = new int[2];
		int[] steer = new int[2];
		for (int j = 0; j< edgeii.size()-1; ++j){
			returnCoordinates[0] = edgeii.get(j);
			returnCoordinates[1] = edgejj.get(j);
			double direction = Math.atan2(edgejj.get(j+1)-returnCoordinates[1],edgeii.get(j+1)-returnCoordinates[0]);
			direction+=Math.PI/4.0;
			for (int i = 0; i< 3; ++i){
				
				steer[0] = (int) Math.round(Math.cos(direction));
				steer[1]= (int) Math.round(Math.sin(direction));
				/*Handle OOB*/
				while ((returnCoordinates[0]+steer[0])<0 || (returnCoordinates[0]+steer[0])>=width ||
						(returnCoordinates[1]+steer[1])<0 || (returnCoordinates[1]+steer[1])>=height){
					direction+=Math.PI/4.0;
					steer[0] = (int) Math.round(Math.cos(direction));
					steer[1]= (int) Math.round(Math.sin(direction));
				}
				if (tempMask[returnCoordinates[0]+steer[0]+(returnCoordinates[1]+steer[1])*width] == 0){
					returnCoordinates[0] +=steer[0];
					returnCoordinates[1] +=steer[1];
					return returnCoordinates;
				}
				direction+=Math.PI/4.0;
			}
		}
		return null;
	}
	
	byte[] fillBorder(byte[] mask, int width,int height){
		ArrayList<Integer> initialI = new ArrayList<Integer>();
		ArrayList<Integer> initialJ= new ArrayList<Integer>();
		initialI.add(0);
		initialJ.add(0);
		int i,j;
		while (initialI.size() >0){
			i =initialI.get(initialI.size()-1);
			j = initialJ.get(initialJ.size()-1);
			initialI.remove( initialI.size()-1);
			initialJ.remove( initialJ.size()-1);

			if (mask[i+j*width] == 0 ){
				mask[i+j*width] = 1;
			}

			if (i-1 >= 0 && mask[i-1+j*width] == 0) {
				initialI.add(i-1);
				initialJ.add(j);
			}

			if (i+1 < width && mask[i+1+j*width] == 0) {
				initialI.add(i+1);
				initialJ.add(j);
			}
			
			if (j-1 >= 0 && mask[i+(j-1)*width] == 0) {
				initialI.add(i);
				initialJ.add(j-1);
			}
			
			if (j+1 < width && mask[i+(j+1)*width] == 0) {
				initialI.add(i);
				initialJ.add(j+1);
			}

		}
		return mask;
	}
	
}
