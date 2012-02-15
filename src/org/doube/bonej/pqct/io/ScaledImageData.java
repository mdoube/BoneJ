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
import java.util.*;	//Vector, Collections

public class ScaledImageData{
	public double[] scaledImage;
	public double[] softScaledImage;
	public double[] justROI;
	public double minimum;
	public double maximum;
	public int width;
	public int height;
	public double pixelSpacing;
	int filterSize;
	//Constructor
	public ScaledImageData(int[] data, int widthIn, int heightIn, double VoxelSize, double scalingFactor, double constant, int filterSize,boolean flipHorizontal, boolean flipVertical,boolean noFiltering){
		height = heightIn;
		width = widthIn;
		pixelSpacing = VoxelSize;
		filterSize = 3;		//filterSize x filterSize median filter will be used
		double[] unFiltered = new double[width*height];

		for (int t = 0;t<width*height;t++){	//Scale the image
			unFiltered[t] = ((double) data[t])*scalingFactor+constant;
		}
		/*Get the min and max values*/
		double[] tempSort = (double[]) unFiltered.clone();
		Arrays.sort(tempSort);
		minimum = tempSort[0];
		maximum = tempSort[tempSort.length-1];
		softScaledImage = medianFilter(unFiltered,width,height,7); //Median filter data
		if (noFiltering){
			scaledImage = (double[]) unFiltered.clone();
		}else{		
			scaledImage = medianFilter(unFiltered,width,height,filterSize); //Median filter data
		}
		
		if (flipHorizontal){//Flip the image around the horizontal axis...
			flipHorizontally();
		}
		if (flipVertical){//Flip the image around the horizontal axis...
			flipVertically();
		}
	}
	
	public void flipHorizontally(){
		double[] temp = (double[]) scaledImage.clone();
		double[] temp2 = (double[]) softScaledImage.clone();
		for (int j = 0;j<height;++j){
			for (int i = 0;i<width;++i){
				scaledImage[i+(height-1-j)*width] = temp[i+j*width];
				softScaledImage[i+(height-1-j)*width] = temp2[i+j*width];
			}
		}
	}
	public void flipVertically(){
		double[] temp = (double[]) scaledImage.clone();
		double[] temp2 = (double[]) softScaledImage.clone();
		for (int j = 0;j<height;++j){
			for (int i = 0;i<width;++i){
				scaledImage[(width-1-i)+j*width] = temp[i+j*width];
				softScaledImage[(width-1-i)+j*width] = temp2[i+j*width];
			}
		}
	}

	//Meadian filter
	public double[] medianFilter(double[] data, int width, int height,int filterSize){
		double[] filtered = new double[width*height];
		double[] toMedian = new double[filterSize*filterSize];
		for (int i = 0; i< filtered.length; ++i) {filtered[i] = minimum;}
		/*Fill filtered with min value to get the frame from messing up with edge detection*/
		
		int noGo = (int)Math.floor(((double)filterSize)/2.0);
		int median = (int)Math.floor(((double)(filterSize*filterSize))/2.0);	//would be ceil, but indexing begins from 0!
		int rowTotal,colTotal;
		for (int row = noGo; row < height-noGo; row++) {
			for (int col = noGo; col < width-noGo; col++) {
				int newPixel = 0;
				for (int rowOffset=-noGo; rowOffset<=noGo; rowOffset++) {
					for (int colOffset=-noGo; colOffset<=noGo; colOffset++) {
						rowTotal = row + rowOffset;
						colTotal = col + colOffset;
						toMedian[newPixel]= data[rowTotal*width+colTotal];
						newPixel++;
					}
				}
				Arrays.sort(toMedian);
				filtered[row*width+col] = toMedian[median];
			}
		}
		
		
		return filtered;
	}
	
}
