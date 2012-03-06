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

	Results image creator for ImageJ density distribution analysis plugin
    Copyright (C) 2011 Timo Rantalainen
*/

package org.doube.bonej.pqct.utils;
import ij.ImagePlus;					//Image creation
import ij.IJ;							//Image rotation
import ij.measure.Calibration;			//Image calibration
import ij.process.FloatProcessor;		//Float Images
import ij.process.ImageConverter;		//Convert float to RGB
import java.awt.Color;					//Color class
import java.util.Vector;				//Vector

public class ResultsImage{
	
	/*Get image into which we'll start adding stuff*/
	public static ImagePlus getRGBResultImage(double[] values,int width,int height, String path){
		ImagePlus tempImage = new ImagePlus(path+"Visual results");
		tempImage.setProcessor(new FloatProcessor(width,height,values));
		new ImageConverter(tempImage).convertToRGB();
		return tempImage;
	}
	
	public static ImagePlus addScale(ImagePlus tempImage, double pixelSpacing){
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
	public static ImagePlus addSoftTissueSieve(ImagePlus tempImage, byte[] sieve){
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
	public static ImagePlus addBoneSieve(ImagePlus tempImage, byte[] sieve,double[] scaledImage, double marrowThreshold){
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
	public static ImagePlus addRadii(ImagePlus tempImage,double alfa,double[] marrowCenter,Vector<Integer> pindColor,
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
	
	public static ImagePlus addMarrowCenter(ImagePlus tempImage,double alfa,double[] marrowCenter){
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
	
	public static ImagePlus addRotate(ImagePlus tempImage,double alfa){
		tempImage.getProcessor().setBackgroundValue(0.0);
		IJ.run(tempImage, "Arbitrarily...", "angle=" + alfa + " grid=1 interpolation=Bilinear enlarge");  
		return tempImage;
	}
	
	/*Concentric rings distribution result image*/
	public static ImagePlus addPeriRadii(ImagePlus tempImage,double[] marrowCenter,Vector<Integer> pindColor, double[] R, double[] Theta){
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
}