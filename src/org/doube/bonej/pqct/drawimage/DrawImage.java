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

package org.doube.bonej.pqct.drawimage;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import java.util.Vector;
import java.awt.image.*;
import java.awt.image.WritableRaster;


public class DrawImage extends JPanel implements MouseListener{

	public BufferedImage preparationBuffer;
	public Image imageToDraw;
	public Vector<Integer> coordx;
	public Vector<Integer> coordy;
	public boolean notReady;
	double width;
	double height;
	public DrawImage(){
		setBackground(new Color(0, 0, 0));
		coordx = new Vector<Integer>();
		coordy = new Vector<Integer>();
		notReady = true;
	}

	public void drawImage(BufferedImage bufferedImage) { 
		imageToDraw = Toolkit.getDefaultToolkit().createImage(bufferedImage.getSource());
		imageToDraw= imageToDraw.getScaledInstance(500, -1, Image.SCALE_SMOOTH);
		repaint();
	}
	
	public void drawImage(short[] imageIn,int widthIn, int heightIn,short min, short max) { 
		int[] image = new int[widthIn*heightIn];
      int pixel;
		for (int x = 0; x < widthIn*heightIn;x++) {
         pixel = (int) (((((double) (imageIn[x] -min))/((double)(max-min)))*255.0));
			image[x]= 255<<24 | pixel <<16| pixel <<8| pixel; 
		}
      imageToDraw = createImage(new MemoryImageSource(widthIn,heightIn,image,0,widthIn));
      imageToDraw= imageToDraw.getScaledInstance(500, -1, Image.SCALE_SMOOTH);
		repaint();
		
	}
	
	

	public void drawImage(byte[] imageIn,int widthIn, int heightIn,short min, short max) {
		int[] image = new int[widthIn*heightIn];
      int pixel;
		for (int x = 0; x < widthIn*heightIn;x++) {
         pixel = (int) (((((double) (imageIn[x] -min))/((double)(max-min)))*255.0));
			image[x]= 255<<24 | pixel <<16| pixel <<8| pixel; 
		}
		imageToDraw = createImage(new MemoryImageSource(widthIn,heightIn,image,0,widthIn));
      imageToDraw= imageToDraw.getScaledInstance(500, -1, Image.SCALE_SMOOTH);
   	repaint();
		
	}
	public void drawImage(double[] imageIn,int widthIn, int heightIn,double min, double max) { 
		int[] image = new int[widthIn*heightIn];
      int pixel;
		for (int x = 0; x < widthIn*heightIn;x++) {
         pixel = (int) (((((double) (imageIn[x] -min))/((double)(max-min)))*255.0));
			image[x]= 255<<24 | pixel <<16| pixel <<8| pixel; 
		}
		imageToDraw = createImage(new MemoryImageSource(widthIn,heightIn,image,0,widthIn));
      imageToDraw= imageToDraw.getScaledInstance(500, -1, Image.SCALE_SMOOTH);
   	repaint();
		
	}
	
	public void drawImage(double[] imageIn,int widthIn, int heightIn,double min, double max,double[] marrowCenter,Vector<Integer> pind, double[] R, double[] R2, double[] Theta2) { 
		int[] image = new int[widthIn*heightIn];
      int pixel;
		for (int x = 0; x < widthIn*heightIn;x++) {
         pixel = (int) (((((double) (imageIn[x] -min))/((double)(max-min)))*255.0));
			image[x]= 255<<24 | pixel <<16| pixel <<8| pixel; 
		}
		for(int i = 0; i< 360;i++) {
			image[((int) (marrowCenter[0]+R[pind.get(i)]*Math.cos(Theta2[i])))+  ((int) (marrowCenter[1]+R[pind.get(i)]*Math.sin(Theta2[i])))*widthIn]= 255<<24 | 255 <<16| 0 <<8| 255;
			image[(int) (marrowCenter[0]+R2[pind.get(i)]*Math.cos(Theta2[i]))+ ((int) (marrowCenter[1]+R2[pind.get(i)]*Math.sin(Theta2[i])))*widthIn]=255<<24 | 0 <<16| 255 <<8| 255;
		}
		imageToDraw = createImage(new MemoryImageSource(widthIn,heightIn,image,0,widthIn));
      imageToDraw= imageToDraw.getScaledInstance(500, -1, Image.SCALE_SMOOTH);
   	repaint();
		
	}
	
	/*
		public BufferedImage getMyImage(int[] image,int[] histo,int peak, int[] peaks) {
			
			 
			 Image imageToDraw = createImage(new MemoryImageSource(1000,300,image,0,1000));
//			 imageToDraw= imageToDraw.getScaledInstance(1000, -1, Image.SCALE_SMOOTH);
			 //System.out.println("Piirrettava skaalattu");
			 BufferedImage bufferedImage = (BufferedImage) showFigure.createImage(imageToDraw.getWidth(null), imageToDraw.getHeight(null));
			 //System.out.println("BI tehty "+bufferedImage);
			 Graphics2D gbuf = bufferedImage.createGraphics();
			 //System.out.println("Graphics tehty");
			 gbuf.drawImage(imageToDraw, 0, 0,null);
			 //System.out.println("Piirretty");
			 return bufferedImage;
		}
	
	*/
	
	public void drawImage(int widthIn, int heightIn,int[] histo,int peak, int[] peaks) { 
		int[] image = new int[widthIn*heightIn];
		for (int x = 1; x < 1000;x++) {
			for (int y = 299;y>=299-(int)((((double) histo[x])/((double)peak))*299.0);--y){
				image[x+y*1000]= 255<<24 | 0 <<16| 0 <<8| 0; 
				if (x > peaks[0] && x < peaks[1]){image[x+y*1000]= 255<<24 | 255 <<16| 0 <<8| 0;} 
				if (x > peaks[2] && x < peaks[3]){image[x+y*1000]= 255<<24 | 0 <<16| 255 <<8| 0;} 
				if (x > peaks[4] && x < peaks[5]){image[x+y*1000]= 255<<24 | 0 <<16| 0 <<8| 255;} 
				if (x > peaks[6] && x < peaks[7]){image[x+y*1000]= 255<<24 | 255 <<16| 255 <<8| 0;} 
				if (x > peaks[8] && x < peaks[9]){image[x+y*1000]= 255<<24 | 0 <<16| 255 <<8| 255;} 
				if (x > peaks[10] && x < peaks[11]){image[x+y*1000]= 255<<24 | 255 <<16| 0 <<8| 255;} 
			}
			for (int y = 299-(int)((((double) histo[x])/((double)peak))*299.0)-1;y>=0;--y){
				image[x+y*1000]= 255<<24 | 255 <<16| 255 <<8| 255; 
			}
		}
		
		imageToDraw = createImage(new MemoryImageSource(widthIn,heightIn,image,0,widthIn));
		imageToDraw= imageToDraw.getScaledInstance(500, 300, Image.SCALE_SMOOTH);
		repaint();
		
	}
	
	public void drawImage(int widthIn, int heightIn,int[] histo,int peak, int[] peaks, int histoWidth) { 
		int[] image = new int[widthIn*heightIn];
		for (int x = 1; x < histoWidth;x++) {
			for (int y = 299;y>=299-(int)((((double) histo[x])/((double)peak))*299.0);--y){
				image[x+y*histoWidth]= 255<<24 | 0 <<16| 0 <<8| 0; 
				if (x > peaks[0] && x < peaks[1]){image[x+y*histoWidth]= 255<<24 | 255 <<16| 0 <<8| 0;} 
				if (x > peaks[2] && x < peaks[3]){image[x+y*histoWidth]= 255<<24 | 0 <<16| 255 <<8| 0;} 
				if (x > peaks[4] && x < peaks[5]){image[x+y*histoWidth]= 255<<24 | 0 <<16| 0 <<8| 255;} 
				if (x > peaks[6] && x < peaks[7]){image[x+y*histoWidth]= 255<<24 | 255 <<16| 255 <<8| 0;} 
				if (x > peaks[8] && x < peaks[9]){image[x+y*histoWidth]= 255<<24 | 0 <<16| 255 <<8| 255;} 
				if (x > peaks[10] && x < peaks[11]){image[x+y*histoWidth]= 255<<24 | 255 <<16| 0 <<8| 255;} 
			}
			for (int y = 299-(int)((((double) histo[x])/((double)peak))*299.0)-1;y>=0;--y){
				image[x+y*histoWidth]= 255<<24 | 255 <<16| 255 <<8| 255; 
			}
		}

		imageToDraw = createImage(new MemoryImageSource(widthIn,heightIn,image,0,widthIn));
		imageToDraw= imageToDraw.getScaledInstance(500, 300, Image.SCALE_SMOOTH);
		repaint();
		
	}
	
	public void drawImage(double[] imageIn,int widthIn, int heightIn,double min, double max,int[] coords, int calibWidth, int calibHeight) { 
		int[] image = new int[widthIn*heightIn];
      int pixel;
		for (int x = 0; x < widthIn*heightIn;x++) {
         pixel = (int) (((((double) (imageIn[x] -min))/((double)(max-min)))*255.0));
			image[x]= 255<<24 | pixel <<16| pixel <<8| pixel; 
		}
		for (int j = coords[1]; j<coords[1]+calibHeight;j++){ 
			for (int i = coords[0]; i<coords[0]+calibWidth;i++){
				pixel = (int) (((((double) (imageIn[i+j*widthIn] -min))/((double)(max-min)))*255.0));
				image[i+j*widthIn] = 255<<24 | 128 <<16| pixel <<8| 0; 

			}
		}
		
		imageToDraw = createImage(new MemoryImageSource(widthIn,heightIn,image,0,widthIn));
      imageToDraw= imageToDraw.getScaledInstance(500, -1, Image.SCALE_SMOOTH);
   	repaint();
		
	}
	
	public void drawImage(double[] imageIn,int widthIn, int heightIn,double min, double max,int[] coords, int calibWidth, int calibHeight,int[] histoSieve,int[] peaks) { 
		int[] image = new int[widthIn*heightIn];
      int pixel;
		for (int x = 0; x < widthIn*heightIn;x++) {
         pixel = (int) (((((double) (imageIn[x] -min))/((double)(max-min)))*255.0));
			image[x]= 255<<24 | pixel <<16| pixel <<8| pixel; 
		}
		for (int j = coords[1]; j<coords[1]+calibHeight;j++){
			for (int i = coords[0]; i<coords[0]+calibWidth;i++){
				pixel = (int) (((((double) (imageIn[i+j*widthIn] -min))/((double)(max-min)))*255.0));
				image[i+j*widthIn] = 255<<24 | 128 <<16| pixel <<8| 0; 

			}
		}
		/*Mark calibration phantoms*/
		for (int x = 0; x < widthIn*heightIn;x++){
			
			if (histoSieve[x] > peaks[0] && histoSieve[x] < peaks[1]){image[x]= 255<<24 | 255 <<16| 0 <<8| 0;} 
			if (histoSieve[x] > peaks[2] && histoSieve[x] < peaks[3]){image[x]= 255<<24 | 0 <<16| 255 <<8| 0;} 
			if (histoSieve[x] > peaks[4] && histoSieve[x] < peaks[5]){image[x]= 255<<24 | 0 <<16| 0 <<8| 255;} 
			if (histoSieve[x] > peaks[6] && histoSieve[x] < peaks[7]){image[x]= 255<<24 | 255 <<16| 255 <<8| 0;} 
			if (histoSieve[x] > peaks[8] && histoSieve[x] < peaks[9]){image[x]= 255<<24 | 0 <<16| 255<<8| 255;} 
			if (histoSieve[x] > peaks[10] && histoSieve[x] < peaks[11]){image[x]= 255<<24 | 255 <<16| 0 <<8| 255;} 
		}
		
		
		imageToDraw = createImage(new MemoryImageSource(widthIn,heightIn,image,0,widthIn));
      imageToDraw= imageToDraw.getScaledInstance(500, -1, Image.SCALE_SMOOTH);
   	repaint();
		
	}
	
	

   	public void drawImageSpline(double[] imageIn,int widthIn, int heightIn,double min, double max,Vector<Integer> xx,Vector<Integer> yy) {
		int[] image = new int[widthIn*heightIn];
      int pixel;
		for (int x = 0; x < widthIn*heightIn;x++) {
         pixel = (int) (((((double) (imageIn[x] -min))/((double)(max-min)))*255.0));
			image[x]= 255<<24 | pixel <<16| pixel <<8| pixel; 
		}
      for (int i = 0;i<xx.size();i++){
					image[xx.get(i)+yy.get(i)*widthIn] = 255<<24| 255;
      }
		imageToDraw = createImage(new MemoryImageSource(widthIn,heightIn,image,0,widthIn));
      imageToDraw= imageToDraw.getScaledInstance(500, -1, Image.SCALE_SMOOTH);
   	repaint();

	}

		
	public void paint(Graphics g) {
		g.drawImage(imageToDraw,0,0,null);
	}
	
	
	
	
	//Implement MouseListener
	public void mouseClicked(MouseEvent me) {
			//int screenX = me.getXOnScreen();
			//int screenY = me.getYOnScreen();
			//System.out.println("screen(X,Y) = " + screenX + "," + screenY);

	}
	public void mousePressed(MouseEvent me) {
//            maybeShowPopup(e);
	}
	public void mouseExited(MouseEvent me) {
//            maybeShowPopup(e);
	}
	public void mouseEntered(MouseEvent me) {
//            maybeShowPopup(e);
	}
	public void mouseReleased(MouseEvent me) {
//          maybeShowPopup(e);
		if (me.getButton() == MouseEvent.BUTTON1){
			coordx.add((int) ((double) me.getX()/500.0*width));				
			coordy.add((int) ((double) me.getY()/500.0*height));
			preparationBuffer.setRGB(coordx.lastElement(),coordy.lastElement(),(byte) 255);
			imageToDraw= preparationBuffer.getScaledInstance(500, -1, Image.SCALE_SMOOTH);
			repaint();
			System.out.println("screen(X,Y) = " + coordx.lastElement() + "," + coordy.lastElement());
		}else{
			notReady = false;
		}
	}
	

} 
	