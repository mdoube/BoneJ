package org.doube.bonej.pqct;

import java.io.*;
import ij.*;
import ij.io.*;
import ij.gui.*;
import ij.plugin.*;

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

	ImageJ Stratec file reader plugin
    Copyright (C) 2011 Timo Rantalainen
*/

public class Read_Stratec_File extends ImagePlus implements PlugIn {
	//Global variables
	ImagePlus img;
	//Stratec header stuff
	public byte[] PName = new byte[40]; //
	public String PatName;
	public long PatNo;
	public int PatMeasNo;
	public long PatBirth; //
	public long MeasDate; //
	public double VoxelSize; //
	public int PicX0; //
	public int PicY0; //
	public int PicMatrixX; //
	public int PicMatrixY; //
	public byte[] MInfo = new byte[324]; //
	public String MeasInfo;
	public byte[] Dev = new byte[13];
	public String Device;
	public byte[] PID = new byte[13];
	public String PatID;
	int offset;
	public double ObjLen; //
	public short[] data;
	public short min,max;
	
	
	public Read_Stratec_File() { //Constructor
		img = null;
	}
	
	//Overriding the abstract runnable run method. Apparently plugins run in threads
	public void run(String arg) { 
		String fileName;
		String directory;
		File file;
		//ij.IJ.showStatus("In Stratec Reader "+arg);
		//try{Thread.sleep(5000);}catch (Exception err){}
		
		if (!arg.isEmpty()){//Called by HandleExtraFileTypes
			File theFile = new File(arg);
			directory = theFile.getParent()+"/";
			fileName = theFile.getName();
			//ij.IJ.showStatus(directory+" "+fileName);
			//try{Thread.sleep(500);}catch (Exception err){}
		
		}else{//select file manually
			OpenDialog od = new OpenDialog("Select stratec image (I*.M*)", arg);
			directory = od.getDirectory();
			fileName = od.getFileName();
		}
		if (fileName==null) return;
		read(directory,fileName);

		if (img==null) return;
		img.show();
	}
	
	
	
	public void read(String directory, String fileName){
		File fileIn = new File(directory+fileName);
		long fileLength = fileIn.length();
		try{BufferedInputStream inputStream = new BufferedInputStream(new FileInputStream(fileIn));
			DataInputStream dataInputStream = new DataInputStream(inputStream);
			byte[] fileData = new byte[(int) fileLength];		//Allocate memory for reading the file into memory
			dataInputStream.read(fileData,0,(int) fileLength);		//Read the data to memory
			dataInputStream.close();	//Close the file after reading
			//Read some data from the file Header
			readHeader(fileData);
			//Create ImageJ image
			img = NewImage.createFloatImage(fileName+" "+Double.toString(VoxelSize), PicMatrixX, PicMatrixY, 1, NewImage.FILL_BLACK);
			//Set ImageJ image properties
			setProperties(img,fileName);
			float[] pixels = (float[]) img.getProcessor().getPixels();
			float min = (float) Math.pow(2,15)-1;
			float max = (float) -Math.pow(2,15);
			for (int j = 0;j < PicMatrixY; ++j){
				for (int i = 0;i < PicMatrixX; ++i){
					offset = 1609+2*(i+j*PicMatrixX);
					float value =(float) ((short) (((fileData[offset+1] & 0xFF)) <<8 | ((short) (fileData[offset] & 0xFF))<<0));
					if (value < min){min = value;}
					if (value > max){max = value;}
					pixels[i+j*PicMatrixX] = value;
				}
			}
			//System.out.println(min +" "+max+" "+voxelSize+" "+PicMatrixX+" "+PicMatrixY);
			img.setDisplayRange((double) min,(double) max);
			//Properties props = img.getProperties();
			//System.out.println(props);
		 }catch (Exception e) {
			IJ.error("Stratec file read failed ", e.getMessage());
			return;
		}
	}
	
	public void readHeader(byte[] fileData){
		byte stringLength;
		offset =12;
		VoxelSize = Double.longBitsToDouble((long) ( ((long) (fileData[offset+7] & 0xFF))<<56 | ((long) (fileData[offset+6] & 0xFF))<<48 | ((long) (fileData[offset+5] & 0xFF))<<40 | ((long) (fileData[offset+4] & 0xFF))<<32 | ((long) (fileData[offset+3] & 0xFF))<<24 | ((long) (fileData[offset+2] & 0xFF))<<16 | ((long) (fileData[offset+1] & 0xFF))<<8 | ((long) (fileData[offset+0] & 0xFF))<<0));

		offset = 318;
		ObjLen = Double.longBitsToDouble((long) ( ((long) (fileData[offset+7] & 0xFF))<<56 | ((long) (fileData[offset+6] & 0xFF))<<48 | ((long) (fileData[offset+5] & 0xFF))<<40 | ((long) (fileData[offset+4] & 0xFF))<<32 | ((long) (fileData[offset+3] & 0xFF))<<24 | ((long) (fileData[offset+2] & 0xFF))<<16 | ((long) (fileData[offset+1] & 0xFF))<<8 | ((long) (fileData[offset+0] & 0xFF))<<0));

		offset = 662;
		stringLength = fileData[offset];
		for (int i=0;i<stringLength;++i){MInfo[i] = fileData[offset+1+i];}
		MeasInfo = new String(MInfo);

		offset = 986;
		MeasDate =  (long) (((long) (fileData[offset+3] & 0xFF)) <<24 | ((long) (fileData[offset+2] & 0xFF))<<16 | ((long) (fileData[offset+1] & 0xFF)) <<8 | ((long) (fileData[offset+0] & 0xFF)));

		offset = 1050;
		stringLength = fileData[offset];
		for (int i=0;i<stringLength;++i){Dev[i] = fileData[offset+1+i];}
		Device = new String(Dev);		//DeviceType

		offset = 1085;
		PatMeasNo = ((int) ((int) (fileData[offset+1] & 0xFF)) <<8 | ((int) (fileData[offset+0] & 0xFF)));

		offset = 1087;
		PatNo =  (long) (((long) (fileData[offset+3] & 0xFF)) <<24 | ((long) (fileData[offset+2] & 0xFF))<<16 | ((long) (fileData[offset+1] & 0xFF)) <<8 | ((long) (fileData[offset+0] & 0xFF)));

		offset = 1091;
		PatBirth =  (long) (((long) (fileData[offset+3] & 0xFF)) <<24 | ((long) (fileData[offset+2] & 0xFF))<<16 | ((long) (fileData[offset+1] & 0xFF)) <<8 | ((long) (fileData[offset+0] & 0xFF)));
		
		offset = 1099;
		stringLength = fileData[offset];
		for (int i=0;i<stringLength;++i){PName[i] = fileData[offset+1+i];}
		PatName = new String(PName);
		
		offset = 1282;
		stringLength = fileData[offset];
		for (int i=0;i<stringLength;++i){PID[i] = fileData[offset+1+i];}
		PatID = new String(PID);

		offset = 1525;
		PicX0 = ((int) ((int) (fileData[offset+1] & 0xFF)) <<8 | ((int) (fileData[offset+0] & 0xFF)));
		
		offset = 1527;
		PicY0 = ((int) ((int) (fileData[offset+1] & 0xFF)) <<8 | ((int) (fileData[offset+0] & 0xFF)));
		
		offset = 1529;
		PicMatrixX = ((int) ((int) (fileData[offset+1] & 0xFF)) <<8 | ((int) (fileData[offset+0] & 0xFF)));
		
		offset = 1531;
		PicMatrixY = ((int) ((int) (fileData[offset+1] & 0xFF)) <<8 | ((int) (fileData[offset+0] & 0xFF)));
	}
	
	public void setProperties(ImagePlus img,String FileName){
		img.setProperty(new String("FileName"),FileName);
		img.setProperty(new String("VoxelSize"),VoxelSize);
		img.setProperty(new String("ObjLen"),ObjLen);
		img.setProperty(new String("MeasInfo"),MeasInfo);
		img.setProperty(new String("MeasDate"),MeasDate);
		img.setProperty(new String("Device"),Device);
		img.setProperty(new String("PatMeasNo"),PatMeasNo);
		img.setProperty(new String("PatNo"),PatNo);
		img.setProperty(new String("PatBirth"),PatBirth);
		img.setProperty(new String("PatName"),PatName);
		img.setProperty(new String("PatID"),PatID);
		img.setProperty(new String("PicX0"),PicX0);
		img.setProperty(new String("PicY0"),PicY0);
		img.setProperty(new String("PicMatrixX"),PicMatrixX);
		img.setProperty(new String("PicMatrixY"),PicMatrixY);
	}

}
