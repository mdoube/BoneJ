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
	public double ObjLen; //
	public short[] data;
	public short min,max;
	public String fileName;
	public String properties;
	
	public Read_Stratec_File() { //Constructor
		//this = null;
	}
	
	//Overriding the abstract runnable run method. Apparently plugins run in threads
	public void run(String arg) { 
		String directory;
	
		if (!arg.isEmpty()){//Called by HandleExtraFileTypes
			File theFile = new File(arg);
			directory = theFile.getParent()+"/";
			fileName = theFile.getName();
		}else{//select file manually
			OpenDialog od = new OpenDialog("Select stratec image (I*.M*)", arg);
			directory = od.getDirectory();
			fileName = od.getFileName();
		}
		if (fileName==null) return;
		read(directory);
		FileInfo fi = this.getFileInfo(); 
		fi.pixelWidth = VoxelSize;
		fi.pixelHeight = VoxelSize;
		fi.valueUnit = "mm";
		fi.fileName = fileName;
		fi.info		= properties;
		fi.fileType = ij.io.FileInfo.GRAY16_SIGNED;	//
        this.setFileInfo(fi);  
		if (this.getHeight()<1) return;
	}
	
	public void read(String directory){
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
			ImagePlus tempImage = NewImage.createFloatImage(fileName+" "+Double.toString(VoxelSize), PicMatrixX, PicMatrixY, 1, NewImage.FILL_BLACK);
			this.setImage(tempImage.getImage());
			this.setProcessor(fileName+" "+Double.toString(VoxelSize),tempImage.getProcessor());
			//Set ImageJ image properties
			setProperties();
			float[] pixels = (float[]) this.getProcessor().getPixels();
			float min = (float) Math.pow(2,15)-1;
			float max = (float) -Math.pow(2,15);
			for (int j = 0;j < PicMatrixY; ++j){
				for (int i = 0;i < PicMatrixX; ++i){
					final int offset = 1609+2*(i+j*PicMatrixX);
					float value =(float) ((short) (((fileData[offset+1] & 0xFF)) <<8 | ((short) (fileData[offset] & 0xFF))<<0));
					if (value < min){min = value;}
					if (value > max){max = value;}
					pixels[i+j*PicMatrixX] = value;
				}
			}
			this.setDisplayRange((double) min,(double) max);
		}catch (Exception e) {
			IJ.error("Stratec file read failed ", e.getMessage());
			return;
		}
	}
	
	public void readHeader(byte[] fileData){
		byte stringLength;
		int offset =12;
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
	
	public void setProperties(){
		String[] propertyNames = {"File Name","Pixel Spacing","ObjLen","MeasInfo","Acquisition Date",
									"Device","PatMeasNo","PatNo","Patient's Birth Date","Patient's Name",
									"Patient ID","PicX0","PicY0",
									"Width","Height"};
		String[] propertyValues = {fileName,Double.toString(VoxelSize),Double.toString(ObjLen),MeasInfo,Long.toString(MeasDate),
									Device,Integer.toString(PatMeasNo),Long.toString(PatNo),Long.toString(PatBirth),PatName,
									PatID,Integer.toString(PicX0),Integer.toString(PicY0),
									Integer.toString(PicMatrixX),Integer.toString(PicMatrixY)};
		properties = new String();
		for (int i = 0;i<propertyNames.length;++i){
			properties += propertyNames[i]+": "+propertyValues[i]+"\n";
		}
		this.setProperty("Info", properties);
	}
}
