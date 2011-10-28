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

package org.doube.bonej.pqct;

import java.io.*;
import ij.*;
import ij.io.*;
import ij.gui.*;
import ij.plugin.*;
import ij.measure.*;						//Calibration

public class Read_Stratec_File extends ImagePlus implements PlugIn {
	//Global variables
	//Stratec header stuff
	private byte[] PName = new byte[40]; //
	private String PatName;
	private long PatNo;
	private int PatMeasNo;
	private long PatBirth; //
	private long MeasDate; //
	private double VoxelSize; //
	private int PicX0; //
	private int PicY0; //
	private int PicMatrixX; //
	private int PicMatrixY; //
	private byte[] MInfo = new byte[324]; //
	private String MeasInfo;
	private byte[] Dev = new byte[13];
	private String Device;
	private byte[] PID = new byte[13];
	private String PatID;
	private double ObjLen; //
	private String fileName;
	private String properties;
	
	public Read_Stratec_File() { //Constructor
		//this = null;
	}
	
	/**
	 * Opens a Stratec pQCT file
	 * 
	 * @param path full path to the file
	 * @return ImagePlus containing the calibrated pQCT image
	 * @throws IllegalArgumentException if the path is null or the file is not normal
	 */
	public ImagePlus open(String path) throws IllegalArgumentException {
		if (path == null)
			throw new IllegalArgumentException("Path cannot be null");
		File theFile = new File(path);
		if (!theFile.isFile())
			throw new IllegalArgumentException("Path is not a normal file");
		String directory = theFile.getParent()+"/";
		fileName = theFile.getName();
		read(directory);
		fileInfo();
		return this;
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
		fileInfo();
		if (arg.isEmpty() && this.getHeight()>0){
			this.show();
			return;
		}
		if (this.getHeight()<1) return;
	}
	
	private void fileInfo() {
		FileInfo fi = this.getFileInfo(); 
		fi.pixelWidth = VoxelSize;
		fi.pixelHeight = VoxelSize;
		fi.valueUnit = "mm";
		fi.fileName = fileName;
		fi.info		= properties;
		fi.fileType = ij.io.FileInfo.GRAY16_SIGNED;	//
        this.setFileInfo(fi);
	}

	private void read(String directory){
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
			ImagePlus tempImage = NewImage.createShortImage(fileName+" "+Double.toString(VoxelSize), PicMatrixX, PicMatrixY, 1, NewImage.FILL_BLACK);
			this.setImage(tempImage.getImage());
			this.setProcessor(fileName,tempImage.getProcessor());
			//Set ImageJ image properties
			setProperties();
			short[] pixels = (short[]) this.getProcessor().getPixels();
			int min = (int) Math.pow(2,16);
			int max = 0;
			for (int j = 0;j < PicMatrixY; ++j){
				for (int i = 0;i < PicMatrixX; ++i){
					final int offset = 1609+2*(i+j*PicMatrixX);
					int value = (int) ((short) (((fileData[offset+1] & 0xFF)) <<8 | ((short) (fileData[offset] & 0xFF))<<0));

					if (value >=0){
						value = (int) -Math.pow(2,15)+value;
					}else{
						value = (int) Math.pow(2,15)-1+value;
					}
					int tempVal = value & 0xFFFF;
					if (tempVal < min){min =  tempVal;}
					if (tempVal > max){max = tempVal;}					
					pixels[i+j*PicMatrixX] = (short) value;
				}
			}
			this.setDisplayRange( min, max);
			Calibration cal = this.getCalibration();
			double[] coefficients = {-32.768, 0.001};
			cal.setFunction(Calibration.STRAIGHT_LINE, coefficients, "1/cm");
			cal.setUnit("mm");
			cal.pixelWidth = cal.pixelHeight = cal.pixelDepth = VoxelSize;
		}catch (Exception e) {
			IJ.error("Stratec file read failed ", e.getMessage());
			return;
		}
	}
	
	private void readHeader(byte[] fileData){
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
	
	private void setProperties(){
		String[] propertyNames = {"File Name","Pixel Spacing","ObjLen","MeasInfo","Acquisition Date",
									"Device","PatMeasNo","PatNo","Patient's Birth Date","Patient's Name",
									"Patient ID","PicX0","PicY0",
									"Width","Height","Stratec File"};
		String[] propertyValues = {fileName,Double.toString(VoxelSize),Double.toString(ObjLen),MeasInfo,Long.toString(MeasDate),
									Device,Integer.toString(PatMeasNo),Long.toString(PatNo),Long.toString(PatBirth),PatName,
									PatID,Integer.toString(PicX0),Integer.toString(PicY0),
									Integer.toString(PicMatrixX),Integer.toString(PicMatrixY),"1"};
		properties = new String();
		for (int i = 0;i<propertyNames.length;++i){
			properties += propertyNames[i]+": "+propertyValues[i]+"\n";
		}
		this.setProperty("Info", properties);
	}
}
