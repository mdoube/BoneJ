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
import java.io.*;

public class ReadStratecFile {
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
		private byte[] bytesFromFile = new byte[2];
		int[] bits = new int[2];
		public double ObjLen; //
		public short[] data;
		public short min,max;
	public ReadStratecFile(DataInputStream inFile){
		try{inFile.mark(1609);}catch (Exception err){System.err.println("Error Mark: "+err.getMessage());}
		try{inFile.skip(12);}catch (Exception err){System.err.println("Error: "+err.getMessage());} //Offsets and value representations obtained from Stratec
		try{VoxelSize = Double.longBitsToDouble(Long.reverseBytes(inFile.readLong()));}catch (Exception err){System.err.println("Error: "+err.getMessage());}
		try{inFile.reset();}catch (Exception err){System.err.println("Error Reset1: "+err.getMessage());}
		try{inFile.skip(318);}catch (Exception err){System.err.println("Error: "+err.getMessage());} //Offsets and value representations obtained from Stratec
		try{ObjLen = Double.longBitsToDouble(Long.reverseBytes(inFile.readLong()));}catch (Exception err){System.err.println("Error: "+err.getMessage());}		
		try{inFile.reset();}catch (Exception err){System.err.println("Error Reset2: "+err.getMessage());}
		try{inFile.skip(662);}catch (Exception err){System.err.println("Error: "+err.getMessage());} //Offsets and value representations obtained from Stratec
		
		byte[] stringLength = new byte[1];
		try{inFile.read(stringLength);}catch (Exception err){System.err.println("Error: "+err.getMessage());}		//Check string length
		try{inFile.read(MInfo,0,stringLength[0]);}catch (Exception err){System.err.println("Error: "+err.getMessage());}
		MeasInfo = new String(MInfo);
		//inFile.read((char*) &mInfo,324);
		try{inFile.reset();}catch (Exception err){System.err.println("Error Reset3: "+err.getMessage());}
		try{inFile.skip(986);}catch (Exception err){System.err.println("Error: "+err.getMessage());} //Offsets and value representations obtained from Stratec
		try{MeasDate = Integer.reverseBytes(inFile.readInt());}catch (Exception err){System.err.println("Error: "+err.getMessage());}
		//MeasDate = (short) bytesFromFile[0];
		try{inFile.reset();}catch (Exception err){System.err.println("Error: "+err.getMessage());}
		try{inFile.skip(1050);}catch (Exception err){System.err.println("Error: "+err.getMessage());}
		try{inFile.read(stringLength);}catch (Exception err){System.err.println("Error: "+err.getMessage());}		//Tarkistetaan deviceStringin stringLength...
		try{inFile.read(Dev,0,stringLength[0]);}catch (Exception err){System.err.println("Error: "+err.getMessage());}
		Device = new String(Dev);		//DeviceType

		try{inFile.reset();}catch (Exception err){System.err.println("Error: "+err.getMessage());}
		try{inFile.skip(1085);}catch (Exception err){System.err.println("Error: "+err.getMessage());} //Offsets and value representations obtained from Stratec
		try{inFile.read(bytesFromFile,0,2);}catch (Exception err){System.err.println("Error: "+err.getMessage());}
		bits[0] = 0;
		bits[1] = 0;
		bits[0] = (0xFF & ((int)bytesFromFile[1]));
		bits[1] = (0xFF & ((int)bytesFromFile[0]));
		PatMeasNo = (int) (bits[0] << 8 | bits[1]);	//BigEndian encoding...
		try{inFile.reset();}catch (Exception err){System.err.println("Error: "+err.getMessage());}
		try{inFile.skip(1087);}catch (Exception err){System.err.println("Error: "+err.getMessage());} //Offsets and value representations obtained from Stratec
		try{PatNo = Integer.reverseBytes(inFile.readInt());}catch (Exception err){System.err.println("Error: "+err.getMessage());}
		
		try{inFile.reset();}catch (Exception err){System.err.println("Error: "+err.getMessage());}
		try{inFile.skip(1091);}catch (Exception err){System.err.println("Error: "+err.getMessage());} //Offsets and value representations obtained from Stratec
		try{PatBirth = Integer.reverseBytes(inFile.readInt());}catch (Exception err){System.err.println("Error: "+err.getMessage());}
		try{inFile.reset();}catch (Exception err){System.err.println("Error: "+err.getMessage());}
		try{inFile.skip(1099);}catch (Exception err){System.err.println("Error: "+err.getMessage());} //Offsets and value representations obtained from Stratec
		try{inFile.read(stringLength);}catch (Exception err){System.err.println("Error: "+err.getMessage());}		//Check string length
		try{inFile.read(PName,0,stringLength[0]);}catch (Exception err){System.err.println("Error: "+err.getMessage());}
		PatName = new String(PName);
		try{inFile.reset();}catch (Exception err){System.err.println("Error: "+err.getMessage());}
		try{inFile.skip(1282);}catch (Exception err){System.err.println("Error: "+err.getMessage());} //Offsets and value representations obtained from Stratec
		try{inFile.read(stringLength);}catch (Exception err){System.err.println("Error: "+err.getMessage());}		//Check string length
		try{inFile.read(PID,0,stringLength[0]);}catch (Exception err){System.err.println("Error: "+err.getMessage());}
		PatID = new String(PID);
		try{inFile.reset();}catch (Exception err){System.err.println("Error: "+err.getMessage());}
		try{inFile.skip(1525);}catch (Exception err){System.err.println("Error: "+err.getMessage());} //Offsets and value representations obtained from Stratec
		try{inFile.read(bytesFromFile,0,2);}catch (Exception err){System.err.println("Error: "+err.getMessage());}
		bits[0] = 0;
		bits[1] = 0;
		bits[0] = (0xFF & ((int)bytesFromFile[1]));
		bits[1] = (0xFF & ((int)bytesFromFile[0]));
		PicX0 = (int) (bits[0] << 8 | bits[1]);	//BigEndian encoding...
		try{inFile.read(bytesFromFile,0,2);}catch (Exception err){System.err.println("Error: "+err.getMessage());}
		bits[0] = 0;
		bits[1] = 0;
		bits[0] = (0xFF & ((int)bytesFromFile[1]));
		bits[1] = (0xFF & ((int)bytesFromFile[0]));
		PicY0 = (int) (bits[0] << 8 | bits[1]);	//BigEndian encoding...
		try{inFile.read(bytesFromFile,0,2);}catch (Exception err){System.err.println("Error: "+err.getMessage());}
		bits[0] = 0;
		bits[1] = 0;
		bits[0] = (0xFF & ((int)bytesFromFile[1]));
		bits[1] = (0xFF & ((int)bytesFromFile[0]));
		PicMatrixX = (int) (bits[0] << 8 | bits[1]);	//BigEndian encoding...
		try{inFile.read(bytesFromFile,0,2);}catch (Exception err){System.err.println("Error: "+err.getMessage());}
		bits[0] = 0;
		bits[1] = 0;
		bits[0] = (0xFF & ((int)bytesFromFile[1]));
		bits[1] = (0xFF & ((int)bytesFromFile[0]));
		PicMatrixY = (int) (bits[0] << 8 | bits[1]);	//BigEndian encoding...
		try{inFile.reset();}catch (Exception err){System.err.println("Error: "+err.getMessage());}
		try{inFile.skip(1609);}catch (Exception err){System.err.println("Error: "+err.getMessage());} //Offsets and value representations obtained from Stratec
		//Read the actual DATA
		data = new short[PicMatrixX*PicMatrixY];
		min = 0;
		max = 0;
		for (int i = 0;i<PicMatrixX*PicMatrixY;i++){
			try{data[i] = Short.reverseBytes(inFile.readShort());}catch (Exception err){System.err.println("Error: "+err.getMessage());}
			if (data[i] < min){min = data[i];}
			if (data[i] > max){max = data[i];}
		}
		try{inFile.close();}catch (Exception err){System.err.println("Error: "+err.getMessage());}
	}
	public ReadStratecFile(DataInputStream inFile, int empty){
		try{inFile.mark(1609);}catch (Exception err){System.err.println("Error Mark: "+err.getMessage());}
		try{inFile.skip(12);}catch (Exception err){System.err.println("Error: "+err.getMessage());} //Offsets and value representations obtained from Stratec
		try{VoxelSize = Double.longBitsToDouble(Long.reverseBytes(inFile.readLong()));}catch (Exception err){System.err.println("Error: "+err.getMessage());}
		try{inFile.reset();}catch (Exception err){System.err.println("Error Reset1: "+err.getMessage());}
		try{inFile.skip(318);}catch (Exception err){System.err.println("Error: "+err.getMessage());} //Offsets and value representations obtained from Stratec
		try{ObjLen = Double.longBitsToDouble(Long.reverseBytes(inFile.readLong()));}catch (Exception err){System.err.println("Error: "+err.getMessage());}		
		try{inFile.reset();}catch (Exception err){System.err.println("Error Reset2: "+err.getMessage());}
		try{inFile.skip(662);}catch (Exception err){System.err.println("Error: "+err.getMessage());} //Offsets and value representations obtained from Stratec
		
		byte[] stringLength = new byte[1];
		try{inFile.read(stringLength);}catch (Exception err){System.err.println("Error: "+err.getMessage());}		//Check string length
		try{inFile.read(MInfo,0,stringLength[0]);}catch (Exception err){System.err.println("Error: "+err.getMessage());}
		MeasInfo = new String(MInfo);
		//inFile.read((char*) &mInfo,324);
		try{inFile.reset();}catch (Exception err){System.err.println("Error Reset3: "+err.getMessage());}
		try{inFile.skip(986);}catch (Exception err){System.err.println("Error: "+err.getMessage());} //Offsets and value representations obtained from Stratec
		try{MeasDate = Integer.reverseBytes(inFile.readInt());}catch (Exception err){System.err.println("Error: "+err.getMessage());}
		//MeasDate = (short) bytesFromFile[0];
		try{inFile.reset();}catch (Exception err){System.err.println("Error: "+err.getMessage());}
		try{inFile.skip(1050);}catch (Exception err){System.err.println("Error: "+err.getMessage());}
		try{inFile.read(stringLength);}catch (Exception err){System.err.println("Error: "+err.getMessage());}		//Tarkistetaan deviceStringin stringLength...
		try{inFile.read(Dev,0,stringLength[0]);}catch (Exception err){System.err.println("Error: "+err.getMessage());}
		Device = new String(Dev);		//DeviceType

		try{inFile.reset();}catch (Exception err){System.err.println("Error: "+err.getMessage());}
		try{inFile.skip(1085);}catch (Exception err){System.err.println("Error: "+err.getMessage());} //Offsets and value representations obtained from Stratec
		try{inFile.read(bytesFromFile,0,2);}catch (Exception err){System.err.println("Error: "+err.getMessage());}
		bits[0] = 0;
		bits[1] = 0;
		bits[0] = (0xFF & ((int)bytesFromFile[1]));
		bits[1] = (0xFF & ((int)bytesFromFile[0]));
		PatMeasNo = (int) (bits[0] << 8 | bits[1]);	//BigEndian encoding...
		try{inFile.reset();}catch (Exception err){System.err.println("Error: "+err.getMessage());}
		try{inFile.skip(1087);}catch (Exception err){System.err.println("Error: "+err.getMessage());} //Offsets and value representations obtained from Stratec
		try{PatNo = Integer.reverseBytes(inFile.readInt());}catch (Exception err){System.err.println("Error: "+err.getMessage());}
		
		try{inFile.reset();}catch (Exception err){System.err.println("Error: "+err.getMessage());}
		try{inFile.skip(1091);}catch (Exception err){System.err.println("Error: "+err.getMessage());} //Offsets and value representations obtained from Stratec
		try{PatBirth = Integer.reverseBytes(inFile.readInt());}catch (Exception err){System.err.println("Error: "+err.getMessage());}
		try{inFile.reset();}catch (Exception err){System.err.println("Error: "+err.getMessage());}
		try{inFile.skip(1099);}catch (Exception err){System.err.println("Error: "+err.getMessage());} //Offsets and value representations obtained from Stratec
		try{inFile.read(stringLength);}catch (Exception err){System.err.println("Error: "+err.getMessage());}		//Check string length
		try{inFile.read(PName,0,stringLength[0]);}catch (Exception err){System.err.println("Error: "+err.getMessage());}
		PatName = new String(PName);
		try{inFile.reset();}catch (Exception err){System.err.println("Error: "+err.getMessage());}
		try{inFile.skip(1282);}catch (Exception err){System.err.println("Error: "+err.getMessage());} //Offsets and value representations obtained from Stratec
		try{inFile.read(stringLength);}catch (Exception err){System.err.println("Error: "+err.getMessage());}		//Check string length
		try{inFile.read(PID,0,stringLength[0]);}catch (Exception err){System.err.println("Error: "+err.getMessage());}
		PatID = new String(PID);
		
		try{inFile.close();}catch (Exception err){System.err.println("Error: "+err.getMessage());}
	}
}
