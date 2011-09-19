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


package ImageReading;
import java.util.*;	//Vector
import java.io.*;	//File IO
import javax.swing.*;
import javax.imageio.*;
/*Java advanced imaging imageIO deactivated for now*/
/*
import com.sun.media.imageioimpl.plugins.jpeg.*;
import com.sun.medialib.codec.jpeg.*;
*/
//import com.sun.media.jai.codec.*;
//import javax.media.jai.*;

import java.awt.image.*;
import java.awt.*;
import java.util.*;	//Vector, Collections


public class DicomReader{
	
	public byte[] transferSyntax;
	public byte[] acquisitionDate;		//0008 0022
	public byte[] studyDescription;		//0008 1030
	public byte[] seriesDescription;		//0008 103E
	public byte[] patientName;		//0010 0010 
	public byte[] patientID;			//0010 0020
	public byte[] patientBirthDate;		//0010 0030
	public byte[] imageComments;		//0020 4000
	public int samplesPerPixel;				//0028 0002
	public int rows;				//0028 0010
	public int columns;				//0028 0011
	public byte[] pixelSpacing;		//0028 0030
	public int bitsAllocated;		//0028 0100
	public int bitsStored;			//0028 0101
	public int highBit;				//0028 0102
	public int pixelRepresentation;	//0028 0103
	public int lowestImagePixel;
	public int highestImagePixel;
	public byte[] windowCenter;
	public byte[] windowWidth;
	public byte[] rescaleIntercept;		//0028 1052
	public byte[] rescaleSlope;		//0028 1053
	public byte[] rescaleType;		//0028 1054
	public byte[] dataRepresentation;		//0028 9108
	public byte[] calibrationImage;		//0050 0004
	public byte[] imageID;		//0050 0400
	public byte[] units;		//0050 1001
	//public int numberOfSlices;				//0054 0081
	public short[] data;	//7FE0 0010
	public int[] data2;	//7FE0 0010
	long position;
	public short[] eTagRead;
	public double[] scaled;
	public double max;
	public double min;
	short[][] elementTag;
	int vars;
	int muuttujia;
	DataInputStream inFile;
	long fileLength;
	byte[] fileData;
	public boolean papyrus;
	boolean readDataToMemory;
	public Vector<DicomFile> files;
	//Functions..
	

	//Constructor
	
	public DicomReader(DataInputStream fileIn){
		inFile = fileIn;
		readDataToMemory = false;
		eTags();
		}
	
	public DicomReader(File fIn){
		fileLength = fIn.length();
		readDataToMemory = true;
		try{BufferedInputStream tied = new BufferedInputStream(new FileInputStream(fIn));
			inFile = new DataInputStream(tied);
			if (!inFile.markSupported()) {
				throw new RuntimeException("Mark/Reset not supported!");
			}
			eTags();
		 }catch (Exception err){System.err.println("Couldn't read Memoryfile: "+err.getMessage());}
	}
	
	public void eTags(){
		eTagRead = new short[4];	//Automatically initialized to zero
		position=0;
			//Read Dicom Header
		vars =0;
		muuttujia = 50;
		
		elementTag = new short[muuttujia][5];
		/*
		for (int j = 0;j<muuttujia;j++){
			elementTag[j] = new byte[5];
		}
		*/
		/*Transfer Syntax -always explicit VR!*/
		elementTag[vars][0]=(short)0x02;
		elementTag[vars][1]=(short)0x00;
		elementTag[vars][2]=(short)0x10;
		elementTag[vars][3]=(short)0x00;
		++vars;
		
		elementTag[vars][0]=(short)0x08;
		elementTag[vars][1]=(short)0x00;
		elementTag[vars][2]=(short)0x22;
		elementTag[vars][3]=(short)0x00;
		////System.out.println(elementTag[vars][0]+" "+elementTag[vars][1]+" "+elementTag[vars][2]+" "+elementTag[vars][3]);
		++vars;
		elementTag[vars][0]=(short)0x08;
		elementTag[vars][1]=(short)0x00;
		elementTag[vars][2]=(short)0x30;
		elementTag[vars][3]=(short)0x10;
		++vars;
		elementTag[vars][0]=(short)0x08;
		elementTag[vars][1]=(short)0x00;
		elementTag[vars][2]=(short)0x3E;
		elementTag[vars][3]=(short)0x10;
		++vars;
		elementTag[vars][0]=(short)0x10;
		elementTag[vars][1]=(short)0x00;
		elementTag[vars][2]=(short)0x10;
		elementTag[vars][3]=(short)0x00;
		++vars;
		elementTag[vars][0]=(short)0x10;
		elementTag[vars][1]=(short)0x00;
		elementTag[vars][2]=(short)0x20;
		elementTag[vars][3]=(short)0x00;
		++vars;
		elementTag[vars][0]=(short)0x10;
		elementTag[vars][1]=(short)0x00;
		elementTag[vars][2]=(short)0x30;
		elementTag[vars][3]=(short)0x00;
		++vars;
		elementTag[vars][0]=(short)0x20;
		elementTag[vars][1]=(short)0x00;
		elementTag[vars][2]=(short)0x00;
		elementTag[vars][3]=(short)0x40;
		++vars;
		elementTag[vars][0]=(short)0x28;
		elementTag[vars][1]=(short)0x00;
		elementTag[vars][2]=(short)0x02;
		elementTag[vars][3]=(short)0x00;
		++vars;
		elementTag[vars][0]=(short)0x28;
		elementTag[vars][1]=(short)0x00;
		elementTag[vars][2]=(short)0x10;
		elementTag[vars][3]=(short)0x00;
		++vars;
		elementTag[vars][0]=(short)0x28;
		elementTag[vars][1]=(short)0x00;
		elementTag[vars][2]=(short)0x11;
		elementTag[vars][3]=(short)0x00;
		++vars;
		elementTag[vars][0]=(short)0x28;
		elementTag[vars][1]=(short)0x00;
		elementTag[vars][2]=(short)0x30;
		elementTag[vars][3]=(short)0x00;
		++vars;
		elementTag[vars][0]=(short)0x28;
		elementTag[vars][1]=(short)0x00;
		elementTag[vars][2]=(short)0x00;
		elementTag[vars][3]=(short)0x01;
		++vars;
		elementTag[vars][0]=(short)0x28;
		elementTag[vars][1]=(short)0x00;
		elementTag[vars][2]=(short)0x01;
		elementTag[vars][3]=(short)0x01;
		++vars;
		elementTag[vars][0]=(short)0x28;
		elementTag[vars][1]=(short)0x00;
		elementTag[vars][2]=(short)0x02;
		elementTag[vars][3]=(short)0x01;
		++vars;
		elementTag[vars][0]=(short)0x28;
		elementTag[vars][1]=(short)0x00;
		elementTag[vars][2]=(short)0x03;
		elementTag[vars][3]=(short)0x01;
		++vars;
		elementTag[vars][0]=(short)0x28;
		elementTag[vars][1]=(short)0x00;
		elementTag[vars][2]=(short)0x52;
		elementTag[vars][3]=(short)0x10;
		++vars;
		elementTag[vars][0]=(short)0x28;
		elementTag[vars][1]=(short)0x00;
		elementTag[vars][2]=(short)0x53;
		elementTag[vars][3]=(short)0x10;
		++vars;
		elementTag[vars][0]=(short)0xE0;
		elementTag[vars][1]=(short)0x7F;
		elementTag[vars][2]=(short)0x10;
		elementTag[vars][3]=(short)0x00;
		++vars;
		
		long headerLength=0;
		try{inFile.mark(400);}catch (Exception err){System.err.println("Error Mark: "+err.getMessage());}	//Mark the file in case we need to rewind...
		try{inFile.skip(128+4+4+4);}catch (Exception err){System.err.println("Error: "+err.getMessage());} //DICM+0200+UL04	fseek(inFile,128+4+4+4, SEEK_SET);
		position+=140;
		try{headerLength = Integer.reverseBytes(inFile.readInt());}catch (Exception err){System.err.println("Error: "+err.getMessage());}
		position+=4;
		vars = 0;
		//System.out.println(elementTag[vars][0]+" "+elementTag[vars][1]+" "+elementTag[vars][2]+" "+elementTag[vars][3]);
		long temp_pos=position;
		transferSyntax = searchHeaderTag(inFile,elementTag,vars,1);
		try{inFile.skip((int) (headerLength-(position-temp_pos)));}catch (Exception err){System.err.println("Error: "+err.getMessage());} //DICM+0200+UL04	fseek(inFile,128+4+4+4, SEEK_SET);
		position = position +headerLength-(position-temp_pos);
		String implicit = new String("1.2.840.10008.1.2");
		String explicit = new String("1.2.840.10008.1.2.1");
		String losslessJPEG = new String("1.2.840.10008.1.2.4.70");
		String tSyntax = new String(transferSyntax).trim();
		
		/*Read the whole file to memory at once*/
		if (readDataToMemory){
			try{inFile.reset();}catch (Exception err){System.err.println("Error Mark: "+err.getMessage());}	//Mark the file in case we need to rewind...
			fileData = new byte[(int) fileLength];		//Allocate memory for reading the file into memory
			int dataRead = 0;
			try{dataRead = inFile.read(fileData,0,(int) fileLength);}catch (Exception err){System.err.println("Error: "+err.getMessage());} //Read the data to memory
		}
		
		papyrus = false;
		if (tSyntax.equalsIgnoreCase(implicit)==true){
			System.out.println("Implicit VR "+position);
			if (readDataToMemory){
				readImplicitMemory();
			}else{
				readImplicit();
			}
		}
		if (tSyntax.equalsIgnoreCase(explicit)==true){
			System.out.println("Explicit VR");
			if (readDataToMemory){
				readExplicitMemory();
			}else{
				readExplicit();
			}
			
		}
		if (tSyntax.equalsIgnoreCase(losslessJPEG)==true){
			System.out.println("Lossless JPEG (made for PAPYRUS file, which uses Explicit VR)");
			papyrus = true;
			readJPEG();
		}
		try{inFile.close();}catch (Exception err){System.out.println("Couldn't close inFile: "+err);}
	}
	
	void readJPEG(){ //Reads Papyrus files, which use the lossless JPEG explicit transfer syntax..
		readPapyrus();	//Reads papyrusfile Structure
	} 
	
	void readPapyrus(){
		files = new Vector<DicomFile>();

		vars =0;
		/*
		muuttujia = 19;
		
		elementTag = new short[muuttujia][5];
		*/

		//Transfer Syntax -always explicit VR!
		vars =0;
		/*Papyrus File item offsets...*/
		elementTag[vars][0]=(short)0x41;
		elementTag[vars][1]=(short)0x00;
		elementTag[vars][2]=(short)0x10;
		elementTag[vars][3]=(short)0x10;
		++vars;
		/*Papyrus No of Items = how many slices*/
		elementTag[vars][0]=(short)0x41;
		elementTag[vars][1]=(short)0x00;
		elementTag[vars][2]=(short)0x15;
		elementTag[vars][3]=(short)0x10;
		
		vars = 0;
		byte[] papyrusFilelist = searchDicomTagMemory(fileData,elementTag,vars,1);
//		System.out.println(papyrusFilelist.length);
		++vars;
		int papyrusSliceNo = searchDicomTagMemory(fileData,elementTag,vars);
	//	System.out.println(papyrusSliceNo);
		/*Get file offsets from papyrusFilelist*/
		vars =0;
		/*Papyrus File item offsets...*/
		elementTag[vars][0]=(short)0x41;
		elementTag[vars][1]=(short)0x00;
		elementTag[vars][2]=(short)0x11;
		elementTag[vars][3]=(short)0x10;
		++vars;
		/*Item delimitation tag*/
		elementTag[vars][0]=(short)0xFE;
		elementTag[vars][1]=(short)0xFF;
		elementTag[vars][2]=(short)0x00;
		elementTag[vars][3]=(short)0xE0;

		long[] offsets = new long[papyrusSliceNo];
		position = 8;
		/*Find the offsets to the data slices*/
		for (int i = 0; i<papyrusSliceNo; i++){
			vars = 0;
			//System.out.println(position);
			offsets[i] = searchDicomTagMemory2(papyrusFilelist,elementTag,vars);
			
			//System.out.println(offsets[i]);
			++vars;
			if (i <papyrusSliceNo-1){
				long useless = searchDicomTagMemory2(papyrusFilelist,elementTag,vars); //Read until item delimitation
			}
		
		}
		/*Set elementTags to regular DICOM file type...*/
		setElementTagsPapyrus();
		
		/*Go through the items*/
		//System.out.println("Go through");
		for (int i = 0; i<papyrusSliceNo; i++){
			files.add(new DicomFile());
		//for (int i = 0; i<1; i++){
			vars = 0;
			position = offsets[i];
			/*Go through dicom header here*/
			//System.out.println("DicomFilea tayttamaan");
			files.get(i).acquisitionDate = searchDicomTagMemory(fileData,elementTag,vars,1);
			//System.out.println("aDate "+new String(acquisitionDate) +" pos "+position);
			
			++vars;
			files.get(i).studyDescription = searchDicomTagMemory(fileData,elementTag,vars,1);
			//System.out.println("studyDesc "+new String(studyDescription) +" pos "+position);
			
			++vars;
			files.get(i).seriesDescription = searchDicomTagMemory(fileData,elementTag,vars,1);
			//System.out.println("seriesDesc "+new String(seriesDescription) +" pos "+position);
			
			++vars;
			files.get(i).patientName = searchDicomTagMemory(fileData,elementTag,vars,1);
			//System.out.println("pName "+new String(patientName) +" pos "+position);
			
			++vars;
			files.get(i).patientID = searchDicomTagMemory(fileData,elementTag,vars,1);
			//System.out.println("pId "+new String(patientID) +" pos "+position);
			
			++vars;
			files.get(i).patientBirthDate = searchDicomTagMemory(fileData,elementTag,vars,1);
			//System.out.println("pBirth "+new String(patientBirthDate) +" pos "+position);
			++vars;
			files.get(i).imageComments = searchDicomTagMemory(fileData,elementTag,vars,1);
			
			++vars;	
			files.get(i).samplesPerPixel = searchDicomTagMemory(fileData,elementTag,vars);
			++vars;
			files.get(i).rows = searchDicomTagMemory(fileData,elementTag,vars);
			++vars;
			files.get(i).columns = searchDicomTagMemory(fileData,elementTag,vars);
			++vars;
			files.get(i).pixelSpacing = searchDicomTagMemory(fileData,elementTag,vars,1);
			++vars;
			files.get(i).bitsAllocated = searchDicomTagMemory(fileData,elementTag,vars);
			++vars;
			files.get(i).bitsStored = searchDicomTagMemory(fileData,elementTag,vars);
			++vars;
			files.get(i).highBit = searchDicomTagMemory(fileData,elementTag,vars);
			++vars;
			files.get(i).pixelRepresentation = searchDicomTagMemory(fileData,elementTag,vars);
			++vars;
			files.get(i).lowestImagePixel = searchDicomTagMemory(fileData,elementTag,vars);
			++vars;
			files.get(i).highestImagePixel = searchDicomTagMemory(fileData,elementTag,vars);
			++vars;
			files.get(i).windowCenter = searchDicomTagMemory(fileData,elementTag,vars,1);
			++vars;
			files.get(i).windowWidth = searchDicomTagMemory(fileData,elementTag,vars,1);
			++vars;
			files.get(i).rescaleIntercept = searchDicomTagMemory(fileData,elementTag,vars,1);
			++vars;
			files.get(i).rescaleSlope = searchDicomTagMemory(fileData,elementTag,vars,1);
			++vars;
			//System.out.println("Lo "+ files.get(i).lowestImagePixel+" Hi "+ files.get(i).highestImagePixel+" Center "+new String(files.get(i).windowCenter)+" Width "+new String(files.get(i).windowWidth));
				/*Header Done*/
				/*Read lossless JPEG*/
			//System.out.println("Item delimitation item");
			byte[] useless2 = searchDicomTagMemory(fileData,elementTag,vars,1);
			//System.out.println("Read JPEG");
			byte[] compressedJPEG = readCompressedDataMemory(fileData,elementTag,vars);
			//System.out.println("Decompress JPEG "+compressedJPEG.length);
			ByteArrayInputStream in = new ByteArrayInputStream(compressedJPEG);
			try{BufferedImage img = javax.imageio.ImageIO.read(in);
				//System.out.println(img.getWidth()+" "+img.getHeight()+" "+img.getType());
				//System.out.println(BufferedImage.TYPE_USHORT_GRAY);
				files.get(i).columns =img.getWidth();
				files.get(i).rows = img.getHeight();
				Raster raster = img.getData();
				DataBuffer dataBuffer = raster.getDataBuffer();
				//System.out.println("Datatyyppi "+dataBuffer.getDataType()+" UShort "+DataBuffer.TYPE_USHORT);
				//System.out.println("Banks "+dataBuffer.getNumBanks()+" size "+dataBuffer.getSize());
				files.get(i).data2 = new int[dataBuffer.getSize()];
				for (int j=0; j<dataBuffer.getSize();++j){
					files.get(i).data2[j]= (dataBuffer.getElem(j))>>1;
				}
				files.get(i).min = (int) files.get(i).lowestImagePixel;
				files.get(i).max = (int) files.get(i).highestImagePixel;
				for (int j=0; j<dataBuffer.getSize();++j){
					if (files.get(i).data2[j] < files.get(i).min | files.get(i).data2[j] > files.get(i).max){
						files.get(i).data2[j]= 0;
					}
				}
				/*
					int[] temp = (int[]) files.get(i).data2.clone();
					Arrays.sort(temp);
					files.get(i).max = temp[temp.length-1-img.getWidth()*2];
					files.get(i).min = temp[img.getWidth()*2];
				*/
			} catch (Exception err){System.err.println("Lossless JPEG: "+err.getMessage());}
		}
		
		////System.out.println(elementTag[vars][0]+" "+elementTag[vars][1]+" "+elementTag[vars][2]+" "+elementTag[vars][3]);
	}
	
	void setElementTagsPapyrus(){
		vars = 0;
		elementTag[vars][0]=(short)0x08;
		elementTag[vars][1]=(short)0x00;
		elementTag[vars][2]=(short)0x20;
		elementTag[vars][3]=(short)0x00;
		++vars;
		elementTag[vars][0]=(short)0x08;
		elementTag[vars][1]=(short)0x00;
		elementTag[vars][2]=(short)0x30;
		elementTag[vars][3]=(short)0x10;
		++vars;
		elementTag[vars][0]=(short)0x08;
		elementTag[vars][1]=(short)0x00;
		elementTag[vars][2]=(short)0x3E;
		elementTag[vars][3]=(short)0x10;
		++vars;
		elementTag[vars][0]=(short)0x10;
		elementTag[vars][1]=(short)0x00;
		elementTag[vars][2]=(short)0x10;
		elementTag[vars][3]=(short)0x00;
		++vars;
		elementTag[vars][0]=(short)0x10;
		elementTag[vars][1]=(short)0x00;
		elementTag[vars][2]=(short)0x20;
		elementTag[vars][3]=(short)0x00;
		++vars;
		elementTag[vars][0]=(short)0x10;
		elementTag[vars][1]=(short)0x00;
		elementTag[vars][2]=(short)0x30;
		elementTag[vars][3]=(short)0x00;
		++vars;
		elementTag[vars][0]=(short)0x20;
		elementTag[vars][1]=(short)0x00;
		elementTag[vars][2]=(short)0x00;
		elementTag[vars][3]=(short)0x40;
		++vars;
		elementTag[vars][0]=(short)0x28;
		elementTag[vars][1]=(short)0x00;
		elementTag[vars][2]=(short)0x02;
		elementTag[vars][3]=(short)0x00;
		++vars;
		elementTag[vars][0]=(short)0x28;
		elementTag[vars][1]=(short)0x00;
		elementTag[vars][2]=(short)0x10;
		elementTag[vars][3]=(short)0x00;
		++vars;
		elementTag[vars][0]=(short)0x28;
		elementTag[vars][1]=(short)0x00;
		elementTag[vars][2]=(short)0x11;
		elementTag[vars][3]=(short)0x00;
		++vars;
		elementTag[vars][0]=(short)0x28;
		elementTag[vars][1]=(short)0x00;
		elementTag[vars][2]=(short)0x30;
		elementTag[vars][3]=(short)0x00;
		++vars;
		elementTag[vars][0]=(short)0x28;
		elementTag[vars][1]=(short)0x00;
		elementTag[vars][2]=(short)0x00;
		elementTag[vars][3]=(short)0x01;
		++vars;
		elementTag[vars][0]=(short)0x28;
		elementTag[vars][1]=(short)0x00;
		elementTag[vars][2]=(short)0x01;
		elementTag[vars][3]=(short)0x01;
		++vars;
		elementTag[vars][0]=(short)0x28;
		elementTag[vars][1]=(short)0x00;
		elementTag[vars][2]=(short)0x02;
		elementTag[vars][3]=(short)0x01;
		++vars;
		elementTag[vars][0]=(short)0x28;
		elementTag[vars][1]=(short)0x00;
		elementTag[vars][2]=(short)0x03;
		elementTag[vars][3]=(short)0x01;
		++vars;
		/*Lowest image pixel value US/SS*/
		elementTag[vars][0]=(short)0x28;
		elementTag[vars][1]=(short)0x00;
		elementTag[vars][2]=(short)0x06;
		elementTag[vars][3]=(short)0x01;
		++vars;
		/*Highest image pixel value US/SS*/
		elementTag[vars][0]=(short)0x28;
		elementTag[vars][1]=(short)0x00;
		elementTag[vars][2]=(short)0x07;
		elementTag[vars][3]=(short)0x01;
		++vars;
		/*Window center DS*/
		elementTag[vars][0]=(short)0x28;
		elementTag[vars][1]=(short)0x00;
		elementTag[vars][2]=(short)0x50;
		elementTag[vars][3]=(short)0x10;
		++vars;
		/*Window Width DS*/
		elementTag[vars][0]=(short)0x28;
		elementTag[vars][1]=(short)0x00;
		elementTag[vars][2]=(short)0x51;
		elementTag[vars][3]=(short)0x10;
		
		++vars;
		elementTag[vars][0]=(short)0x28;
		elementTag[vars][1]=(short)0x00;
		elementTag[vars][2]=(short)0x52;
		elementTag[vars][3]=(short)0x10;
		++vars;
		elementTag[vars][0]=(short)0x28;
		elementTag[vars][1]=(short)0x00;
		elementTag[vars][2]=(short)0x53;
		elementTag[vars][3]=(short)0x10;
		++vars;
		/*Set to Item delimitation item*/
		elementTag[vars][0]=(short)0xFE;
		elementTag[vars][1]=(short)0xFF;
		elementTag[vars][2]=(short)0x00;
		elementTag[vars][3]=(short)0xE0;
	}
	
	void memorySearchTag(byte[] fileData,short[][] elementTag,int vars){
		byte[] swapBytes = new byte[4];
		byte VR[] = new byte[2];
		String stringVR;
		int elementLength=0;
		long	elementLengthVariant=0;
			for (int i = 0; i<4;i++){swapBytes[i]=fileData[(int)position+i];}	//Element Name
		position+=4;
		for (int i =0;i<4;i++){eTagRead[i] = (short) (0x000000FF & (int)swapBytes[i]);}
		while ((eTagRead[0] != elementTag[vars][0] || eTagRead[1] != elementTag[vars][1] || eTagRead[2] != elementTag[vars][2]  || eTagRead[3] != elementTag[vars][3]) //As long as the tag does not match
		&& (((long) (eTagRead[1])<<24 | (long) (eTagRead[0])<<16 | (long) (eTagRead[3])<<8 | (long) (eTagRead[2])) < ((long) (elementTag[vars][1])<<24 | (long) (elementTag[vars][0])<<16 | (long) (elementTag[vars][3])<<8 | (long) (elementTag[vars][2])))
		 //&&//As long as the tag in file is smaller than the tag we are looking for
		//ftell(inFile) < fileLength
		){ 
			for (int i =0;i<2;i++){VR[i] = fileData[(int)position+i];} //Value representation
			position+=2;
			//Check whether VR is one of the special ones
			stringVR = new String(VR);
			if (stringVR.equals(new String("OB")) == false &&	stringVR.equals(new String("OW")) == false &&	stringVR.equals(new String("OF")) == false &&	stringVR.equals(new String("SQ")) == false &&	stringVR.equals(new String("UT")) == false &&	stringVR.equals(new String("UN")) == false){
				for (int i = 0; i<2;i++){swapBytes[i]=fileData[(int)position+i];}	//Element length
				position+=2;
				elementLength = (int) swapBytes[1]<<8|swapBytes[0];
				if (elementLength > 0){
					position+=elementLength;	//Skip element
				}
			} else {
				position+=2;	//Skip empty
				for (int i = 0; i<4;i++){swapBytes[i]=fileData[(int)position+i];}	//Element length variant
				position+=4;
				elementLengthVariant = (long) ((0x000000ff &(int) swapBytes[3])<<24 | (0x000000ff &(int) swapBytes[2])<<16 | (0x000000ff &(int) swapBytes[1])<<8 | (0x000000ff &(int) swapBytes[0]));
				if (elementLengthVariant > 0){
					position+=elementLengthVariant; //Skip element
				}
			}
			

			for (int i = 0; i<4;i++){swapBytes[i]=fileData[(int)position+i];}	//Element Nam
			position+=4;
			for (int i =0;i<4;i++){eTagRead[i] = (short) (0x000000FF & (int)swapBytes[i]);}
		}
	}
	
	long searchDicomTagMemory2(byte[] fileData,short[][] elementTag,int vars){
		long returnValue=0;
		byte[] buffer = new byte[500];
		byte[] swapBytes = new byte[4];
		short[] swapShort = new short[4];
		
		byte VR[] = new byte[2];
		String stringVR;
		int elementLength=0;
		long	elementLengthVariant=0;
		//Read Data Section  0008
		//Data elements are in increasing order -> read the file until finding what we're looking for
		//Look for Acquisition date...
		memorySearchTag(fileData,elementTag,vars);
				
		if(eTagRead[0] == elementTag[vars][0] && eTagRead[1] == elementTag[vars][1] && eTagRead[2] == elementTag[vars][2]  && eTagRead[3] == elementTag[vars][3]){
			for (int i =0;i<2;i++){VR[i] = fileData[(int)position+i];} //Value representation
			position+=2;
			//Check whether VR is one of the special ones
			stringVR = new String(VR);
			if (stringVR.equals(new String("OB")) == false &&	stringVR.equals(new String("OW")) == false &&	stringVR.equals(new String("OF")) == false &&	stringVR.equals(new String("SQ")) == false &&	stringVR.equals(new String("UT")) == false &&	stringVR.equals(new String("UN")) == false){
					for (int i = 0; i<2;i++){swapBytes[i]=fileData[(int)position+i];}	//Element length
					position+=2;
					elementLength = (int) swapBytes[1]<<8|swapBytes[0];
					if (elementLength > 0){
						//Find out whether it is unsigned or signed short or integer....
						for (int i = 0; i<elementLength;i++){swapBytes[i]=fileData[(int)position+i];}	//Read Element
						position+=elementLength;
						for (int i =0;i<4;i++){swapShort[i] = (short) (0x000000FF & (int)swapBytes[i]);}
						returnValue = (long) ((long) (swapShort[3])<<24|(long) (swapShort[2])<<16 | (long) (swapShort[1])<<8|(long) (swapShort[0]));
					}
				} else {
					position+=2;	//Skip empty element length
					for (int i = 0; i<4;i++){swapBytes[i]=fileData[(int)position+i];}	//Element length variant
					position+=4;
					elementLengthVariant = (long) ((0x000000ff &(int) swapBytes[3])<<24 | (0x000000ff &(int) swapBytes[2])<<16 | (0x000000ff &(int) swapBytes[1])<<8 | (0x000000ff &(int) swapBytes[0]));
					for (int i = 0; i<(int) elementLengthVariant;i++){swapBytes[i]=fileData[(int)position+i];}	//Read Element
					position+=elementLengthVariant;
					for (int i =0;i<4;i++){swapShort[i] = (short) (0x000000FF & (int)swapBytes[i]);}
						returnValue = (long) ((long) (swapShort[3])<<24|(long) (swapShort[2])<<16 | (long) (swapShort[1])<<8|(long) (swapShort[0]));
				}
		}else{
			position-=4;
		}
		return returnValue;
	}
	
	
	int searchDicomTagMemory(byte[] fileData,short[][] elementTag,int vars){
		int returnValue=0;
		byte[] buffer = new byte[500];
		byte[] swapBytes = new byte[4];
		short[] swapShort = new short[4];
		byte VR[] = new byte[2];
		String stringVR;
		int elementLength=0;
		long	elementLengthVariant=0;
		//Read Data Section  0008
		//Data elements are in increasing order -> read the file until finding what we're looking for
		//Look for Acquisition date...
		memorySearchTag(fileData,elementTag,vars);
		
		if(eTagRead[0] == elementTag[vars][0] && eTagRead[1] == elementTag[vars][1] && eTagRead[2] == elementTag[vars][2]  && eTagRead[3] == elementTag[vars][3]){
			for (int i =0;i<2;i++){VR[i] = fileData[(int)position+i];} //Value representation
			position+=2;
			//Check whether VR is one of the special ones
			stringVR = new String(VR);
			if (stringVR.equals(new String("OB")) == false &&	stringVR.equals(new String("OW")) == false &&	stringVR.equals(new String("OF")) == false &&	stringVR.equals(new String("SQ")) == false &&	stringVR.equals(new String("UT")) == false &&	stringVR.equals(new String("UN")) == false){
					for (int i = 0; i<2;i++){swapBytes[i]=fileData[(int)position+i];}	//Element length
					position+=2;
					elementLength = (int) swapBytes[1]<<8|swapBytes[0];
					if (elementLength > 0){
						//Find out whether it is unsigned or signed short or integer....
						for (int i = 0; i<elementLength;i++){swapBytes[i]=fileData[(int)position+i];}	//Read Element
						position+=elementLength;
						returnValue = (int) swapBytes[1]<<8|swapBytes[0];
						if (stringVR.equals(new String("US")) ==true){ //If returning US..
							for (int i =0;i<4;i++){swapShort[i] = (short) (0x000000FF & (int)swapBytes[i]);}
							int temp=(int) ((int) (swapShort[1])<<8|(int) (swapShort[0]));
							returnValue = (short) temp;
						}
					}
				} else {
					position+=2;	//Skip empty element length
					for (int i = 0; i<4;i++){swapBytes[i]=fileData[(int)position+i];}	//Element length variant
					position+=4;
					elementLengthVariant = (long) ((0x000000ff &(int) swapBytes[3])<<24 | (0x000000ff &(int) swapBytes[2])<<16 | (0x000000ff &(int) swapBytes[1])<<8 | (0x000000ff &(int) swapBytes[0]));
					if (elementLengthVariant > 0){
						for (int i = 0; i<(int) elementLengthVariant;i++){swapBytes[i]=fileData[(int)position+i];}	//Read Element
						position+=elementLengthVariant;
						returnValue = (int) ((0x000000ff &(int) swapBytes[1])<<8 | (0x000000ff &(int) swapBytes[0]));
					}
				}
		}else{
			position-=4;
		}
		return returnValue;
	}
	
	byte[] searchDicomTagMemory(byte[] fileData,short[][] elementTag, int vars, int reserveMemory){
		byte[] returnValue =new byte[0];
		byte[] buffer = new byte[500];
		byte[] swapBytes = new byte[4];
		byte VR[] = new byte[2];
		String stringVR;
		int elementLength=0;
		long	elementLengthVariant=0;
		//Read Data Section  0008
		//Data elements are in increasing order -> read the file until finding what we're looking for
		//Look for Acquisition date...
		//System.out.println("Starting search "+ elementTag[vars][0]+" "+elementTag[vars][1]+" "+elementTag[vars][2]+" "+elementTag[vars][3]+" pos "+position);
		memorySearchTag(fileData,elementTag,vars);
		//System.out.println("Done MemorySearching "+ eTagRead[0]+" "+eTagRead[1]+" "+eTagRead[2]+" "+eTagRead[3]+ " eT "+ elementTag[vars][0]+" "+elementTag[vars][1]+" "+elementTag[vars][2]+" "+elementTag[vars][3]);
		if(eTagRead[0] == elementTag[vars][0] && eTagRead[1] == elementTag[vars][1] && eTagRead[2] == elementTag[vars][2]  && eTagRead[3] == elementTag[vars][3]){
			for (int i =0;i<2;i++){VR[i] = fileData[(int)position+i];} //Value representation
			position+=2;
			//Check whether VR is one of the special ones
			stringVR = new String(VR);
			//System.out.println(stringVR);
			if (stringVR.equals(new String("OB")) == false &&	stringVR.equals(new String("OW")) == false &&	stringVR.equals(new String("OF")) == false &&	stringVR.equals(new String("SQ")) == false &&	stringVR.equals(new String("UT")) == false &&	stringVR.equals(new String("UN")) == false){
				for (int i = 0; i<2;i++){swapBytes[i]=fileData[(int)position+i];}	//Element length
				position+=2;
				elementLength = (int) swapBytes[1]<<8|swapBytes[0];
				if (elementLength > 0){
						if (reserveMemory==1){
							returnValue = new byte[elementLength+1];		//0008 0022
						}
					//Find out whether it is unsigned or signed short or integer....
					for (int i = 0; i<elementLength;i++){returnValue[i]=fileData[(int)position+i];}	//Read Element
					position+=elementLength;
					//System.out.println("ReserveMem "+new String(returnValue));
				}
			} else {
				position+=2;	//Skip empty element length
				for (int i = 0; i<4;i++){swapBytes[i]=fileData[(int)position+i];}	//Element length variant
				position+=4;
				elementLengthVariant = (long) ((0x000000ff &(int) swapBytes[3])<<24 | (0x000000ff &(int) swapBytes[2])<<16 | (0x000000ff &(int) swapBytes[1])<<8 | (0x000000ff &(int) swapBytes[0]));
				if (elementLengthVariant > 0){
					if (reserveMemory==1){
						returnValue = new byte[(int) elementLengthVariant+1];		//0008 0022
					}
					for (int i = 0; i<(int) elementLengthVariant;i++){returnValue[i]=fileData[(int)position+i];}	//Read Element
					position+=elementLengthVariant;
				}
			}
		}else{
			position-=4;	//Rewind the file to the beginning of the last Element Name
		}
		return returnValue;
	}
	
	byte[] readCompressedDataMemory(byte[] fileData,short[][] elementTag,int vars){

		byte[] buffer = new byte[500];
		byte[] swapBytes = new byte[4];
		byte VR[] = new byte[2];
		String stringVR;
		int elementLength=0;
		long	elementLengthVariant=0;
		byte[] readData=null;
		memorySearchTag(fileData,elementTag,vars); /*Search to the proper position*/
		//Read Data Section  0008
		//Data elements are in increasing order -> read the file until finding what we're looking for
		//Look for Acquisition date...
		
		if(eTagRead[0] == elementTag[vars][0] && eTagRead[1] == elementTag[vars][1] && eTagRead[2] == elementTag[vars][2]  && eTagRead[3] == elementTag[vars][3]){
			for (int i = 0; i<4;i++){swapBytes[i]=fileData[(int)position+i];}	//Element length variant
			position+=4;
			elementLengthVariant = (long) ((0x000000ff &(int) swapBytes[3])<<24 | (0x000000ff &(int) swapBytes[2])<<16 | (0x000000ff &(int) swapBytes[1])<<8 | (0x000000ff &(int) swapBytes[0]));
			readData = new byte[(int) elementLengthVariant];
			//System.out.println("Position for reading "+position);
			for (int i = 0; i<(int) elementLengthVariant;i++){readData[i]=fileData[(int)position+i];}	//Read Element
			position+=elementLengthVariant;
		}
		return readData;
	}

	
	/*Implicit value representation*/
	void readImplicitMemory(){
	
		//Data elements are in increasing order -> read the file until finding what we're looking for
		++vars;
		acquisitionDate = searchDicomTagMemoryImplicit(fileData,elementTag,vars,1);
		//System.out.println("aDate "+new String(acquisitionDate));
		
		++vars;
		studyDescription = searchDicomTagMemoryImplicit(fileData,elementTag,vars,1);
		++vars;
		seriesDescription = searchDicomTagMemoryImplicit(fileData,elementTag,vars,1);
		++vars;
		patientName = searchDicomTagMemoryImplicit(fileData,elementTag,vars,1);
		++vars;
		patientID = searchDicomTagMemoryImplicit(fileData,elementTag,vars,1);
		++vars;
		patientBirthDate = searchDicomTagMemoryImplicit(fileData,elementTag,vars,1);
		//System.out.println("pBirth "+new String(patientBirthDate));
		++vars;
		imageComments = searchDicomTagMemoryImplicit(fileData,elementTag,vars,1);
		
		++vars;	
		samplesPerPixel = searchDicomTagMemoryImplicit(fileData,elementTag,vars);
		++vars;
		rows = searchDicomTagMemoryImplicit(fileData,elementTag,vars);
		++vars;
		columns = searchDicomTagMemoryImplicit(fileData,elementTag,vars);
		++vars;
		pixelSpacing = searchDicomTagMemoryImplicit(fileData,elementTag,vars,1);
		++vars;
		bitsAllocated = searchDicomTagMemoryImplicit(fileData,elementTag,vars);
		++vars;
		bitsStored = searchDicomTagMemoryImplicit(fileData,elementTag,vars);
		++vars;
		highBit = searchDicomTagMemoryImplicit(fileData,elementTag,vars);
		++vars;
		pixelRepresentation = searchDicomTagMemoryImplicit(fileData,elementTag,vars);
		++vars;
		rescaleIntercept = searchDicomTagMemoryImplicit(fileData,elementTag,vars,1);
		++vars;
		rescaleSlope = searchDicomTagMemoryImplicit(fileData,elementTag,vars,1);
		++vars;
		//Header info read, read DATA
		//System.out.println("Next step is reading data");
		data = new short[rows*columns];
		
		//System.out.println("Dataa lukemaan "+rows*columns);
		readDataMemoryImplicit(fileData,elementTag,vars, data);
		short[] temp = (short[]) data.clone();
		Arrays.sort(temp);
		max = temp[temp.length-1];
		min = temp[0];
	}
	
	void readImplicit(){	
		//Data elements are in increasing order -> read the file until finding what we're looking for
		++vars;
		acquisitionDate = searchDicomTagImplicit(inFile,elementTag,vars,1);
		//System.out.println("aDate "+new String(acquisitionDate) +" pos "+position);
		
		++vars;
		studyDescription = searchDicomTagImplicit(inFile,elementTag,vars,1);
		//System.out.println("sDesc "+new String(studyDescription) +" pos "+position);
		
		++vars;
		seriesDescription = searchDicomTagImplicit(inFile,elementTag,vars,1);
		//System.out.println("sDesc "+new String(seriesDescription) +" pos "+position);
		
		++vars;
		patientName = searchDicomTagImplicit(inFile,elementTag,vars,1);
		//System.out.println("pName "+new String(patientName) +" pos "+position);
		
		++vars;
		patientID = searchDicomTagImplicit(inFile,elementTag,vars,1);
		//System.out.println("pId "+new String(patientID) +" pos "+position);
		
		++vars;
		patientBirthDate = searchDicomTagImplicit(inFile,elementTag,vars,1);
		//System.out.println("pBirth "+new String(patientBirthDate));
		++vars;
		imageComments = searchDicomTagImplicit(inFile,elementTag,vars,1);
		
		++vars;	
		samplesPerPixel = searchDicomTagImplicit(inFile,elementTag,vars);
		++vars;
		rows = searchDicomTagImplicit(inFile,elementTag,vars);
		++vars;
		columns = searchDicomTagImplicit(inFile,elementTag,vars);
		++vars;
		pixelSpacing = searchDicomTagImplicit(inFile,elementTag,vars,1);
		++vars;
		bitsAllocated = searchDicomTagImplicit(inFile,elementTag,vars);
		++vars;
		bitsStored = searchDicomTagImplicit(inFile,elementTag,vars);
		++vars;
		highBit = searchDicomTagImplicit(inFile,elementTag,vars);
		++vars;
		pixelRepresentation = searchDicomTagImplicit(inFile,elementTag,vars);
		++vars;
		rescaleIntercept = searchDicomTagImplicit(inFile,elementTag,vars,1);
		++vars;
		rescaleSlope = searchDicomTagImplicit(inFile,elementTag,vars,1);
		++vars;
		//Header info read, read DATA
		
		data = new short[rows*columns];
		
		//System.out.println("Dataa lukemaan "+rows*columns);
		readDataImplicit(inFile,elementTag,vars, data);
		short[] temp = (short[]) data.clone();
		Arrays.sort(temp);
		max = temp[temp.length-1];
		min = temp[0];
	}
	
	/*Implicit VR*/
	int searchDicomTagImplicit(DataInputStream inFile,short[][] elementTag,int vars){
		int returnValue=0;
		byte[] buffer = new byte[500];
		byte[] swapBytes = new byte[4];
		long elementLength=0;
		//Read Data Section  0008
		//Data elements are in increasing order -> read the file until finding what we're looking for
		//Look for Acquisition date...
		try{inFile.read(swapBytes,0,4);}catch (Exception err){System.err.println("Error: "+err.getMessage());}//Element Name
		position+=4;
		for (int i =0;i<4;i++){eTagRead[i] = (short) (0x000000FF & (int)swapBytes[i]);}
		while ((eTagRead[0] != elementTag[vars][0] || eTagRead[1] != elementTag[vars][1] || eTagRead[2] != elementTag[vars][2]  || eTagRead[3] != elementTag[vars][3]) //As long as the tag does not match
		&& (((long) (eTagRead[1])<<24 | (long) (eTagRead[0])<<16 | (long) (eTagRead[3])<<8 | (long) (eTagRead[2])) < ((long) (elementTag[vars][1])<<24 | (long) (elementTag[vars][0])<<16 | (long) (elementTag[vars][3])<<8 | (long) (elementTag[vars][2])))		 //&&//As long as the tag in file is smaller than the tag we are looking for
		){ 
			try{inFile.read(swapBytes,0,4);}catch (Exception err){System.err.println("Error: "+err.getMessage());} //Element length
			position+=4;
			elementLength = (long) ((0x000000ff &(int) swapBytes[3])<<24 | (0x000000ff &(int) swapBytes[2])<<16 | (0x000000ff &(int) swapBytes[1])<<8 | (0x000000ff &(int) swapBytes[0]));
			if (elementLength > 0){// && (int) (eTagRead[3]<<8 | eTagRead[2]) > 0 ){ //Don't skip over whole sections..
				if(elementLength > 500){buffer = new byte[(int) elementLength];}
				try{inFile.read(buffer,0,(int) elementLength);}catch (Exception err){System.err.println("Error: "+err.getMessage());}	//Read Element
				position+=elementLength;
			}
			try{inFile.mark(4);}catch (Exception err){System.err.println("Error Mark: "+err.getMessage());}	//Mark the file in case we need to rewind...
			try{inFile.read(swapBytes,0,4);}catch (Exception err){System.err.println("Error: "+err.getMessage());}//Element Name
			position+=4;
			for (int i =0;i<4;i++){eTagRead[i] = (short) (0x000000FF & (int)swapBytes[i]);}
		}
				
		if(eTagRead[0] == elementTag[vars][0] && eTagRead[1] == elementTag[vars][1] && eTagRead[2] == elementTag[vars][2]  && eTagRead[3] == elementTag[vars][3]){
			try{inFile.read(swapBytes,0,4);}catch (Exception err){System.err.println("Error: "+err.getMessage());} //Element length
			position+=4;
			elementLength = (long) ((0x000000ff &(int) swapBytes[3])<<24 | (0x000000ff &(int) swapBytes[2])<<16 | (0x000000ff &(int) swapBytes[1])<<8 | (0x000000ff &(int) swapBytes[0]));
			if (elementLength > 0){
				try{inFile.read(swapBytes,0,(int) elementLength);}catch (Exception err){System.err.println("Error: "+err.getMessage());}	//Read Element
				position+=elementLength;
				returnValue = (int) swapBytes[1]<<8|swapBytes[0];
			}
		}else{
			try{inFile.reset();}catch (Exception err){System.err.println("Error Reset1: "+err.getMessage());}	//Rewind the file to the beginning of the last Element Name
		}
		return returnValue;
	}
	
	
	byte[] searchDicomTagImplicit(DataInputStream inFile,short[][] elementTag, int vars, int reserveMemory){
		byte[] returnValue =new byte[0];
		byte[] buffer = new byte[500];
		byte[] swapBytes = new byte[4];
		long elementLength=0;

		//Read Data Section  0008
		//Data elements are in increasing order -> read the file until finding what we're looking for
		//Look for Acquisition date...
		try{inFile.read(swapBytes,0,4);}catch (Exception err){System.err.println("Error: "+err.getMessage());}//Element Name
		position+=4;
		for (int i =0;i<4;i++){eTagRead[i] = (short) (0<<8|swapBytes[i]);}
		while ((eTagRead[0] != elementTag[vars][0] || eTagRead[1] != elementTag[vars][1] || eTagRead[2] != elementTag[vars][2]  || eTagRead[3] != elementTag[vars][3]) //As long as the tag does not match
		&& (((long) (eTagRead[1])<<24 | (long) (eTagRead[0])<<16 | (long) (eTagRead[3])<<8 | (long) (eTagRead[2])) < ((long) (elementTag[vars][1])<<24 | (long) (elementTag[vars][0])<<16 | (long) (elementTag[vars][3])<<8 | (long) (elementTag[vars][2])))		 //&&//As long as the tag in file is smaller than the tag we are looking for
		){ 
			try{inFile.read(swapBytes,0,4);}catch (Exception err){System.err.println("Error: "+err.getMessage());} //Element length
			position+=4;
			elementLength = (long) ((0x000000ff &(int) swapBytes[3])<<24 | (0x000000ff &(int) swapBytes[2])<<16 | (0x000000ff &(int) swapBytes[1])<<8 | (0x000000ff &(int) swapBytes[0]));
			if (elementLength > 0){// && (int) (eTagRead[3]<<8 | eTagRead[2]) > 0 ){ //Don't skip over whole sections..
				if(elementLength > 500){buffer = new byte[(int) elementLength];}
				try{inFile.read(buffer,0,(int) elementLength);}catch (Exception err){System.err.println("Error: "+err.getMessage());}	//Read Element
				position+=elementLength;
			}
			try{inFile.mark(4);}catch (Exception err){System.err.println("Error Mark: "+err.getMessage());}	//Mark the file in case we need to rewind...
			try{inFile.read(swapBytes,0,4);}catch (Exception err){System.err.println("Error: "+err.getMessage());}//Element Name
			position+=4;
			for (int i =0;i<4;i++){eTagRead[i] = (short) (0x000000FF & (int)swapBytes[i]);}
		}
				
		if(eTagRead[0] == elementTag[vars][0] && eTagRead[1] == elementTag[vars][1] && eTagRead[2] == elementTag[vars][2]  && eTagRead[3] == elementTag[vars][3]){
			try{inFile.read(swapBytes,0,4);}catch (Exception err){System.err.println("Error: "+err.getMessage());} //Element length
			position+=4;
			elementLength = (long) ((0x000000ff &(int) swapBytes[3])<<24 | (0x000000ff &(int) swapBytes[2])<<16 | (0x000000ff &(int) swapBytes[1])<<8 | (0x000000ff &(int) swapBytes[0]));
			if (elementLength > 0){
					if (reserveMemory==1){
						returnValue = new byte[(int) elementLength+1];		//0008 0022
					}
				try{inFile.read(returnValue,0,(int) elementLength);}catch (Exception err){System.err.println("Error: "+err.getMessage());}	//Read Element
				position+=elementLength;
			}
		}else{
			try{inFile.reset();}catch (Exception err){System.err.println("Error Reset1: "+err.getMessage());}	//Rewind the file to the beginning of the last Element Name
		}
		return returnValue;
	}
	
	/*Implicit memory function*/
	void memorySearchTagImplicit(byte[] fileData,short[][] elementTag,int vars){
		byte[] swapBytes = new byte[4];
		long	elementLengthVariant=0;
		for (int i = 0; i<4;i++){swapBytes[i]=fileData[(int)position+i];}	//Element Name
		position+=4;
		for (int i =0;i<4;i++){eTagRead[i] = (short) (0x000000FF & (int)swapBytes[i]);}
		while ((eTagRead[0] != elementTag[vars][0] || eTagRead[1] != elementTag[vars][1] || eTagRead[2] != elementTag[vars][2]  || eTagRead[3] != elementTag[vars][3]) //As long as the tag does not match
		&& (((long) (eTagRead[1])<<24 | (long) (eTagRead[0])<<16 | (long) (eTagRead[3])<<8 | (long) (eTagRead[2])) < ((long) (elementTag[vars][1])<<24 | (long) (elementTag[vars][0])<<16 | (long) (elementTag[vars][3])<<8 | (long) (elementTag[vars][2])))
		 //&&//As long as the tag in file is smaller than the tag we are looking for
		//ftell(inFile) < fileLength
		){ 
			for (int i = 0; i<4;i++){swapBytes[i]=fileData[(int)position+i];}	//Element length variant
			position+=4;
			elementLengthVariant = (long) ((0x000000ff &(int) swapBytes[3])<<24 | (0x000000ff &(int) swapBytes[2])<<16 | (0x000000ff &(int) swapBytes[1])<<8 | (0x000000ff &(int) swapBytes[0]));
			if (elementLengthVariant > 0){
				position+=elementLengthVariant; //Skip element
			}
			for (int i = 0; i<4;i++){swapBytes[i]=fileData[(int)position+i];}	//Element Name
			position+=4;
			for (int i =0;i<4;i++){eTagRead[i] = (short) (0x000000FF & (int)swapBytes[i]);}
		}
	}

	int searchDicomTagMemoryImplicit(byte[] fileData,short[][] elementTag,int vars){
		int returnValue=0;
		byte[] buffer = new byte[500];
		byte[] swapBytes = new byte[4];
		short[] swapShort = new short[4];
		long	elementLengthVariant=0;
		//Read Data Section  0008
		//Data elements are in increasing order -> read the file until finding what we're looking for
		//Look for Acquisition date...
		memorySearchTagImplicit(fileData,elementTag,vars);
		if(eTagRead[0] == elementTag[vars][0] && eTagRead[1] == elementTag[vars][1] && eTagRead[2] == elementTag[vars][2]  && eTagRead[3] == elementTag[vars][3]){
			for (int i = 0; i<4;i++){swapBytes[i]=fileData[(int)position+i];}	//Element length variant
			position+=4;
			elementLengthVariant = (long) ((0x000000ff &(int) swapBytes[3])<<24 | (0x000000ff &(int) swapBytes[2])<<16 | (0x000000ff &(int) swapBytes[1])<<8 | (0x000000ff &(int) swapBytes[0]));
			if (elementLengthVariant > 0){
				for (int i = 0; i<(int) elementLengthVariant;i++){swapBytes[i]=fileData[(int)position+i];}	//Read Element
				position+=elementLengthVariant;
				returnValue = (int) ((0x000000ff &(int) swapBytes[1])<<8 | (0x000000ff &(int) swapBytes[0]));
			}
		}else{
			position-=4;
		}
		return returnValue;
	}
	
	byte[] searchDicomTagMemoryImplicit(byte[] fileData,short[][] elementTag, int vars, int reserveMemory){
		byte[] returnValue =new byte[0];
		byte[] buffer = new byte[500];
		byte[] swapBytes = new byte[4];
		long	elementLengthVariant=0;
		//Read Data Section  0008
		//Data elements are in increasing order -> read the file until finding what we're looking for
		//Look for Acquisition date...
		//System.out.println("Starting search "+ elementTag[vars][0]+" "+elementTag[vars][1]+" "+elementTag[vars][2]+" "+elementTag[vars][3]+" pos "+position);
		memorySearchTagImplicit(fileData,elementTag,vars);
		//System.out.println("Done MemorySearching "+ eTagRead[0]+" "+eTagRead[1]+" "+eTagRead[2]+" "+eTagRead[3]+ " pos "+ position);
		if(eTagRead[0] == elementTag[vars][0] && eTagRead[1] == elementTag[vars][1] && eTagRead[2] == elementTag[vars][2]  && eTagRead[3] == elementTag[vars][3]){
			for (int i = 0; i<4;i++){swapBytes[i]=fileData[(int)position+i];}	//Element length variant
			position+=4;
			elementLengthVariant = (long) ((0x000000ff &(int) swapBytes[3])<<24 | (0x000000ff &(int) swapBytes[2])<<16 | (0x000000ff &(int) swapBytes[1])<<8 | (0x000000ff &(int) swapBytes[0]));
			if (elementLengthVariant > 0){
				if (reserveMemory==1){
					returnValue = new byte[(int) elementLengthVariant+1];		//0008 0022
				}
				for (int i = 0; i<(int) elementLengthVariant;i++){returnValue[i]=fileData[(int)position+i];}	//Read Element
				position+=elementLengthVariant;
			}
		}else{
			position-=4;	//Rewind the file to the beginning of the last Element Name
		}
		return returnValue;
	}

	void readDataMemoryImplicit(byte[] fileData,short[][] elementTag,int vars, short[] data){

		byte[] buffer = new byte[500];
		byte[] swapBytes = new byte[4];
		long	elementLengthVariant=0;
		memorySearchTagImplicit(fileData,elementTag,vars); /*Search to the proper position*/
		//Read Data Section  0008
		//Data elements are in increasing order -> read the file until finding what we're looking for
		//Look for Acquisition date...
		
		if(eTagRead[0] == elementTag[vars][0] && eTagRead[1] == elementTag[vars][1] && eTagRead[2] == elementTag[vars][2]  && eTagRead[3] == elementTag[vars][3]){
			for (int i = 0; i<4;i++){swapBytes[i]=fileData[(int)position+i];}	//Element length variant
			position+=4;
			elementLengthVariant = (long) ((0x000000ff &(int) swapBytes[3])<<24 | (0x000000ff &(int) swapBytes[2])<<16 | (0x000000ff &(int) swapBytes[1])<<8 | (0x000000ff &(int) swapBytes[0]));
			if (elementLengthVariant > 0){
				for (int i = 0; i<elementLengthVariant/2;i++){data[i]=(short) ((0x000000ff & (int) fileData[(int)position+i*2+1])<<8 | (0x000000ff & (int) fileData[(int)position+i*2]));}	//Element length variant
				position+=elementLengthVariant;
			}
		}
	}
	/*Implicit memory function*/
	
	
	void readDataImplicit(DataInputStream inFile,short[][] elementTag,int vars, short[] data){

		byte[] buffer = new byte[500];
		byte[] swapBytes = new byte[4];
		long elementLength=0;
		//Read Data Section  0008
		//Data elements are in increasing order -> read the file until finding what we're looking for
		//Look for Acquisition date...
		try{inFile.read(swapBytes,0,4);}catch (Exception err){System.err.println("Error: "+err.getMessage());}//Element Name
		position+=4;
		for (int i =0;i<4;i++){eTagRead[i] = (short) (0x000000FF & (int)swapBytes[i]);}
		while ((eTagRead[0] != elementTag[vars][0] || eTagRead[1] != elementTag[vars][1] || eTagRead[2] != elementTag[vars][2]  || eTagRead[3] != elementTag[vars][3]) //As long as the tag does not match
		&& (((long) (eTagRead[1])<<24 | (long) (eTagRead[0])<<16 | (long) (eTagRead[3])<<8 | (long) (eTagRead[2])) < ((long) (elementTag[vars][1])<<24 | (long) (elementTag[vars][0])<<16 | (long) (elementTag[vars][3])<<8 | (long) (elementTag[vars][2])))		 //&&//As long as the tag in file is smaller than the tag we are looking for
		){ 
			try{inFile.read(swapBytes,0,4);}catch (Exception err){System.err.println("Error: "+err.getMessage());} //Element length
			position+=4;
			elementLength = (long) ((0x000000ff &(int) swapBytes[3])<<24 | (0x000000ff &(int) swapBytes[2])<<16 | (0x000000ff &(int) swapBytes[1])<<8 | (0x000000ff &(int) swapBytes[0]));
			if (elementLength > 0){// && (int) (eTagRead[3]<<8 | eTagRead[2]) > 0 ){ //Don't skip over whole sections..
				if(elementLength > 500){buffer = new byte[(int) elementLength];}
				try{inFile.read(buffer,0,(int) elementLength);}catch (Exception err){System.err.println("Error: "+err.getMessage());}	//Read Element
				position+=elementLength;
			}
			try{inFile.mark(4);}catch (Exception err){System.err.println("Error Mark: "+err.getMessage());}	//Mark the file in case we need to rewind...
			try{inFile.read(swapBytes,0,4);}catch (Exception err){System.err.println("Error: "+err.getMessage());}//Element Name
			position+=4;
			for (int i =0;i<4;i++){eTagRead[i] = (short) (0x000000FF & (int)swapBytes[i]);}
		}
		if(eTagRead[0] == elementTag[vars][0] && eTagRead[1] == elementTag[vars][1] && eTagRead[2] == elementTag[vars][2]  && eTagRead[3] == elementTag[vars][3]){
try{inFile.read(swapBytes,0,4);}catch (Exception err){System.err.println("Error: "+err.getMessage());} //Element length
			position+=4;
			elementLength = (long) ((0x000000ff &(int) swapBytes[3])<<24 | (0x000000ff &(int) swapBytes[2])<<16 | (0x000000ff &(int) swapBytes[1])<<8 | (0x000000ff &(int) swapBytes[0]));
			if (elementLength > 0){
				//Find out whether it is unsigned or signed short or integer....
				for (int i=0;i<(int) elementLength/2;i++){
					try{data[i]=Short.reverseBytes(inFile.readShort());}catch (Exception err){System.err.println("Error: "+err.getMessage());} //Read Element
					position+=2;
				}
			}
		}else{
			try{inFile.reset();}catch (Exception err){System.err.println("Error Reset1: "+err.getMessage());}	//Rewind the file to the beginning of the last Element Name
		}
	}
	/*IMPLICIT VR ENDS*/
	
	/*EXPLICIT VR*/
	void readExplicitMemory(){
	
		//Data elements are in increasing order -> read the file until finding what we're looking for
		++vars;
		acquisitionDate = searchDicomTagMemory(fileData,elementTag,vars,1);
		//System.out.println("aDate "+new String(acquisitionDate));
		
		++vars;
		studyDescription = searchDicomTagMemory(fileData,elementTag,vars,1);
		++vars;
		seriesDescription = searchDicomTagMemory(fileData,elementTag,vars,1);
		++vars;
		patientName = searchDicomTagMemory(fileData,elementTag,vars,1);
		++vars;
		patientID = searchDicomTagMemory(fileData,elementTag,vars,1);
		++vars;
		patientBirthDate = searchDicomTagMemory(fileData,elementTag,vars,1);
		//System.out.println("pBirth "+new String(patientBirthDate));
		++vars;
		imageComments = searchDicomTagMemory(fileData,elementTag,vars,1);
		
		++vars;	
		samplesPerPixel = searchDicomTagMemory(fileData,elementTag,vars);
		++vars;
		rows = searchDicomTagMemory(fileData,elementTag,vars);
		++vars;
		columns = searchDicomTagMemory(fileData,elementTag,vars);
		++vars;
		pixelSpacing = searchDicomTagMemory(fileData,elementTag,vars,1);
		++vars;
		bitsAllocated = searchDicomTagMemory(fileData,elementTag,vars);
		++vars;
		bitsStored = searchDicomTagMemory(fileData,elementTag,vars);
		++vars;
		highBit = searchDicomTagMemory(fileData,elementTag,vars);
		++vars;
		pixelRepresentation = searchDicomTagMemory(fileData,elementTag,vars);
		++vars;
		rescaleIntercept = searchDicomTagMemory(fileData,elementTag,vars,1);
		++vars;
		rescaleSlope = searchDicomTagMemory(fileData,elementTag,vars,1);
		++vars;
		//Header info read, read DATA
		
		data = new short[rows*columns];
		
		//System.out.println("Dataa lukemaan "+rows*columns);
		readDataMemory(fileData,elementTag,vars, data);
		short[] temp = (short[]) data.clone();
		Arrays.sort(temp);
		max = temp[temp.length-1];
		min = temp[0];
	}

	void readExplicit(){	
		//Data elements are in increasing order -> read the file until finding what we're looking for
		++vars;
		acquisitionDate = searchDicomTag(inFile,elementTag,vars,1);
		//System.out.println("aDate "+new String(acquisitionDate));
		
		++vars;
		studyDescription = searchDicomTag(inFile,elementTag,vars,1);
		++vars;
		seriesDescription = searchDicomTag(inFile,elementTag,vars,1);
		++vars;
		patientName = searchDicomTag(inFile,elementTag,vars,1);
		++vars;
		patientID = searchDicomTag(inFile,elementTag,vars,1);
		++vars;
		patientBirthDate = searchDicomTag(inFile,elementTag,vars,1);
		//System.out.println("pBirth "+new String(patientBirthDate));
		++vars;
		imageComments = searchDicomTag(inFile,elementTag,vars,1);
		
		++vars;	
		samplesPerPixel = searchDicomTag(inFile,elementTag,vars);
		++vars;
		rows = searchDicomTag(inFile,elementTag,vars);
		++vars;
		columns = searchDicomTag(inFile,elementTag,vars);
		++vars;
		pixelSpacing = searchDicomTag(inFile,elementTag,vars,1);
		++vars;
		bitsAllocated = searchDicomTag(inFile,elementTag,vars);
		++vars;
		bitsStored = searchDicomTag(inFile,elementTag,vars);
		++vars;
		highBit = searchDicomTag(inFile,elementTag,vars);
		++vars;
		pixelRepresentation = searchDicomTag(inFile,elementTag,vars);
		++vars;
		rescaleIntercept = searchDicomTag(inFile,elementTag,vars,1);
		++vars;
		rescaleSlope = searchDicomTag(inFile,elementTag,vars,1);
		++vars;
		//Header info read, read DATA
		
		data = new short[rows*columns];
		
		//System.out.println("Dataa lukemaan "+rows*columns);
		readData(inFile,elementTag,vars, data);
		short[] temp = (short[]) data.clone();
		Arrays.sort(temp);
		max = temp[temp.length-1];
		min = temp[0];
	}	

	byte[] searchHeaderTag(DataInputStream inFile,short[][] elementTag, int vars, int reserveMemory){
		byte[] returnValue =new byte[0];
		byte[] buffer = new byte[500];
		byte[] swapBytes = new byte[4];
		byte VR[] = new byte[2];
		String stringVR;
		int elementLength=0;
		long	elementLengthVariant=0;
		//Read Data Section  0008
		//Data elements are in increasing order -> read the file until finding what we're looking for
		//Look for Acquisition date...
		try{inFile.read(swapBytes,0,4);}catch (Exception err){System.err.println("Error: "+err.getMessage());}//Element Name
		position+=4;
		for (int i =0;i<4;i++){eTagRead[i] = (short) (0<<8|swapBytes[i]);}
		while ((eTagRead[0] != elementTag[vars][0] || eTagRead[1] != elementTag[vars][1] || eTagRead[2] != elementTag[vars][2]  || eTagRead[3] != elementTag[vars][3]) //As long as the tag does not match
		&& (((long) (eTagRead[1])<<24 | (long) (eTagRead[0])<<16 | (long) (eTagRead[3])<<8 | (long) (eTagRead[2])) < ((long) (elementTag[vars][1])<<24 | (long) (elementTag[vars][0])<<16 | (long) (elementTag[vars][3])<<8 | (long) (elementTag[vars][2])))		 //&&//As long as the tag in file is smaller than the tag we are looking for
		){ 
			try{inFile.read(VR,0,2);}catch (Exception err){System.err.println("Error: "+err.getMessage());}//Value representation (OB OW OF SQ UT & UN different representation)
			position+=2;
			//Check whether VR is one of the special ones
			stringVR = new String(VR);
			if (stringVR.equals(new String("OB")) == false &&	stringVR.equals(new String("OW")) == false &&	stringVR.equals(new String("OF")) == false &&	stringVR.equals(new String("SQ")) == false &&	stringVR.equals(new String("UT")) == false &&	stringVR.equals(new String("UN")) == false){
				try{inFile.read(swapBytes,0,2);}catch (Exception err){System.err.println("Error: "+err.getMessage());} //Element length
				position+=2;
				elementLength = (int) swapBytes[1]<<8|swapBytes[0];
				if (elementLength > 0){
					if(elementLength > 500){buffer = new byte[elementLength];}
					try{inFile.read(buffer,0,elementLength);}catch (Exception err){System.err.println("Error: "+err.getMessage());}	//Read Element
					position+=elementLength;
				}
			} else {
				try{elementLength = Short.reverseBytes(inFile.readShort());}catch (Exception err){System.err.println("Error: "+err.getMessage());} //Element length 0x00
				position+=2;
				try{inFile.read(swapBytes,0,4);}catch (Exception err){System.err.println("Error: "+err.getMessage());} //Element length
				position+=4;
				elementLengthVariant = (long) ((0x000000ff &(int) swapBytes[3])<<24 | (0x000000ff &(int) swapBytes[2])<<16 | (0x000000ff &(int) swapBytes[1])<<8 | (0x000000ff &(int) swapBytes[0]));
				if(elementLengthVariant > 500){buffer = new byte[(int)elementLengthVariant];}
				try{inFile.read(buffer,0,(int) elementLengthVariant);}catch (Exception err){System.err.println("Error: "+err.getMessage());}	//Read Element
				position+=elementLengthVariant;
			}

			try{inFile.read(swapBytes,0,4);}catch (Exception err){System.err.println("Error: "+err.getMessage());}//Element Name
			position+=4;
			for (int i =0;i<4;i++){eTagRead[i] = (short) (0x000000FF & (int)swapBytes[i]);}
		}
				
		if(eTagRead[0] == elementTag[vars][0] && eTagRead[1] == elementTag[vars][1] && eTagRead[2] == elementTag[vars][2]  && eTagRead[3] == elementTag[vars][3]){
			try{inFile.read(VR,0,2);}catch (Exception err){System.err.println("Error: "+err.getMessage());}//Value representation (OB OW OF SQ UT & UN different representation)
				position+=2;
				stringVR = new String(VR);
				if (stringVR.equals(new String("OB")) == false &&	stringVR.equals(new String("OW")) == false &&	stringVR.equals(new String("OF")) == false &&	stringVR.equals(new String("SQ")) == false &&	stringVR.equals(new String("UT")) == false &&	stringVR.equals(new String("UN")) == false){
					try{inFile.read(swapBytes,0,2);}catch (Exception err){System.err.println("Error: "+err.getMessage());} //Element length
					position+=2;
					elementLength = (int) swapBytes[1]<<8|swapBytes[0];
					if (elementLength > 0){
						if (reserveMemory==1){
							returnValue = new byte[elementLength+1];		//0008 0022
						}
						//Find out whether it is unsigned or signed short or integer....
						try{inFile.read(returnValue,0,elementLength);}catch (Exception err){System.err.println("Error: "+err.getMessage());} //Read Element
						position+=elementLength;
					}
				} else {
					try{elementLength = Short.reverseBytes(inFile.readShort());}catch (Exception err){System.err.println("Error: "+err.getMessage());} //Element length 0x00
					position+=2;
					try{inFile.read(swapBytes,0,4);}catch (Exception err){System.err.println("Error: "+err.getMessage());} //Element length
					position+=4;
					elementLengthVariant = (long) ((0x000000ff &(int) swapBytes[3])<<24 | (0x000000ff &(int) swapBytes[2])<<16 | (0x000000ff &(int) swapBytes[1])<<8 | (0x000000ff &(int) swapBytes[0]));
					if (reserveMemory==1){
						returnValue = new byte[(int) elementLengthVariant+1];		//0008 0022
					}
					try{inFile.read(returnValue,0,(int) elementLengthVariant);}catch (Exception err){System.err.println("Error: "+err.getMessage());} //Read Element
					position+=elementLengthVariant;
				}
		}
		return returnValue;
	}	
	
	int searchDicomTag(DataInputStream inFile,short[][] elementTag,int vars){
		int returnValue=0;
		byte[] buffer = new byte[500];
		byte[] swapBytes = new byte[4];
		byte VR[] = new byte[2];
		String stringVR;
		int elementLength=0;
		long	elementLengthVariant=0;
		//Read Data Section  0008
		//Data elements are in increasing order -> read the file until finding what we're looking for
		//Look for Acquisition date...
		try{inFile.read(swapBytes,0,4);}catch (Exception err){System.err.println("Error: "+err.getMessage());}//Element Name
		position+=4;
		for (int i =0;i<4;i++){eTagRead[i] = (short) (0x000000FF & (int)swapBytes[i]);}
		while ((eTagRead[0] != elementTag[vars][0] || eTagRead[1] != elementTag[vars][1] || eTagRead[2] != elementTag[vars][2]  || eTagRead[3] != elementTag[vars][3]) //As long as the tag does not match
		&& (((long) (eTagRead[1])<<24 | (long) (eTagRead[0])<<16 | (long) (eTagRead[3])<<8 | (long) (eTagRead[2])) < ((long) (elementTag[vars][1])<<24 | (long) (elementTag[vars][0])<<16 | (long) (elementTag[vars][3])<<8 | (long) (elementTag[vars][2])))		 //&&//As long as the tag in file is smaller than the tag we are looking for
		){ 
			try{inFile.read(VR,0,2);}catch (Exception err){System.err.println("Error: "+err.getMessage());}//Value representation (OB OW OF SQ UT & UN different representation)
			position+=2;
			//Check whether VR is one of the special ones
			stringVR = new String(VR);
			if (stringVR.equals(new String("OB")) == false &&	stringVR.equals(new String("OW")) == false &&	stringVR.equals(new String("OF")) == false &&	stringVR.equals(new String("SQ")) == false &&	stringVR.equals(new String("UT")) == false &&	stringVR.equals(new String("UN")) == false){
				try{inFile.read(swapBytes,0,2);}catch (Exception err){System.err.println("Error: "+err.getMessage());} //Element length
				position+=2;
				elementLength = (int) swapBytes[1]<<8|swapBytes[0];
				if (elementLength > 0){
					if(elementLength > 500){buffer = new byte[elementLength];}
					try{inFile.read(buffer,0,elementLength);}catch (Exception err){System.err.println("Error: "+err.getMessage());}	//Read Element
					position+=elementLength;
				}
			} else {
				try{elementLength = Short.reverseBytes(inFile.readShort());}catch (Exception err){System.err.println("Error: "+err.getMessage());} //Element length 0x00
				position+=2;
				try{inFile.read(swapBytes,0,4);}catch (Exception err){System.err.println("Error: "+err.getMessage());} //Element length
				position+=4;
				elementLengthVariant = (long) ((0x000000ff &(int) swapBytes[3])<<24 | (0x000000ff &(int) swapBytes[2])<<16 | (0x000000ff &(int) swapBytes[1])<<8 | (0x000000ff &(int) swapBytes[0]));
				if(elementLengthVariant > 500){buffer = new byte[(int)elementLengthVariant];}
				try{inFile.read(buffer,0,(int) elementLengthVariant);}catch (Exception err){System.err.println("Error: "+err.getMessage());}	//Read Element
				position+=elementLengthVariant;
			}
			
			try{inFile.mark(4);}catch (Exception err){System.err.println("Error Mark: "+err.getMessage());}	//Mark the file in case we need to rewind...
			try{inFile.read(swapBytes,0,4);}catch (Exception err){System.err.println("Error: "+err.getMessage());}//Element Name
			position+=4;
			for (int i =0;i<4;i++){eTagRead[i] = (short) (0x000000FF & (int)swapBytes[i]);}
		}
				
		if(eTagRead[0] == elementTag[vars][0] && eTagRead[1] == elementTag[vars][1] && eTagRead[2] == elementTag[vars][2]  && eTagRead[3] == elementTag[vars][3]){
			try{inFile.read(VR,0,2);}catch (Exception err){System.err.println("Error: "+err.getMessage());}//Value representation (OB OW OF SQ UT & UN different representation)
				position+=2;
				stringVR = new String(VR);
				if (stringVR.equals(new String("OB")) == false &&	stringVR.equals(new String("OW")) == false &&	stringVR.equals(new String("OF")) == false &&	stringVR.equals(new String("SQ")) == false &&	stringVR.equals(new String("UT")) == false &&	stringVR.equals(new String("UN")) == false){
					try{inFile.read(swapBytes,0,2);}catch (Exception err){System.err.println("Error: "+err.getMessage());} //Element length
					position+=2;
					elementLength = (int) swapBytes[1]<<8|swapBytes[0];
					if (elementLength > 0){
						//Find out whether it is unsigned or signed short or integer....
						try{inFile.read(swapBytes,0,elementLength);}catch (Exception err){System.err.println("Error: "+err.getMessage());} //Read Element
						position+=elementLength;
						returnValue = (int) swapBytes[1]<<8|swapBytes[0];
					}
				} else {
					try{elementLength = Short.reverseBytes(inFile.readShort());}catch (Exception err){System.err.println("Error: "+err.getMessage());} //Element length 0x00
					try{inFile.read(swapBytes,0,4);}catch (Exception err){System.err.println("Error: "+err.getMessage());} //Element length
					position+=4;
					elementLengthVariant = (long) ((0x000000ff &(int) swapBytes[3])<<24 | (0x000000ff &(int) swapBytes[2])<<16 | (0x000000ff &(int) swapBytes[1])<<8 | (0x000000ff &(int) swapBytes[0]));
					try{inFile.read(swapBytes,0,(int) elementLengthVariant);}catch (Exception err){System.err.println("Error: "+err.getMessage());} //Read Element
					position+=elementLengthVariant;
					returnValue = (int) swapBytes[1]<<8|swapBytes[0];
				}
		}else{
			try{inFile.reset();}catch (Exception err){System.err.println("Error Reset1: "+err.getMessage());}	//Rewind the file to the beginning of the last Element Name
		}
		return returnValue;
	}
	
	byte[] searchDicomTag(DataInputStream inFile,short[][] elementTag, int vars, int reserveMemory){
		byte[] returnValue =new byte[0];
		byte[] buffer = new byte[500];
		byte[] swapBytes = new byte[4];
		byte VR[] = new byte[2];
		String stringVR;
		int elementLength=0;
		long	elementLengthVariant=0;
		//Read Data Section  0008
		//Data elements are in increasing order -> read the file until finding what we're looking for
		//Look for Acquisition date...
		try{inFile.read(swapBytes,0,4);}catch (Exception err){System.err.println("Error: "+err.getMessage());}//Element Name
		position+=4;
		for (int i =0;i<4;i++){eTagRead[i] = (short) (0<<8|swapBytes[i]);}
		while ((eTagRead[0] != elementTag[vars][0] || eTagRead[1] != elementTag[vars][1] || eTagRead[2] != elementTag[vars][2]  || eTagRead[3] != elementTag[vars][3]) //As long as the tag does not match
		&& (((long) (eTagRead[1])<<24 | (long) (eTagRead[0])<<16 | (long) (eTagRead[3])<<8 | (long) (eTagRead[2])) < ((long) (elementTag[vars][1])<<24 | (long) (elementTag[vars][0])<<16 | (long) (elementTag[vars][3])<<8 | (long) (elementTag[vars][2])))		 //&&//As long as the tag in file is smaller than the tag we are looking for
		){ 
			try{inFile.read(VR,0,2);}catch (Exception err){System.err.println("Error: "+err.getMessage());}//Value representation (OB OW OF SQ UT & UN different representation)
			position+=2;
			//Check whether VR is one of the special ones
			stringVR = new String(VR);
			if (stringVR.equals(new String("OB")) == false &&	stringVR.equals(new String("OW")) == false &&	stringVR.equals(new String("OF")) == false &&	stringVR.equals(new String("SQ")) == false &&	stringVR.equals(new String("UT")) == false &&	stringVR.equals(new String("UN")) == false){
				try{inFile.read(swapBytes,0,2);}catch (Exception err){System.err.println("Error: "+err.getMessage());} //Element length
				position+=2;
				elementLength = (int) swapBytes[1]<<8|swapBytes[0];
				if (elementLength > 0){
					if(elementLength > 500){buffer = new byte[elementLength];}
					try{inFile.read(buffer,0,elementLength);}catch (Exception err){System.err.println("Error: "+err.getMessage());}	//Read Element
					position+=elementLength;
				}
			} else {
				try{elementLength = Short.reverseBytes(inFile.readShort());}catch (Exception err){System.err.println("Error: "+err.getMessage());} //Element length 0x00
				position+=2;
				try{inFile.read(swapBytes,0,4);}catch (Exception err){System.err.println("Error: "+err.getMessage());} //Element length
				position+=4;
				elementLengthVariant = (long) ((0x000000ff &(int) swapBytes[3])<<24 | (0x000000ff &(int) swapBytes[2])<<16 | (0x000000ff &(int) swapBytes[1])<<8 | (0x000000ff &(int) swapBytes[0]));
				if(elementLengthVariant > 500){buffer = new byte[(int)elementLengthVariant];}
				try{inFile.read(buffer,0,(int) elementLengthVariant);}catch (Exception err){System.err.println("Error: "+err.getMessage());}	//Read Element
				position+=elementLengthVariant;
			}
			
			try{inFile.mark(4);}catch (Exception err){System.err.println("Error Mark: "+err.getMessage());}	//Mark the file in case we need to rewind...
			try{inFile.read(swapBytes,0,4);}catch (Exception err){System.err.println("Error: "+err.getMessage());}//Element Name
			position+=4;
			for (int i =0;i<4;i++){eTagRead[i] = (short) (0x000000FF & (int)swapBytes[i]);}
		}
				
		if(eTagRead[0] == elementTag[vars][0] && eTagRead[1] == elementTag[vars][1] && eTagRead[2] == elementTag[vars][2]  && eTagRead[3] == elementTag[vars][3]){
			try{inFile.read(VR,0,2);}catch (Exception err){System.err.println("Error: "+err.getMessage());}//Value representation (OB OW OF SQ UT & UN different representation)
				position+=2;
				stringVR = new String(VR);
				if (stringVR.equals(new String("OB")) == false &&	stringVR.equals(new String("OW")) == false &&	stringVR.equals(new String("OF")) == false &&	stringVR.equals(new String("SQ")) == false &&	stringVR.equals(new String("UT")) == false &&	stringVR.equals(new String("UN")) == false){
					try{inFile.read(swapBytes,0,2);}catch (Exception err){System.err.println("Error: "+err.getMessage());} //Element length
					position+=2;
					elementLength = (int) swapBytes[1]<<8|swapBytes[0];
					if (elementLength > 0){
						if (reserveMemory==1){
							returnValue = new byte[elementLength+1];		//0008 0022
						}
						//Find out whether it is unsigned or signed short or integer....
						try{inFile.read(returnValue,0,elementLength);}catch (Exception err){System.err.println("Error: "+err.getMessage());} //Read Element
						position+=elementLength;
					}
				} else {
					try{elementLength = Short.reverseBytes(inFile.readShort());}catch (Exception err){System.err.println("Error: "+err.getMessage());} //Element length 0x00
					position+=2;
					try{inFile.read(swapBytes,0,4);}catch (Exception err){System.err.println("Error: "+err.getMessage());} //Element length
					position+=4;
					elementLengthVariant = (long) ((0x000000ff &(int) swapBytes[3])<<24 | (0x000000ff &(int) swapBytes[2])<<16 | (0x000000ff &(int) swapBytes[1])<<8 | (0x000000ff &(int) swapBytes[0]));
					if (reserveMemory==1){
						returnValue = new byte[(int) elementLengthVariant+1];		//0008 0022
					}
					try{inFile.read(returnValue,0,(int) elementLengthVariant);}catch (Exception err){System.err.println("Error: "+err.getMessage());} //Read Element
					position+=elementLengthVariant;
				}
		}else{
			try{inFile.reset();}catch (Exception err){System.err.println("Error Reset1: "+err.getMessage());}	//Rewind the file to the beginning of the last Element Name
		}
		return returnValue;
	}
	
	void readData(DataInputStream inFile,short[][] elementTag,int vars, short[] data){

		byte[] buffer = new byte[500];
		byte[] swapBytes = new byte[4];
		byte VR[] = new byte[2];
		String stringVR;
		int elementLength=0;
		long	elementLengthVariant=0;
		//Read Data Section  0008
		//Data elements are in increasing order -> read the file until finding what we're looking for
		//Look for Acquisition date...
		try{inFile.read(swapBytes,0,4);}catch (Exception err){System.err.println("Error: "+err.getMessage());}//Element Name
		position+=4;
		for (int i =0;i<4;i++){eTagRead[i] = (short) (0x000000FF & (int)swapBytes[i]);}
		while ((eTagRead[0] != elementTag[vars][0] || eTagRead[1] != elementTag[vars][1] || eTagRead[2] != elementTag[vars][2]  || eTagRead[3] != elementTag[vars][3]) //As long as the tag does not match
		&& (((long) (eTagRead[1])<<24 | (long) (eTagRead[0])<<16 | (long) (eTagRead[3])<<8 | (long) (eTagRead[2])) < ((long) (elementTag[vars][1])<<24 | (long) (elementTag[vars][0])<<16 | (long) (elementTag[vars][3])<<8 | (long) (elementTag[vars][2])))		 //&&//As long as the tag in file is smaller than the tag we are looking for
		){ 
			try{inFile.read(VR,0,2);}catch (Exception err){System.err.println("Error: "+err.getMessage());}//Value representation (OB OW OF SQ UT & UN different representation)
			position+=2;
			//Check whether VR is one of the special ones
			stringVR = new String(VR);
			if (stringVR.equals(new String("OB")) == false &&	stringVR.equals(new String("OW")) == false &&	stringVR.equals(new String("OF")) == false &&	stringVR.equals(new String("SQ")) == false &&	stringVR.equals(new String("UT")) == false &&	stringVR.equals(new String("UN")) == false){
				try{inFile.read(swapBytes,0,2);}catch (Exception err){System.err.println("Error: "+err.getMessage());} //Element length
				position+=2;
				elementLength = (int) swapBytes[1]<<8|swapBytes[0];
				//try{elementLength = Short.reverseBytes(inFile.readUnsignedShort());}catch (Exception err){System.err.println("Error: "+err.getMessage());} //Element length
				if (elementLength > 0){
					if(elementLength > 500){buffer = new byte[elementLength];}
					try{inFile.read(buffer,0,elementLength);}catch (Exception err){System.err.println("Error: "+err.getMessage());}	//Read Element
					position+=elementLength;
				}
			} else {
				////System.out.println("VAR "+stringVR);
				try{elementLength = Short.reverseBytes(inFile.readShort());}catch (Exception err){System.err.println("Error: "+err.getMessage());} //Element length 0x00
				position+=2;
				try{inFile.read(swapBytes,0,4);}catch (Exception err){System.err.println("Error: "+err.getMessage());} //Element length
				position+=4;
				elementLengthVariant = (long) ((0x000000ff &(int) swapBytes[3])<<24 | (0x000000ff &(int) swapBytes[2])<<16 | (0x000000ff &(int) swapBytes[1])<<8 | (0x000000ff &(int) swapBytes[0]));
				if(elementLengthVariant > 500){buffer = new byte[(int)elementLengthVariant];}
				try{inFile.read(buffer,0,(int) elementLengthVariant);}catch (Exception err){System.err.println("Error: "+err.getMessage());}	//Read Element
				position+=elementLengthVariant;
			}
			
			try{inFile.mark(4);}catch (Exception err){System.err.println("Error Mark: "+err.getMessage());}	//Mark the file in case we need to rewind...
			try{inFile.read(swapBytes,0,4);}catch (Exception err){System.err.println("Error: "+err.getMessage());}//Element Name
			position+=4;
			for (int i =0;i<4;i++){eTagRead[i] = (short) (0x000000FF & (int)swapBytes[i]);}
		}
		if(eTagRead[0] == elementTag[vars][0] && eTagRead[1] == elementTag[vars][1] && eTagRead[2] == elementTag[vars][2]  && eTagRead[3] == elementTag[vars][3]){
			try{inFile.read(VR,0,2);}catch (Exception err){System.err.println("Error: "+err.getMessage());}//Value representation (OB OW OF SQ UT & UN different representation)
				position+=2;
				stringVR = new String(VR);
				if (stringVR.equals(new String("OB")) == false &&	stringVR.equals(new String("OW")) == false &&	stringVR.equals(new String("OF")) == false &&	stringVR.equals(new String("SQ")) == false &&	stringVR.equals(new String("UT")) == false &&	stringVR.equals(new String("UN")) == false){
					try{inFile.read(swapBytes,0,2);}catch (Exception err){System.err.println("Error: "+err.getMessage());} //Element length
					position+=2;
					elementLength = (int) swapBytes[1]<<8|swapBytes[0];
					if (elementLength > 0){
						//Find out whether it is unsigned or signed short or integer....
						for (int i=0;i<elementLength/2;i++){
							try{data[i]=Short.reverseBytes(inFile.readShort());}catch (Exception err){System.err.println("Error: "+err.getMessage());} //Read Element
							position+=2;
						}
					}
				} else {
					try{elementLength = Short.reverseBytes(inFile.readShort());}catch (Exception err){System.err.println("Error: "+err.getMessage());} //Element length 0x00
					position+=2;
					try{inFile.read(swapBytes,0,4);}catch (Exception err){System.err.println("Error: "+err.getMessage());} //Element length
					position+=4;
					elementLengthVariant = (long) ((0x000000ff &(int) swapBytes[3])<<24 | (0x000000ff &(int) swapBytes[2])<<16 | (0x000000ff &(int) swapBytes[1])<<8 | (0x000000ff &(int) swapBytes[0]));
					for (int i=0;i<elementLengthVariant/2;i++){
						try{data[i]=Short.reverseBytes(inFile.readShort());}catch (Exception err){System.err.println("Error: "+err.getMessage());} //Read Element
						position+=2;
					}
				}
		}else{
			try{inFile.reset();}catch (Exception err){System.err.println("Error Reset1: "+err.getMessage());}	//Rewind the file to the beginning of the last Element Name
		}
	}
	
	void readDataMemory(byte[] fileData,short[][] elementTag,int vars, short[] data){

		byte[] buffer = new byte[500];
		byte[] swapBytes = new byte[4];
		byte VR[] = new byte[2];
		String stringVR;
		int elementLength=0;
		long	elementLengthVariant=0;

		memorySearchTag(fileData,elementTag,vars); /*Search to the proper position*/
		//Read Data Section  0008
		//Data elements are in increasing order -> read the file until finding what we're looking for
		//Look for Acquisition date...
		
		if(eTagRead[0] == elementTag[vars][0] && eTagRead[1] == elementTag[vars][1] && eTagRead[2] == elementTag[vars][2]  && eTagRead[3] == elementTag[vars][3]){
		
			for (int i =0;i<2;i++){VR[i] = fileData[(int)position+i];} //Value representation
			position+=2;
			//Check whether VR is one of the special ones
			stringVR = new String(VR);
			//System.out.println(stringVR);
			if (stringVR.equals(new String("OB")) == false &&	stringVR.equals(new String("OW")) == false &&	stringVR.equals(new String("OF")) == false &&	stringVR.equals(new String("SQ")) == false &&	stringVR.equals(new String("UT")) == false &&	stringVR.equals(new String("UN")) == false){
				for (int i = 0; i<2;i++){swapBytes[i]=fileData[(int)position+i];}	//Element length
				position+=2;
				elementLength = (int) swapBytes[1]<<8|swapBytes[0];
				if (elementLength > 0){
					for (int i = 0; i<elementLength/2;i++){data[i]=(short) ((0x000000ff & (int) fileData[(int)position+i*2+1])<<8 | (0x000000ff & (int) fileData[(int)position+i*2]));}	//Element length variant
					position+=elementLength;
				}
			} else {
				position+=2;	//Skip empty element length
				for (int i = 0; i<4;i++){swapBytes[i]=fileData[(int)position+i];}	//Element length variant
				position+=4;
				elementLengthVariant = (long) ((0x000000ff &(int) swapBytes[3])<<24 | (0x000000ff &(int) swapBytes[2])<<16 | (0x000000ff &(int) swapBytes[1])<<8 | (0x000000ff &(int) swapBytes[0]));
				if (elementLengthVariant > 0){
					for (int i = 0; i<elementLengthVariant/2;i++){data[i]=(short) ((0x000000ff & (int) fileData[(int)position+i*2+1])<<8 | (0x000000ff & (int) fileData[(int)position+i*2]));}	//Element length variant
					position+=elementLengthVariant;
				}
			}
		
		}
	}
	
	/*EXPLICIT VR*/
}
