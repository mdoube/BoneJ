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

	ImageJ Stratec file header stripper. Files stripped with this plug-in still contain a header,
	just not data identifying the patient without having some additional knowledge.
    Copyright (C) 2012 Timo Rantalainen
*/

package org.doube.bonej.pqct;
import ij.plugin.PlugIn;
import ij.Prefs;
import java.io.*;
import java.util.*;
import ij.gui.*;
import ij.IJ;
import javax.activation.*;					//UnsupportedDataTypeException

import org.doube.util.UsageReporter;

public class Strip_Stratec_File_Header implements PlugIn {

	//Overriding the abstract runnable run method. Apparently plugins run in threads
	public void run(String arg) {
		GenericDialog dialog = new GenericDialog("Strip Stratec Header");
		/*MeasInfo,PatBirth,PatMenoAge,PatName,PatTitle&Comment*/
		dialog.addCheckbox("Strip_MeasInfo",false);
		dialog.addCheckbox("Strip_PatBirth",false);
		dialog.addCheckbox("Strip_PatMenoAge",false);
		dialog.addCheckbox("Strip_PatName",true);
		dialog.addCheckbox("Strip_PatTitleAndComment",false);
		dialog.addStringField("Stratec_file_to_strip",Prefs.getDefaultDirectory()+"I0020001.M01",60);
		dialog.addStringField("File_save_name",Prefs.getDefaultDirectory()+"I0020001.M01",60);
		dialog.showDialog();
		
		if (dialog.wasOKed()){ //Stop in case of cancel..
			boolean[] toStrip = new boolean[5];
			for (int i = 0; i<toStrip.length; ++i){
				toStrip[i]			= dialog.getNextBoolean();
			}
			String fileIn 		= dialog.getNextString();
			String fileOut 		= dialog.getNextString();
			if (fileIn==null || fileOut==null) {
				IJ.log("Give both input and output file as parameters");
				return;
			}
			File test = new File(fileIn);
			if (!test.exists()){
				IJ.error("Input file didn't exist");
				return;
			}
			
			try{
				stripFile(fileIn,fileOut,toStrip);
			}catch (Exception err){
				IJ.error("Stratec file header stripping failed", err.getMessage());
			}
		}
		UsageReporter.reportEvent(this).send();
	}

	private void stripFile(String fileNameIn,String fileNameOut,boolean[] toStrip) throws Exception {
		File fileIn = new File(fileNameIn);
		long fileLength = fileIn.length();
		byte[] fileData;
		try{
			BufferedInputStream inputStream = new BufferedInputStream(new FileInputStream(fileIn));
			DataInputStream dataInputStream = new DataInputStream(inputStream);
			fileData = new byte[(int) fileLength];		//Allocate memory for reading the file into memory
			dataInputStream.read(fileData,0,(int) fileLength);		//Read the data to memory
			dataInputStream.close();	//Close the file after reading
		}catch (Exception e) {
			throw new UnsupportedDataTypeException("Could not read input file.");
		}
		fileData = stripHeader(fileData,toStrip); //Strip header
		writeFile(fileNameOut,fileData);
	}
	
	
	public void writeFile(String fileName,byte[] fileData){
		try{
			FileOutputStream writer = new FileOutputStream(fileName);
			writer.write(fileData);
			writer.close();
		}catch (Exception err){System.out.println("Saving failed");}
	}
	
	//Writing dummy header containing sufficient details for Distribution_Analysis imageJ plugin, might not suffice for Geanie or Stractec software
	byte[] stripHeader(byte[] data,boolean[] toStrip){
		int[] offsetsToStrip = {662,1091,1095,1099,1141};	/*MeasInfo,PatBirth,PatMenoAge,PatName,PatTitle&Comment*/
		int[] stripLengths = {324,4,4,41,124};
		for (int s= 0;s<offsetsToStrip.length;++s){
			if (toStrip[s]){
				data = fillWithZero(data,offsetsToStrip[s],stripLengths[s]);
			}
		}
		return data;
	}
	byte[] fillWithZero(byte[] data,int offset,int zeros){
		for (int z = 0;z<zeros;++z){
			data[offset+z] = (byte) 0;
		}
		return data;
	}
}
