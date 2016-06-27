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

	ImageJ Stratec CT file header reader plugin
    Copyright (C) 2016 Timo Rantalainen
 */

package org.doube.bonej.pqct;

import java.io.*;
import javax.activation.*; //UnsupportedDataTypeException

public class CTHeaderReader{
	File fileIn;
	public double percent;
	public String site;
	CTHeaderReader(String fileName){
		fileIn = new File(fileName);
	}
	
	//Read the file to memory
	public void read() throws Exception{
		long fileLength = fileIn.length();
		byte[] fileData;
		try {
			BufferedInputStream inputStream = new BufferedInputStream(
					new FileInputStream(fileIn));
			DataInputStream dataInputStream = new DataInputStream(inputStream);
			// Allocate memory for reading the file into memory
			fileData = new byte[(int) fileLength];
			// Read the data to memory
			dataInputStream.read(fileData, 0, (int) fileLength);
			// Close the file after reading
			dataInputStream.close();
			readHeader(fileData);
		} catch (Exception e) {
			throw new UnsupportedDataTypeException("Could not read input file or file did not exist.");
		}
	}
	
	//Extract the header
	private void readHeader(byte[] fileData){
		int offset = 334;
		percent = Double
				.longBitsToDouble((long) (((long) (fileData[offset + 7] & 0xFF)) << 56
						| ((long) (fileData[offset + 6] & 0xFF)) << 48
						| ((long) (fileData[offset + 5] & 0xFF)) << 40
						| ((long) (fileData[offset + 4] & 0xFF)) << 32
						| ((long) (fileData[offset + 3] & 0xFF)) << 24
						| ((long) (fileData[offset + 2] & 0xFF)) << 16
						| ((long) (fileData[offset + 1] & 0xFF)) << 8 | ((long) (fileData[offset + 0] & 0xFF)) << 0));
		offset = 1316;
		site = new String(fileData, offset + 1, (int) fileData[offset]);
	}
}