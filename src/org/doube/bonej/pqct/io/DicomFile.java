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

package ImageReading;
/*For saving several slices from a single file*/
public class DicomFile{
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
	public int max;
	public int min;
	public DicomFile(){
		//Do Nothing
	}
}