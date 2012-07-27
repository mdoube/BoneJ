package org.bonej.io;

//This file is a conversion utility to downsize the ISQ files. It keeps the ISQ header and crops the area outside the
//user defined ROI.

//I did not look at the code for a very long time.
//

/**
 *
 * @author Prof. Dr. Karl-Heinz Kunzelmann
 * @author karl-heinz@kunzelmann.de
 */

import ij.plugin.PlugIn;
import ij.*;
import ij.io.*;
import ij.gui.*;
import java.io.*;

//import ij.process.*;

/** This plugin implements the Acquire/ISQ command. */
public class ISQCropper implements PlugIn {

	private int width, height;
	private long skipCount;
	private int bytesPerPixel, bufferSize, byteCount, nPixels;

	private int eofErrorCount;
	private String path;

	int i, offset, xdimension, ydimension, zdimension;
	long skipSlices = 0;
	int tmpInt;

	float el_size_mm_x, el_size_mm_y, el_size_mm_z;

	// necessary for the clip ROI

	private int upperLeftX, upperLeftY, lowerRightX, lowerRightY;
	int startROI, endROI, gapBetweenLines, heightROI, widthROI, nFirstSlice;

	short[] pixels;

	private static String defaultDirectory = null;

	private int bufferCounter = 0; // buffercounter := number of buffers = no.
									// of zdimension plus 1 x offset;

	public void run(String arg) {

		// the ISQ-File is selected

		OpenDialog od = new OpenDialog("Open ISQ...", arg);
		String directory = od.getDirectory();
		String fileName = od.getFileName();
		path = directory + fileName;
		if (fileName == null)
			return;

		// a FileInputStream is opened with the file, which was selected above

		try {
			File iFile = new File(directory + fileName);
			FileInputStream p = new FileInputStream(iFile);

			// reading the ASCII File-header

			p.skip(44);
			xdimension = p.read() + p.read() * 256 + p.read() * 65536;
			p.skip(1);
			ydimension = p.read() + p.read() * 256 + p.read() * 65536;
			p.skip(1);
			zdimension = p.read() + p.read() * 256 + p.read() * 65536;
			p.skip(1);
			tmpInt = (p.read() + p.read() * 256 + p.read() * 65536 + p.read() * 256 * 65536);
			el_size_mm_x = tmpInt / xdimension;
			tmpInt = (p.read() + p.read() * 256 + p.read() * 65536 + p.read() * 256 * 65536);
			el_size_mm_y = tmpInt / ydimension;
			tmpInt = (p.read() + p.read() * 256 + p.read() * 65536 + p.read() * 256 * 65536);
			el_size_mm_z = tmpInt / zdimension;
			el_size_mm_x = el_size_mm_x / 1000;
			el_size_mm_y = el_size_mm_y / 1000;
			el_size_mm_z = el_size_mm_z / 1000;
			width = xdimension;
			height = ydimension;

			p.skip(440);

			offset = (p.read() + p.read() * 256 + p.read() * 65536 + 1) * 512;

			p.close();

		} catch (IOException e) {

		}

		// Generic dialog to input the ROI-coordinates
		getRoiCoordinates();

		// Open the file

		IJ.showStatus("Cropping started.");

		try {
			copyISQ(path);
		} catch (IOException e) {
			IJ.showStatus("Cropping made a mistake somewhere.");
		}
		IJ.showStatus("Cropping finished.");

	}

	/**
	 * *************************************end of "run"
	 * ******************************************
	 **/
	/**
	 * *************************************************************************
	 * *******************
	 **/

	// Generic dialog to input the ROI-coordinates
	void getRoiCoordinates() {

		GenericDialog gd = new GenericDialog(
				"Coordinates for the input selection");

		String text1 = new String("offset: " + offset + "\nxdimension: "
				+ xdimension + "\nydimension: " + ydimension + "\nzdimension: "
				+ zdimension);
		String text2 = new String("el_size x (in mm): " + el_size_mm_x
				+ "\nel_size y (in mm): " + el_size_mm_y
				+ "\nel_size z (in mm): " + el_size_mm_z);
		nFirstSlice = 0; // Frage: könnte die Nullpointer exeption hier
							// begründet sein - normal von 0 bis n-1, hier ab
							// 1???

		gd.addMessage(text1);
		gd.addMessage(text2);
		gd.addMessage("\nEnter the coordinates for the \nbounding rectangle to crop the \nmicroCT stack during conversion");

		gd.addNumericField("x-coord. of the upper left corner: ", upperLeftX, 0);
		gd.addNumericField("y-coord. of the upper left corner: ", upperLeftY, 0);
		gd.addNumericField("x-coord. of the lower right corner: ", width - 1, 0);
		gd.addNumericField("y-coord. of the lower right corner: ", height - 1,
				0);
		gd.addNumericField("no. of  z-slices: ", zdimension, 0);
		gd.addNumericField("start import at slice no.: ", nFirstSlice, 0);

		gd.showDialog();
		if (gd.wasCanceled())
			return;

		upperLeftX = (int) gd.getNextNumber();
		upperLeftY = (int) gd.getNextNumber();
		lowerRightX = (int) gd.getNextNumber();
		lowerRightY = (int) gd.getNextNumber();
		zdimension = (int) gd.getNextNumber();
		nFirstSlice = (int) gd.getNextNumber();

		startROI = upperLeftY * width + upperLeftX;

		gapBetweenLines = (width - lowerRightX) + upperLeftX - 1;
		widthROI = lowerRightX - upperLeftX + 1;
		heightROI = lowerRightY - upperLeftY + 1;

		if (nFirstSlice > 0) {

			// kurze Bemerkung nebenbei: diese Zeile und die Variable skipSlices
			// scheint mir eigentlich ueberfluessig !!

			skipSlices = (offset + nFirstSlice * width * height * 2); // 2 is
																		// hardcoded
																		// no.
																		// of
																		// bytesPerPixel
																		// (short)
		}

		if (zdimension > zdimension - nFirstSlice) {
			zdimension = zdimension - nFirstSlice; // Frage: könnte die
													// Nullpointer exeption hier
													// begründet sein - normal
													// von 0 bis n-1, hier ab
													// 1???
		}

	}

	/**
	 * The static method that actually performs the file copy. Before copying
	 * the file, however, it performs a lot of tests to make sure everything is
	 * as it should be.
	 */
	public void copyISQ(String in) throws IOException {
		String to_file = path + "_cropped";
		FileInputStream from = null;
		FileOutputStream to = null;

		try {
			from = new FileInputStream(in); // Create input stream
			to = new FileOutputStream(to_file); // Create output stream

			// 1. Section: read the header and write back a modified version

			/**
			 * Strategie: anhand von offset ersten Puffer so gross wie offset
			 * wählen dann diesen Header patchen mit crop daten Diesen Header
			 * raus schreiben
			 **/

			byte[] headerBuffer = new byte[offset]; // A buffer to hold the
													// header data
			byte[] sliceBuffer = new byte[width * height * 2]; // A buffer to
																// hold the
																// contents of
																// one slice
			byte[] sliceBufferROI = new byte[widthROI * heightROI * 2];

			// patch headerBuffer with the new settings

			while (bufferCounter < 1) {

				from.read(headerBuffer);

				// There are two kinds of division in Java, integer and floating
				// point. They both use the / operator to ensure you will get
				// them confused. It depends on whether the operands surrounding
				// it are int/long or float/double which form is used. Integer
				// division always gives an integer result, no fraction.
				// KH man könnte hieraus leicht eine Methode machen
				// was ich auch noch nicht kontrolliert habe, ist die Java
				// Besonderheit mit den signed bytes
				// beim Debuggen auch noch nicht ausprobiert wenn Zeilenlänge
				// 2048 etc.

				int temp1 = widthROI / 65536;
				int temp2 = (widthROI - temp1) / 256;
				int temp3 = (widthROI - temp1) % 256;
				headerBuffer[44] = (byte) temp3; // p.skip(44);
				headerBuffer[45] = (byte) temp2;
				headerBuffer[46] = (byte) temp1;

				int temp4 = heightROI / 65536;
				int temp5 = (heightROI - temp4) / 256;
				int temp6 = (heightROI - temp4) % 256;
				headerBuffer[48] = (byte) temp6;
				headerBuffer[49] = (byte) temp5;
				headerBuffer[50] = (byte) temp4;

				int temp7 = zdimension / 65536;
				int temp8 = (zdimension - temp7) / 256;
				int temp9 = (zdimension - temp7) % 256;
				headerBuffer[52] = (byte) temp9;
				headerBuffer[53] = (byte) temp8;
				headerBuffer[54] = (byte) temp7;

				to.write(headerBuffer);
				bufferCounter = bufferCounter + 1;
			}

			/**
			 * Strategie:
			 * 
			 * dann den zweiten buffer definieren und bis zum schluss dann wie
			 * gehabt abarbeiten.
			 **/
			int count;

			// Section 2:
			// I need this look to skip the amount of bytes which should not be
			// imported

			for (int z = 0; z <= nFirstSlice; z++) {
				count = from.read(sliceBuffer);
				if (count == -1) {
					eofError();
				}
			}

			// Section 3: I import the ROI data and write them back as an ISQ
			// File.

			// zdimension is the no of slices to be imported

			while (bufferCounter <= zdimension) {

				count = from.read(sliceBuffer);

				if (count == -1) {
					eofError();
				}

				// I copy only the content of the ROI into the buffer which has
				// to be writen to the Output-Stream

				int temp = 0;

				for (int u = 0; u < heightROI; u++) { // Lines in the ROI buffer
					for (int v = 0, j = 0; v < widthROI; v++, j = j + 2) { // columns
																			// in
																			// the
																			// ROI
																			// buffer

						// do not forget: we have short = 2 byte data

						sliceBufferROI[u * widthROI * 2 + j] = sliceBuffer[startROI
								* 2 + u * (widthROI + gapBetweenLines) * 2 + j];
						sliceBufferROI[u * widthROI * 2 + j + 1] = sliceBuffer[startROI
								* 2
								+ u
								* (widthROI + gapBetweenLines)
								* 2
								+ j
								+ 1];

					}
				}

				IJ.showProgress((double) bufferCounter / zdimension);
				to.write(sliceBufferROI);

				bufferCounter++;
			}

		}

		catch (IOException e) {
			IJ.showStatus("Cropping made a mistake in write loop.");
		}

		// Always close the streams, even if exceptions were thrown
		finally {
			if (from != null)
				try {
					from.close();
				} catch (IOException e) {
					;
				}
			if (to != null)
				try {
					to.close();
				} catch (IOException e) {
					;
				}
		}
	}

	void eofError() {
		eofErrorCount++;
	}

}