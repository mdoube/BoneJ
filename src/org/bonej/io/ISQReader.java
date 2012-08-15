package org.bonej.io;

/*
 History:
 14.09.09 added information about the name[] field in the header (nameStringInHeader)
 22.03.08 old and now irrelevant_code_removed
 23.03.08 KHKs_Scanco_ISQ_FileReader : Downsampling by factor 2 in x, y and z direction by calculating the average of the pixels which are reduced to one new pixel
 including adjustment of fi.width, fi.height,fi.depth and fi.unit
 16.02.08 Bugfix for files larger than 2 GB revisited
 29.06.06 pixeldistance is taken from the orginal ISQ file header and used for metric scaling in ImageJ (fi.info)
 02.08.06 Bugfix for files larger than 2 GB (adjustment of offset and longOffset so that Skip will be correct)

 This file is nothing else then copying the modules which are used by imagej to import
 raw data or files into this one single file.
 It helped me a lot to understand imagej
 and it made it easier to devolop this tool

 The usage should be self-explanatory.

 In short you select a ROI for the import by entering the upper left and the lower right coordinate
 of the rectangle.

 In order to obtain these coordinates you can deside how many slices you want to import and with which slice
 you start the import.



 Then you take paper and pencile and note the required coordinates for the second import attempt.

 Note: this tool was writen because we have only limited memory on our Windows PCs. It was my very first Java program!

 Meaning of the checkboxes:
 "Scale for lin. attenuation coeff.(LAC):"
 the ISQ files are signed integer (short) raw data (see Scanco Website for more information).
 To get a meaningful value for the linear attenuation coefficient the short values have to be divided by 4096.
 This makes it necessary to use float (4byte) instead of 2byte images, which consumes more memory
 It should be default to use the scaled images for any grey level evaluation, if the morphology is the only important criterion 
 then scaling is unnecessary and safes memory.

 "Downsample by factor 2 in x,y,z (method = average:)"
 for a view overview on a computer with limited RAM the stack can be downscaled by the factor of 2 in x and y and z direction
 The downsampling is made by averaging 8 pixels into 1 new pixel 

 "Check = distances in metric units, uncheck = pixels"
 Option to decide whether to import the images with metric units as taken from the ISQ-file or alternatively keep the distances in pixels
 The pixel variant is useful to determine the ROI, otherwise it is better to use the metric units

 "8-bit-import (overrules 'Scale for LAC')"
 To save even more RAM then with Short one can decide to use 8-Bit stacks only. 
 The 8-Bit-Stack is a massive data reduction. The values are a simple division of the short values by 256.
 For better visability - auto-adjust brightness and contrast yourself.
 To get a quick idea of the LAC based on the 8-bit-data you can multiply by 256 and divide by 4096 yourself. But be aware you will lose precision this way due to rounding errors. 

 Please consider: this tool was writen by a dentist
 I have no programming experience.
 I am satisfied if the code is running.
 Programming esthetics is nice, but it is beyond my capabilities.


 TODOs: Select the ROI with the mouse
 */

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.measure.Calibration;
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;
import ij.gui.GenericDialog;
import ij.io.FileInfo;
import ij.io.OpenDialog;

import org.doube.util.UsageReporter;

/**
 * This plugin implements the Import > Scanco ISQ command.
 * 
 * @author B. Koller, SCANCO Medical AG, April 2005
 * @author K.-H. Kunzelmann, Operative Dentistry, LMU-München, Ger, April 2006
 */
public class ISQReader implements PlugIn {

	private static final String MAGIC = "CTDATA-HEADER_V1";

	private FileInfo fi;
	private long skipCount;
	private int bytesPerPixel, bufferSize, byteCount, nPixels;
	private int eofErrorCount;

	// Anpassung für Files > 2 GB
	// wäre aber vermutlich gar nicht nötig. änderung bei Zeile 276 (ca) haette
	// vermutlich gereicht

	private boolean downsample = false;
	private boolean eightBitOnly = false;

	// necessary for the clip ROI

	private int upperLeftX, upperLeftY, lowerRightX, lowerRightY;
	private int startROI, gapBetweenLines, heightROI, widthROI, nFirstSlice;

	public void run(String arg) {

		// the ISQ-File is selected

		OpenDialog od = new OpenDialog("Open ISQ...", arg);
		String directory = od.getDirectory();
		String fileName = od.getFileName();
		String path = directory + fileName;
		if (fileName == null)
			return;
		if (!isScancoISQ(path)) {
			IJ.error("ISQ Reader", "Not an ISQ file. Magic number ("
					+ getMagic(path) + ")" + " does not match.");
			return;
		}

		int[] imageSize = getImageSize(path);
		double[] pixelSize = getPixelSize(path);
		int offset = getOffset(path);
		int xdimension = imageSize[0];
		int ydimension = imageSize[1];
		int zdimension = imageSize[2];

		getMagic(path);

		// FileInfo
		fi = new FileInfo();

		fi.fileName = fileName;
		fi.directory = directory;
		fi.width = xdimension;
		fi.height = ydimension;

		// hier Anpassung fuer Files > 2 GB

		if (offset <= Integer.MAX_VALUE && offset > 0) {
			fi.offset = (int) offset;
		}
		if (offset > Integer.MAX_VALUE) {
			fi.longOffset = offset;
		}
		fi.nImages = zdimension;
		fi.gapBetweenImages = 0;
		fi.intelByteOrder = true;
		fi.whiteIsZero = false;
		fi.fileType = FileInfo.GRAY16_SIGNED;
		fi.pixelWidth = pixelSize[0];
		fi.pixelHeight = pixelSize[1];
		fi.pixelDepth = pixelSize[2];
		fi.unit = "mm";

		// Generic dialog to input the ROI-coordinates
		try {
			getRoiCoordinates(path);
		} catch (Exception e) {
			IJ.error("ISQ Reader", e.getMessage());
			return;
		}

		// Open the file
		ImagePlus imp = openScancoISQ(path);
		imp.show();
		UsageReporter.reportEvent(this).send();
	}

	/**
	 * Generic dialog to input the ROI-coordinates
	 * 
	 * @param path
	 *            Path to ISQ file (directory + filename)
	 */
	private void getRoiCoordinates(String path) {

		int[] imageSize = getImageSize(path);
		int width = imageSize[0];
		int height = imageSize[1];
		int depth = imageSize[2];

		GenericDialog gd = new GenericDialog(
				"Kunzelmann: Import Scanco ISQ-Files");

		nFirstSlice = 0;
		String name = getName(path);
		gd.addMessage(name + "\n");

		gd.addMessage("\nEnter the coordinates for the bounding rectangle\nto crop the microCT stack during import");

		gd.addNumericField("Upper_left_X: ", upperLeftX, 0);
		gd.addNumericField("Upper_left_Y: ", upperLeftY, 0);
		gd.addNumericField("Lower_right_X", width - 1, 0);
		gd.addNumericField("Lower_right_Y: ", height - 1, 0);
		gd.addNumericField("First_slice: ", nFirstSlice, 0);
		gd.addNumericField("Number_of_slices: ", depth, 0);
		gd.addCheckbox("Downsample 2x", downsample);

		gd.showDialog();
		if (gd.wasCanceled()) {
			fi.nImages = 0;
			nFirstSlice = 1;
			return;
		}

		upperLeftX = (int) gd.getNextNumber();
		upperLeftY = (int) gd.getNextNumber();
		lowerRightX = (int) gd.getNextNumber();
		lowerRightY = (int) gd.getNextNumber();
		nFirstSlice = (int) gd.getNextNumber();
		fi.nImages = (int) gd.getNextNumber();
		downsample = gd.getNextBoolean();

		if (upperLeftX < 0 || upperLeftX >= width || upperLeftY < 0
				|| upperLeftY >= height || lowerRightX < 0
				|| lowerRightX >= width || lowerRightY < 0
				|| lowerRightY >= height || nFirstSlice < 0
				|| nFirstSlice >= depth || fi.nImages < 1
				|| fi.nImages > depth - nFirstSlice) {
			throw new IllegalArgumentException(
					"Crop parameters fall outside image bounds");
		}

		startROI = upperLeftY * width + upperLeftX;

		gapBetweenLines = (width - lowerRightX) + upperLeftX - 1;
		widthROI = lowerRightX - upperLeftX + 1;
		heightROI = lowerRightY - upperLeftY + 1;

		// ***********************************************
		// Anpassung wegen Files > 2 GB
		// ich habe aus offset ein long statt integer gemacht
		// ich habe Abfrage ergänzt of fi.offset oder fi.longOffset verwendet
		// werden soll
		// ich habe hier änderung ergänzt mit Abfrage zur Grösse

		if (nFirstSlice > 0) {
			long dummy;

			long area;
			long sliceTimesArea;
			long sliceTimesAreaTimes2;

			area = fi.width * fi.height;
			sliceTimesArea = area * nFirstSlice;
			// * 2 wegen Short = 2 Byte
			sliceTimesAreaTimes2 = sliceTimesArea * 2;

			dummy = (long) fi.offset + sliceTimesAreaTimes2;

			if (dummy <= Integer.MAX_VALUE && dummy > 0) {
				// 2 is hardcoded no. of bytesPerPixel (short)
				fi.offset = fi.offset
						+ (nFirstSlice * fi.width * fi.height * 2);
			} else {
				fi.longOffset = (long) (fi.offset + sliceTimesAreaTimes2);
			}
		}
		System.out.println("offset: " + fi.offset);
		System.out.println("longoffset: " + fi.longOffset);
		if (fi.nImages > getImageSize(path)[2] - nFirstSlice) {
			fi.nImages = getImageSize(path)[2] - nFirstSlice;
		}
	}

	/** Opens a stack of images. */
	public ImagePlus openScancoISQ(String path) {
		System.out.println("downsample: " + downsample);

		int widthStack = 0;
		int heightStack = 0;

		if (downsample == true) {

			fi.pixelWidth = fi.pixelWidth * 2;
			fi.pixelHeight = fi.pixelHeight * 2;
			fi.pixelDepth = getPixelSize(path)[2] * 2;
			widthStack = widthROI / 2;
			heightStack = heightROI / 2;
		} else {
			widthStack = widthROI;
			heightStack = heightROI;
		}

		// temp stack, needed later for downsampling
		float[] downsampledPixels32_temp = new float[(widthROI * heightROI)
				/ (2 * 2)];

		// modified to match the size of the ROI
		ImageStack stack = new ImageStack(widthStack, heightStack);
		long skip = fi.longOffset > 0 ? fi.longOffset : fi.offset;
		// System.out.println("longoffset-vor-skip: "+fi.longOffset);
		// System.out.println("offset-vor-skip: "+offset);
		// System.out.println("skip: "+skip);
		System.out.println("after Imagestack -> x: " + fi.pixelWidth + " ; y: "
				+ fi.pixelHeight);

		try {
			FileInputStream is = new FileInputStream(path);

			System.out.println("Try FileInputStream");

			// Obsolet comment: I reduce the no of slices by 1 to
			// avoid a nullpointerexception error
			for (int i = 1; i <= fi.nImages; i++) {
				IJ.showStatus("Reading: " + i + "/" + fi.nImages);
				// System.out.println("fi.nImages: "+fi.nImages);
				short[] pixels = readPixels(is, skip);

				// get pixels for ROI only
				int indexCountPixels = startROI;
				int indexCountROI = 0;

				short[] pixelsROI;
				pixelsROI = new short[widthROI * heightROI];

				for (int u = 0; u < heightROI; u++) {
					System.arraycopy(pixels, indexCountPixels, pixelsROI,
							indexCountROI, widthROI);
					indexCountPixels = indexCountPixels + widthROI
							+ gapBetweenLines;
					indexCountROI = indexCountROI + widthROI;
					// System.out.println(i+"::"+"indexCountPixels:"+indexCountPixels+":"+"IndexCountROI:"+indexCountROI+":"+"Size:"+widthROI+heightROI);
				}

				if (pixels == null)
					break;

				float[] pixels32 = new float[widthROI * heightROI];
				for (int s = 0; s < widthROI * heightROI; s++) {
					pixels32[s] = (pixelsROI[s] & 0xffff);
					pixels32[s] = pixels32[s] - 32768;
					// hier wird nun ISQ bzgl. des lin. att. Coeff. skaliert
					// (nächste Zeile)
					pixels32[s] = pixels32[s] / 4096;
					if (pixels32[s] < 0) {
						pixels32[s] = 0;
					}
				}

				if (downsample == true) {
					// System.out.println("Downsample loop ... ");
					float[] downsampledPixels32 = new float[(widthROI * heightROI)
							/ (2 * 2)];
					// float[] downsampledPixels32_temp = new
					// float[(widthROI*heightROI)/(2*2)];
					short[] downsampledPixels_av = new short[(widthROI * heightROI)
							/ (2 * 2)];

					int index = 0;
					// here we calculate the average in the x,y plane.
					for (int h = 0; h < heightROI - 1; h = h + 2) {
						for (int w = 0; w < widthROI - 1; w = w + 2) {
							downsampledPixels32[index] = ((pixels32[(h * widthROI)
									+ w]
									+ pixels32[(h * widthROI) + w + 1]
									+ pixels32[((h + 1) * widthROI) + w] + pixels32[((h + 1) * widthROI)
									+ w + 1]) / 4);
							index = index + 1;
							if (index >= widthStack * heightStack)
								index = 0;
							// System.out.print(".");
						}
					}
					if (i % 2 > 0) {
						System.arraycopy(downsampledPixels32, 0,
								downsampledPixels32_temp, 0,
								downsampledPixels32.length);
					} else {
						float temp1, temp2, temp3 = 0;
						for (int s = 0; s < heightStack * widthStack; s++) {
							temp1 = downsampledPixels32[s];
							temp2 = downsampledPixels32_temp[s];
							temp3 = ((temp1 + temp2) / 2) * 4096;
							if (temp3 < 0.0)
								temp3 = 0.0f;
							if (temp3 > 65535.0)
								temp3 = 65535.0f;
							downsampledPixels_av[s] = (short) temp3;
						}

						// *** shortTo8bit

						int value;
						byte[] pixels8 = new byte[heightStack * widthStack];
						for (int index2 = 0; index2 < heightStack * widthStack; index2++) {
							value = downsampledPixels_av[index2] / 256;
							if (value > 255)
								value = 255;
							pixels8[index2] = (byte) value;
						}

						if (eightBitOnly == true) {
							stack.addSlice("microCT-Import_by_KH_w_"
									+ widthStack + "_h_" + heightStack
									+ "_slice." + i, pixels8);
						} else {
							stack.addSlice("microCT-Import_by_KH_w_"
									+ widthStack + "_h_" + heightStack
									+ "_slice." + i, downsampledPixels_av);
						}

						// *** shortTo8bit ende

						// stack.addSlice("microCT-Import_by_KH_w_"+widthROI/2+"_h_"+heightROI/2+"_slice."+
						// i, downsampledPixels_av);
					}
				} else {
					// System.out.println("Instead of downsample loop ... ELSE loop");
					for (int index = 0; index < widthROI * heightROI; index++) {
						pixelsROI[index] = (short) (pixelsROI[index] - 32768);
						if (pixelsROI[index] < 0)
							pixelsROI[index] = 0;
					}
					// *** shortTo8bit

					int value;
					byte[] pixels8 = new byte[widthROI * heightROI];
					for (int index2 = 0; index2 < widthROI * heightROI; index2++) {
						value = pixelsROI[index2] / 256;
						if (value > 255)
							value = 255;
						pixels8[index2] = (byte) value;
					}

					// *** shortTo8bit ende
					// System.out.println("Instead of downsample loop ... ELSE loop - just before add stack");

					if (eightBitOnly == true) {
						stack.addSlice("microCT-Import_by_KH_w_" + widthROI
								+ "_h_" + heightROI + "_slice." + i, pixels8);
					} else {
						stack.addSlice("microCT-Import_by_KH_w_" + widthROI
								+ "_h_" + heightROI + "_slice." + i, pixelsROI);
					}
				}

				// *************************************************
				// war orginal
				// stack.addSlice("microCT-Import_by_KH_w_"+widthROI+"_h_"+heightROI+"_slice."+
				// i, pixelsROI);
				// }

				skip = fi.gapBetweenImages;
				IJ.showProgress((double) i / fi.nImages);
			}
			is.close();
		} catch (Exception e) {
			IJ.log("" + e);
		} catch (OutOfMemoryError e) {
			IJ.outOfMemory(fi.fileName);
			stack.trim();
		}
		if (stack.getSize() == 0)
			return null;
		if (fi.sliceLabels != null && fi.sliceLabels.length <= stack.getSize()) {
			for (int i = 0; i < fi.sliceLabels.length; i++)
				stack.setSliceLabel(fi.sliceLabels[i], i + 1);
		}
		ImagePlus imp = new ImagePlus(fi.fileName, stack);
		Calibration cal = imp.getCalibration();

		if (fi.info != null)
			imp.setProperty("Info", fi.info);
		imp.setFileInfo(fi);
		System.out.println("after imp.show() -> x: " + fi.pixelWidth + " ; y: "
				+ fi.pixelHeight);

		double[] pixelSize = getPixelSize(path);

		cal.pixelWidth = (downsample) ? pixelSize[0] * 2 : pixelSize[0];
		cal.pixelHeight = (downsample) ? pixelSize[1] * 2 : pixelSize[1];
		cal.pixelDepth = (downsample) ? pixelSize[2] * 2 : pixelSize[2];
		cal.setUnit("mm");
		cal.setFunction(Calibration.STRAIGHT_LINE, new double[] { 0,
				1.0 / getMuScaling(path) }, "1/cm");
		imp.setCalibration(cal);
		//set display range
		double min = Double.MAX_VALUE;
		double max = -Double.MAX_VALUE;
		final int n = stack.getSize();
		for (int i = 1; i <= n; i++) {
			IJ.showStatus("Calculating stack min and max: " + i + "/" + n);
			ImageProcessor ip = stack.getProcessor(i);
			max = Math.max(max, ip.getMax());
			min = Math.min(min, ip.getMin());
		}
		imp.getProcessor().setMinAndMax(min, max);
		IJ.showProgress(1.0);
		return imp;
	}

	/** *********************************************************************** **/
	/**
	 * aus ImageReader.java: Skips the specified number of bytes, then reads an
	 * image and returns the pixel array (byte, short, int or float). Returns
	 * null if there was an IO exception. Does not close the InputStream.
	 */
	private short[] readPixels(FileInputStream in, long skipCount) {
		this.skipCount = skipCount;
		short[] pixels = readPixels(in);
		if (eofErrorCount > 0)
			return null;
		else
			return pixels;
	}

	/**
	 * Reads the image from the InputStream and returns the pixel array (byte,
	 * short, int or float). Returns null if there was an IO exception. Does not
	 * close the InputStream.
	 */
	private short[] readPixels(FileInputStream in) {
		try {

			bytesPerPixel = 2;
			skip(in);
			return read16bitImage(in);

		} catch (IOException e) {
			IJ.log("" + e);
			return null;
		}
	}

	/**
	 * *********************** this is the central import routine
	 * ***********************************
	 **/
	// there is still room to reduce code bits which are not really necessary
	// for the ISQ-Import.

	/**
	 * Reads a 16-bit image. Signed pixels are converted to unsigned by adding
	 * 32768.
	 */
	private short[] read16bitImage(FileInputStream in) throws IOException {
		int pixelsRead;
		byte[] buffer = new byte[bufferSize];
		short[] pixels = new short[nPixels];
		int totalRead = 0;
		int base = 0;
		int count;
		int bufferCount;

		while (totalRead < byteCount) {
			if ((totalRead + bufferSize) > byteCount)
				bufferSize = byteCount - totalRead;
			bufferCount = 0;

			while (bufferCount < bufferSize) { // fill the buffer
				count = in.read(buffer, bufferCount, bufferSize - bufferCount);
				if (count == -1) {
					eofError();
					if (fi.fileType == FileInfo.GRAY16_SIGNED)
						for (int i = base; i < pixels.length; i++)
							pixels[i] = (short) 32768;
					return pixels;
				}
				bufferCount += count;
			}
			totalRead += bufferSize;
			pixelsRead = bufferSize / bytesPerPixel;
			if (fi.intelByteOrder) {
				for (int i = base, j = 0; i < (base + pixelsRead); i++, j += 2)
					pixels[i] = (short) ((((buffer[j + 1] & 0xff) << 8) | (buffer[j] & 0xff)) + 32768);

			} else {
				for (int i = base, j = 0; i < (base + pixelsRead); i++, j += 2)
					pixels[i] = (short) ((((buffer[j] & 0xff) << 8) | (buffer[j + 1] & 0xff)) + 32768);
			}
			base += pixelsRead;
		}
		return pixels;
	}

	private void skip(FileInputStream in) throws IOException {

		// I count, how often this routine is used:
		// System.out.println("skip_kh called");
		// answer: called for every slice

		if (skipCount > 0) {
			long bytesRead = 0;
			int skipAttempts = 0;
			long count;
			while (bytesRead < skipCount) {
				count = in.skip(skipCount - bytesRead);
				skipAttempts++;
				if (count == -1 || skipAttempts > 5)
					break;
				bytesRead += count;
				// IJ.log("skip: "+skipCount+" "+count+" "+bytesRead+" "+skipAttempts);
			}
		}
		byteCount = fi.width * fi.height * bytesPerPixel;

		nPixels = fi.width * fi.height;
		bufferSize = byteCount / 25;
		if (bufferSize < 8192)
			bufferSize = 8192;
		else
			bufferSize = (bufferSize / 8192) * 8192;
	}

	private void eofError() {
		eofErrorCount++;
	}

	// **----------------------------------------------------------------*/
	/*
	 * Scanco ISQ Header Information:
	 * 
	 * typedef struct { /*--------------------------------------------- 00 char
	 * check[16]; // Char is in Java 2 Byte 16 int data_type; // Int = 4 Byte 20
	 * int nr_of_bytes; /* either one of them 24 int nr_of_blocks; /* or both,
	 * but min. of 1 28 int patient_index; /* 1 block = 512 bytes 32 int
	 * scanner_id; 36 int creation_date[2];
	 * /*--------------------------------------------- 40 int dimx_p; 44 int
	 * dimy_p; 48 int dimz_p; 52 int dimx_um; 56 int dimy_um; 60 int dimz_um; 64
	 * int slice_thickness_um; 68 int slice_increment_um; 72 int slice_1_pos_um;
	 * 76 int min_data_value; 78 int max_data_value; 82 int mu_scaling; /*
	 * p(x,y,z)/mu_scaling = value [1/cm] 86 int nr_of_samples; 90 int
	 * nr_of_projections; 94 int scandist_um; 98 int scanner_type; 102 int
	 * sampletime_us; 106 int index_measurement; 110 int site; /* Coded value
	 * 114 int reference_line_um; 120 int recon_alg; /* Coded value 124 char
	 * name[40]; // KHK char evtl. nur 1 Byte... siehe unten ?????? int energy;
	 * /*V int intensity; /* uA int fill[83];
	 * /*--------------------------------------------- int data_offset; /* in
	 * 512-byte-blocks } ima_data_type, *ima_data_typeP;
	 * 
	 * So the first 16 bytes are a string 'CTDATA-HEADER_V1', used to identify
	 * the type of data. The 'int' are all 4-byte integers.
	 * 
	 * dimx_p is the dimension in pixels, dimx_um the dimension in microns.
	 * 
	 * So dimx_p is at byte-offset 40, then dimy_p at 44, dimz_p (=number of
	 * slices) at 48.
	 * 
	 * The microCT calculates so called 'x-ray linear attenuation' values. These
	 * (float) values are scaled with 'mu_scaling' (see header, e.g. 4096) to
	 * get to the signed 2-byte integers values that we save in the .isq file.
	 * 
	 * e.g. Pixel value 8192 corresponds to lin. att. coeff. of 2.0 [1/cm]
	 * (8192/4096)
	 * 
	 * Following to the headers is the data part. It is in 2-byte short integers
	 * (signed) and starts from the top-left pixel of slice 1 to the left, then
	 * the next line follows, until the last pixel of the last sclice in the
	 * lower right.
	 */

	public boolean isScancoISQ(String path) {
		if (getMagic(path).equals(MAGIC))
			return true;
		else
			return false;
	}

	public String getMagic(String path) {
		if (path == null)
			throw new IllegalArgumentException();
		try {
			File iFile = new File(path);
			FileInputStream p = new FileInputStream(iFile);
			String magic = "";
			for (int kh = 0; kh < 16; kh++) {
				char ch = (char) p.read();
				magic += ch;
			}
			p.close();
			return magic;
		} catch (IOException e) {
			IJ.handleException(e);
		}
		return null;
	}

	public int[] getImageSize(String path) {
		if (path == null)
			throw new IllegalArgumentException();
		try {
			File iFile = new File(path);
			FileInputStream p = new FileInputStream(iFile);
			p.skip(44);
			int width = p.read() + p.read() * 256 + p.read() * 65536;
			p.skip(1);
			int height = p.read() + p.read() * 256 + p.read() * 65536;
			p.skip(1);
			int depth = p.read() + p.read() * 256 + p.read() * 65536;
			int[] dimensions = { width, height, depth };
			p.close();
			return dimensions;
		} catch (IOException e) {
			IJ.handleException(e);
		}
		return null;
	}

	public double[] getRealSize(String path) {
		if (path == null)
			throw new IllegalArgumentException();
		try {
			File iFile = new File(path);
			FileInputStream p = new FileInputStream(iFile);
			p.skip(56);
			double width = (p.read() + p.read() * 256 + p.read() * 65536 + p
					.read() * 256 * 65536);
			double height = (p.read() + p.read() * 256 + p.read() * 65536 + p
					.read() * 256 * 65536);
			double depth = (p.read() + p.read() * 256 + p.read() * 65536 + p
					.read() * 256 * 65536);
			width /= 1000;
			height /= 1000;
			depth /= 1000;
			double[] dimensions = { width, height, depth };
			p.close();
			return dimensions;
		} catch (IOException e) {
			IJ.handleException(e);
		}
		return null;
	}

	public double[] getPixelSize(String path) {
		int[] nPixels = getImageSize(path);
		double[] realSize = getRealSize(path);
		double[] pixelSize = { realSize[0] / nPixels[0],
				realSize[1] / nPixels[1], realSize[2] / nPixels[2] };
		return pixelSize;
	}

	public int getMuScaling(String path) {
		if (path == null)
			throw new IllegalArgumentException();
		try {
			File iFile = new File(path);
			FileInputStream p = new FileInputStream(iFile);
			p.skip(88);
			int muScaling = (p.read() + p.read() * 256 + p.read() * 65536 + p
					.read() * 256 * 65536);
			p.close();
			return muScaling;
		} catch (IOException e) {
			IJ.handleException(e);
		}
		return -1;
	}

	public String getName(String path) {
		if (path == null)
			throw new IllegalArgumentException();
		try {
			File iFile = new File(path);
			FileInputStream p = new FileInputStream(iFile);
			p.skip(128);
			String name = "";
			for (int kh = 0; kh < 40; kh++) {
				char ch = (char) p.read();
				name += ch;
				// System.out.println(nameStringInHeader);
			}
			p.close();
			return name;
		} catch (IOException e) {
			IJ.handleException(e);
		}
		return null;
	}

	public int getOffset(String path) {
		if (path == null)
			throw new IllegalArgumentException();
		try {
			File iFile = new File(path);
			FileInputStream p = new FileInputStream(iFile);
			p.skip(508);
			int offset = (p.read() + p.read() * 256 + p.read() * 65536 + 1) * 512;
			p.close();
			return offset;
		} catch (IOException e) {
			IJ.handleException(e);
		}
		return -1;
	}

}
