package org.bonej.io;

/*
 History:
 1.9.12		Decoder for Creation Time in Scanco header implemented 
 (vms quadword, big endian, converted to unix timestamp and finally date/time).
 30.08.12	corrected errors in the Scanco header information
 updated a few comments
 added the most meaningful file header information to the ImagePlus property "Info" which can be retrieved with 
 the menu command "Show Info"
 26.08.12   removed unnecessary comments
 added missing {} in if-else-statements
 formated Scanco header information to ease readability
 14.09.09 added information about the name[] field in the header (nameStringInHeader)
 22.03.08 old and now irrelevant_code_removed
 23.03.08 KHKs_Scanco_ISQ_FileReader : Downsampling by factor 2 in x, y and z direction by calculating the average of the pixels which are reduced to one new pixel
 including adjustment of fi.width, fi.height,fi.depth and fi.unit
 16.02.08 Bugfix for files larger than 2 GB revisited
 29.06.06 pixeldistance is taken from the orginal ISQ file header and used for metric scaling in ImageJ (fi.info)
 02.08.06 Bugfix for files larger than 2 GB (adjustment of offset and longOffset so that Skip will be correct)

 This file is nothing else then copying the modules which are used by ImageJ to import
 raw data or files into this one single file.
 It helped me a lot to understand ImageJ
 and it made it easier to develop this tool

 The usage should be self-explanatory.

 In short you select a ROI for the import by entering the upper left and the lower right coordinate
 of the rectangle.

 In order to obtain these coordinates you can decide how many slices you want to import and with which slice
 you start the import.

 Then you take paper and pencil and note the required coordinates for the second import attempt.

 Note: this tool was written because we have only limited memory on our Windows PCs. It was my very first Java program!

 Meaning of the checkboxes:

 "Downsample by factor 2 in x,y,z (method = average:)"
 for a view overview on a computer with limited RAM the stack can be downscaled by the factor of 2 in x and y and z direction
 The downsampling is made by averaging 8 pixels into 1 new pixel

 Please consider: this tool was written by a dentist
 -> and modified by a veterinary surgeon ;-)

 /* Scanco ISQ Header Information: - note: Scanco uses OpenVMS on Alpha workstations

 Little endian byte order (the least significant bit occupies the lowest memory position.

 00   char    check[16];              // CTDATA-HEADER_V1
 16   int     data_type;              
 20   int     nr_of_bytes;     
 24   int     nr_of_blocks;           
 28   int     patient_index;          //p.skip(28);
 32   int     scanner_id;				//p.skip(32);
 36   int     creation_date[2];		//P.skip(36);
 44   int     dimx_p;					//p.skip(44);
 48   int     dimy_p;
 52   int     dimz_p;
 56   int     dimx_um;				//p.skip(56);
 60   int     dimy_um;
 64   int     dimz_um;
 68   int     slice_thickness_um;		//p.skip(68);
 72   int     slice_increment_um;		//p.skip(72);
 76   int     slice_1_pos_um;
 80   int     min_data_value;
 84   int     max_data_value;
 88   int     mu_scaling;             //p.skip(88);  /* p(x,y,z)/mu_scaling = value [1/cm]
 92	int     nr_of_samples;
 96	int     nr_of_projections;
 100  int     scandist_um;
 104  int     scanner_type;
 108  int     sampletime_us;
 112  int     index_measurement;
 116  int     site;                   //coded value
 120  int     reference_line_um;
 124  int     recon_alg;              //coded value
 128  char    name[40]; 		 		//p.skip(128);
 168  int     energy;        /* V     //p.skip(168);  
 172  int     intensity;     /* uA    //p.skip(172);

 ...

 508 int     data_offset;     /* in 512-byte-blocks  //p.skip(508);

 * So the first 16 bytes are a string 'CTDATA-HEADER_V1', used to identify
 * the type of data. The 'int' are all 4-byte integers.
 * 
 * dimx_p is the dimension in pixels, dimx_um the dimensions in micrometer
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

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.io.FileInfo;
import ij.io.OpenDialog;
import ij.measure.Calibration;
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;

import ij.process.ShortProcessor;

//import java.awt.Checkbox;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.math.BigInteger;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Arrays;
import java.util.Date;

import org.doube.util.UsageReporter;

/**
 * This plugin implements the Import > Scanco ISQ command.
 * 
 * @author B. Koller, SCANCO Medical AG, April 2005
 * @author K.-H. Kunzelmann, Operative Dentistry, LMU-München, Ger, April 2006
 * @author Michael Doube, RVC, London, UK. Refactoring for BoneJ 2012 2014
 */
public class ISQReader implements PlugIn {

	private static final String MAGIC = "CTDATA-HEADER_V1";

	private long skipCount;
	private int bufferSize, byteCount, nPixels;
	private int eofErrorCount;

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
		int width = imageSize[0];
		int height = imageSize[1];
		int depth = imageSize[2];

		GenericDialog gd = new GenericDialog("Import Scanco ISQ file");
		String name = getPatientName(path);
		gd.addMessage("Patient:" + name + "\n");
		gd.addMessage("\nEnter the coordinates for the bounding rectangle\n"
				+ "to crop the microCT stack during import");

		gd.addNumericField("Upper_left_X: ", 0, 0);
		gd.addNumericField("Upper_left_Y: ", 0, 0);
		gd.addNumericField("Lower_right_X", width - 1, 0);
		gd.addNumericField("Lower_right_Y: ", height - 1, 0);
		gd.addNumericField("First_slice: ", 0, 0);
		gd.addNumericField("Number_of_slices: ", depth, 0);
		gd.addCheckbox("Downsample 2x", false);

		gd.showDialog();
		if (gd.wasCanceled())
			return;

		int startX = (int) gd.getNextNumber();
		int startY = (int) gd.getNextNumber();
		int endX = (int) gd.getNextNumber();
		int endY = (int) gd.getNextNumber();
		int startZ = (int) gd.getNextNumber();
		int nSlices = (int) gd.getNextNumber();
		final boolean downsample = gd.getNextBoolean();
		// Open the file
		try {
			ImagePlus imp = openScancoISQ(path, downsample, startX, startY,
					endX, endY, startZ, nSlices);

			String scancoHeaderdata = getHeaderData(path);
			imp.setProperty(
					"Info",
					appendToFileInfo(imp, scancoHeaderdata));
			imp.show();
			UsageReporter.reportEvent(this).send();
		} catch (IllegalArgumentException e) {
			IJ.error("ISQ Reader", e.getMessage());
			return;
		}
	}

	/**
	 * Opens a Scanco ISQ file as an ImageJ ImagePlus
	 * 
	 * @param path
	 * @param downsample
	 * @param startX
	 * @param startY
	 * @param endX
	 * @param endY
	 * @param startZ
	 * @param nSlices
	 * @return
	 */
	public ImagePlus openScancoISQ(String path, boolean downsample, int startX,
			int startY, int endX, int endY, int startZ, int nSlices) {

		int[] imageSize = getImageSize(path);
		int width = imageSize[0];
		int height = imageSize[1];
		int depth = imageSize[2];
		double[] pixelSize = getPixelSize(path);
		int offset = getOffset(path);
		if (startX < 0 || startX >= width || startY < 0 || startY >= height
				|| endX < 0 || endX >= width || endY < 0 || endY >= height
				|| startZ < 0 || startZ >= depth || nSlices < 1
				|| nSlices > depth - startZ)
			throw new IllegalArgumentException(
					"Crop parameters fall outside image bounds");

		// FileInfo
		FileInfo fi = new FileInfo();
		fi.fileName = new File(path).getName();
		fi.directory = new File(path).getParent()
				+ ((IJ.isWindows()) ? "\\" : "/");
		fi.width = width;
		fi.height = height;
		fi.offset = offset;

		if (startZ > 0) {
			final long area = width * height;
			final long sliceTimesArea = area * startZ;
			// multiplication * 2 because a "short" value is 2 bytes long
			final long sliceTimesAreaTimes2 = sliceTimesArea * 2;
			long dummy = (long) fi.offset + sliceTimesAreaTimes2;

			if (dummy <= Integer.MAX_VALUE && dummy > 0) {
				// 2 is hardcoded no. of bytesPerPixel (short)
				fi.offset += (startZ * width * height * 2);
			} else {
				fi.longOffset = (long) (fi.offset + sliceTimesAreaTimes2);
			}
		}
		if (nSlices > getImageSize(path)[2] - startZ) {
			nSlices = getImageSize(path)[2] - startZ;
		}

		if (offset <= Integer.MAX_VALUE && offset > 0) {
			fi.offset = (int) offset;
		}
		if (offset > Integer.MAX_VALUE) {
			fi.longOffset = offset;
		}
		fi.nImages = nSlices;
		fi.gapBetweenImages = 0;
		fi.intelByteOrder = true;
		fi.whiteIsZero = false;
		fi.fileType = FileInfo.GRAY16_SIGNED;
		fi.pixelWidth = pixelSize[0];
		fi.pixelHeight = pixelSize[1];
		fi.pixelDepth = pixelSize[2];
		fi.unit = "mm";

		int widthStack;
		int heightStack;

		final int widthROI = endX - startX + 1;
		final int heightROI = endY - startY + 1;

		if (downsample == true) {

			fi.pixelWidth = fi.pixelWidth * 2;
			fi.pixelHeight = fi.pixelHeight * 2;
			widthStack = widthROI / 2;
			heightStack = heightROI / 2;
		} else {
			widthStack = widthROI;
			heightStack = heightROI;
		}

		// modified to match the size of the ROI
		ImageStack stack = new ImageStack(widthStack, heightStack);
		long skip = fi.longOffset > 0 ? fi.longOffset : fi.offset;
		skip += startZ * width * height * 2;

		try {
			FileInputStream is = new FileInputStream(path);

			for (int i = 1; i <= nSlices; i++) {
				IJ.showStatus("Reading: " + i + "/" + nSlices);

				//read the whole slice into an array
				//this is really inefficient if only a small ROI is needed
				//as must read all pixels off disk then throw them away
				short[] pixels = readPixels(is, skip, width, height);

				// get pixels for ROI only
				int indexCountPixels = startY * width + startX;
				int indexCountROI = 0;

				short[] pixelsROI = new short[widthROI * heightROI];

				for (int u = 0; u < heightROI; u++) {
					System.arraycopy(pixels, indexCountPixels, pixelsROI,
							indexCountROI, widthROI);
					indexCountPixels += width;
					indexCountROI += widthROI;
				}

				if (pixels == null)
					break;

	
				if (downsample) {
					ImageProcessor ip = new ShortProcessor(widthROI, heightROI);
					ip.setPixels(pixelsROI);
					ip.setInterpolationMethod(ImageProcessor.BICUBIC);
					ImageProcessor downsizedIp = ip.resize(widthROI / 2, heightROI / 2, true);
					stack.addSlice("microCT-Import_by_KHK_w_" + widthROI
							+ "_h_" + heightROI + "_slice." + i, downsizedIp);

				} else {

					for (int index = 0; index < widthROI * heightROI; index++) {
						pixelsROI[index] = (short) (pixelsROI[index] - 32768);
						if (pixelsROI[index] < 0)
							pixelsROI[index] = 0;
					}

					stack.addSlice("microCT-Import_by_KHK_w_" + widthROI
							+ "_h_" + heightROI + "_slice." + i, pixelsROI);
				}

				skip = fi.gapBetweenImages;
				IJ.showProgress((double) i / nSlices);
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

		cal.pixelWidth = (downsample) ? pixelSize[0] * 2 : pixelSize[0];
		cal.pixelHeight = (downsample) ? pixelSize[1] * 2 : pixelSize[1];
		cal.pixelDepth = pixelSize[2];
		cal.setUnit("mm");
		cal.xOrigin = -startX;
		cal.yOrigin = -startY;
		cal.zOrigin = -startZ;
		cal.setFunction(Calibration.STRAIGHT_LINE, new double[] { 0,
				1.0 / getMuScaling(path) }, "1/cm");
		imp.setCalibration(cal);
		// set display range
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
	 * from ImageReader.java: Skips the specified number of bytes, then reads an
	 * image and returns the pixel array (byte, short, int or float). Returns
	 * null if there was an IO exception. Does not close the InputStream.
	 */
	private short[] readPixels(FileInputStream in, long skipCount, int width,
			int height) {
		this.skipCount = skipCount;
		short[] pixels = readPixels(in, width, height);
		if (eofErrorCount > 0)
			return null;
		else
			return pixels;
	}

	/**
	 * Reads the image from the InputStream and returns the pixel array (byte,
	 * short, int or float). Returns null if there was an IO exception. Does not
	 * close the InputStream.
	 * 
	 * @param in
	 * @param width
	 * @param height
	 * @return
	 */
	private short[] readPixels(FileInputStream in, int width, int height) {
		try {
			skip(in, width, height);
			return read16bitImage(in);

		} catch (IOException e) {
			IJ.log("" + e);
			return null;
		}
	}

	/**
	 * Reads a 16-bit image. Signed pixels are converted to unsigned by adding
	 * 32768.
	 * 
	 * @param in
	 * @return
	 * @throws IOException
	 */
	private short[] read16bitImage(FileInputStream in) throws IOException {
		int pixelsRead;
		byte[] buffer = new byte[bufferSize];
		short[] pixels = new short[nPixels];
		int totalRead = 0;
		int base = 0;
		int count;

		while (totalRead < byteCount) {
			if ((totalRead + bufferSize) > byteCount)
				bufferSize = byteCount - totalRead;

				count = in.read(buffer, 0, bufferSize);
				
				if (count == -1) {
					eofErrorCount++;
					Arrays.fill(pixels, Short.MAX_VALUE);
					return pixels;
				}
			totalRead += bufferSize;
			// divide by two because there are 2 bytes per pixel
			pixelsRead = bufferSize / 2;
			for (int i = base, j = 0; i < (base + pixelsRead); i++, j += 2)
				pixels[i] = (short) ((((buffer[j + 1] & 0xff) << 8) | (buffer[j] & 0xff)) + 32768);
			base += pixelsRead;
		}
		return pixels;
	}

	/**
	 * 
	 * @param in
	 * @param width
	 * @param height
	 * @throws IOException
	 */
	private void skip(FileInputStream in, int width, int height)
			throws IOException {

		// This routine is called for every slice

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

			}
		} 
		byteCount = width * height * 2; //16-bit images so 2 bytes per pixel

		nPixels = width * height;
		bufferSize = byteCount / 25;
		if (bufferSize < 8192)
			bufferSize = 8192;
		else
			bufferSize = (bufferSize / 8192) * 8192;
	}

	/**
	 * Check the magic number to determine if a file is a Scanco ISQ
	 * 
	 * @param path
	 * @return true if the file is a Scanco ISQ
	 */
	public boolean isScancoISQ(String path) {
		if (getMagic(path).equals(MAGIC))
			return true;
		else
			return false;
	}

	/**
	 * Get the magic number from a file
	 * 
	 * @param path
	 * @return the first 16 bytes as a character string
	 */
	public String getMagic(String path) {
		return readString(path, 0, 16);
	}

	public int getPatientIndex(String path) {
		return readInt(path, 28);
	}

	public int getScannerID(String path) {
		return readInt(path, 32);
	}

	public Date getCreationDate(String path) {
		int[] quadWord = new int[8];

		if (path == null)
			throw new IllegalArgumentException();
		try {
			File iFile = new File(path);
			FileInputStream p = new FileInputStream(iFile);
			p.skip(36);
			for (int index = 7; index >= 0; index--) {
				quadWord[index] = p.read();
			}
			p.close();
			return vmsQuadwordToTimestamp(quadWord);
		} catch (Exception e) {
			IJ.handleException(e);
		}
		return null;
	}

	/**
	 * Get the ISQ creation date as a string
	 * @param path
	 * @return date in yyyy-MM-dd HH:mm:ss format
	 */
	public String getCreationDateAsString(String path) {
		Date date = getCreationDate(path);
		DateFormat df = new SimpleDateFormat("yyyy-MM-dd' 'HH:mm:ss");
		String reportDate = df.format(date);
		return reportDate;
	}

	/**
	 * Get the pixel size of the image
	 * @param path
	 * @return {x, y, z} image size in pixels
	 */
	public int[] getImageSize(String path) {
		int[] sizes = { readInt(path, 44), readInt(path, 48), readInt(path, 52) };
		return sizes;
	}

	/**
	 * Get the real size of the image
	 * @param path
	 * @return {x, y, z} image size in metric units
	 */
	public double[] getRealSize(String path) {
		double[] sizes = { (double) readInt(path, 56) / 1000.0,
				(double) readInt(path, 60) / 1000.0,
				(double) readInt(path, 64) / 1000.0 };
		return sizes;
	}

	public int getSliceThickness(String path) {
		return readInt(path, 68);
	}

	public int getSliceIncrement(String path) {
		return readInt(path, 72);
	}

	public int getScanDistance(String path) {
		return readInt(path, 76);
	}

	public int getMinDataValue(String path) {
		return readInt(path, 80);
	}

	public int getMaxDataValue(String path) {
		return readInt(path, 84);
	}

	public int getMuScaling(String path) {
		return readInt(path, 88);
	}

	public int getNrSamples(String path) {
		return readInt(path, 92);
	}

	public int getNrProjections(String path) {
		return readInt(path, 96);
	}

	public int getScanDistanceUm(String path) {
		return readInt(path, 100);
	}

	public int getScannerType(String path) {
		return readInt(path, 104);
	}

	public int getSampleTimeUs(String path) {
		return readInt(path, 108);
	}

	public int getMeasurementIndex(String path) {
		return readInt(path, 112);
	}

	public int getSite(String path) {
		return readInt(path, 116);
	}

	public int getReferenceLineUm(String path) {
		return readInt(path, 120);
	}

	public int getReconstructionAlgorithm(String path) {
		return readInt(path, 124);
	}

	public String getPatientName(String path) {
		return readString(path, 128, 40);
	}

	public int getEnergy(String path) {
		return readInt(path, 168);
	}

	public int getIntensity(String path) {
		return readInt(path, 172);
	}

	// what to do with 176 int fill[83]? Skip?

	public int getOffset(String path) {
		return readInt(path, 508) * 512 + 512;
	}

	/**
	 * Get the pixel spacing in real units (mm)
	 * 
	 * @param path location of file in filesystem
	 * @return {x, y, z} pixel spacing in real units (mm) 
	 */
	public double[] getPixelSize(String path) {
		int[] nPixels = getImageSize(path);
		double[] realSize = getRealSize(path);
		double[] pixelSize = { realSize[0] / nPixels[0],
				realSize[1] / nPixels[1], realSize[2] / nPixels[2] };
		return pixelSize;
	}

	private int readInt(String path, int firstByte) {
		if (path == null)
			throw new IllegalArgumentException();
		try {
			File iFile = new File(path);
			FileInputStream p = new FileInputStream(iFile);
			p.skip(firstByte);
			int value = (p.read() + p.read() * 256 + p.read() * 65536 + p
					.read() * 256 * 65536);
			p.close();
			return value;
		} catch (IOException e) {
			IJ.handleException(e);
		}
		return -1;
	}

	private String readString(String path, int firstByte, int length) {
		if (path == null)
			throw new IllegalArgumentException();
		try {
			File iFile = new File(path);
			FileInputStream p = new FileInputStream(iFile);
			p.skip(firstByte);
			String string = "";
			for (int kh = 0; kh < length; kh++) {
				char ch = (char) p.read();
				string += ch;
			}
			p.close();
			return string;
		} catch (IOException e) {
			IJ.handleException(e);
		}
		return null;
	}

	private String getHeaderData(String path) {
		if (path == null) {
			throw new IllegalArgumentException();
		}
		try {
			File iFile = new File(path);
			FileInputStream p = new FileInputStream(iFile);

			String headerData = "Scanco Header Data\n\n"
			+ "Patient Name: " + getPatientName(path) + "\n"
			+ "Patient Index: " + getPatientIndex(path) + "\n"
			+ "Site: " + getSite(path) + "\n"
			+ "Reference Line: " + getReferenceLineUm(path) + " µm\n"
			+ "Scanner-ID: " + getScannerID(path) + "\n"
			+ "Scanner_type: " + getScannerType(path) + "\n"
			+ "Creation Date: "	+ getCreationDateAsString(path) + "\n"
			+ "Slice Thickness: " + getSliceThickness(path) + " µm\n"
			+ "Slice Increment: " + getSliceIncrement(path) + " µm\n"
			+ "Min Value: " + getMinDataValue(path) + "\n"
			+ "Max Value: " + getMaxDataValue(path) + "\n"
			+ "µ-Scaling: " + getMuScaling(path) + "\n"
			+ "Scan-Distance: " + getScanDistanceUm(path) + " µm\n"
			+ "Sampletime: " + getSampleTimeUs(path) + " µs\n"
			+ "Samples: " + getNrSamples(path) + "\n"
			+ "Projections: " + getNrProjections(path) + "\n"
			+ "Reconstruction Algorithm: " + getReconstructionAlgorithm(path) + "\n"
			+ "Measurement Index: "	+ getMeasurementIndex(path) + "\n"
			+ "Energy : " + getEnergy(path) + " V\n"
			+ "Intensity : " + getIntensity(path) + " µA";

			p.close();
			return headerData;

		} catch (IOException e) {
			IJ.handleException(e);
		}
		return null;
	}

	/**
	 * converts the VMS quadword timestamp to date/time
	 * 
	 * OpenVMS and Unix Date and Time Conversions
	 * 
	 * The creation date of the Scanco ISQ files is coded as an 8 Byte sequence
	 * which is called a quadword. The byte order of the creation date is
	 * "big endian".
	 * 
	 * In summary, the 8 Bytes encode a large number, which represents the
	 * numbers of 100 nanosecond intervals since 00:00 on November 17, 1858
	 * local time; the OpenVMS Epoch.
	 * 
	 * To convert this number to a date/time just follow the recommendations of
	 * Stephen Hoffman:
	 * 
	 * "To get from the OpenVMS quadword to the C quadword, subtract the OpenVMS
	 * quadword value containing the Unix epoch value for 1-Jan-1970:00:00
	 * (0x007c95674beb4000) from the OpenVMS quadword value, and then divide by
	 * 10000000 to get from the 100ns-unit to the seconds longword."
	 * 
	 * For Java the divisor is a bit different as the date functions of Java are
	 * based on ms (milliseconds) and not on s. Therefore we just have to divide
	 * by 10000.
	 * 
	 * A word of caution has to be added here:
	 * 
	 * It took me a few hours to figure out that the use of the recommended
	 * Calender object and its methods either was not correctly initialized by
	 * me or it is simply buggy. I could not get the correct creation date. The
	 * date was always 1 month off the original value.
	 * 
	 * When I changed to the depreciated approach to use the "Date" object
	 * everything worked fine.
	 * 
	 * Read the invaluable background information from Stephen Hoffman, Hoffman
	 * Labs. Especially the conversion utility was very helpful for debugging
	 * this routine.
	 * 
	 * @see http://labs.hoffmanlabs.com/node/735
	 * @see http://labs.hoffmanlabs.com/node/282
	 * @see http://www.mpp.mpg.de/~huber/util/main/cvdate.html
	 * @param vmsQuadWord
	 *            A VMS quadword - they are 64 bit unsigned integers. The system
	 *            base date is: November 17, 1858 00:00:00.00
	 * @return Date object containing the creation date of the ISQ file
	 */
	private Date vmsQuadwordToTimestamp(int[] vmsQuadWord) {

		String hexString = "";

		// OpenVMS VAX, OpenVMS Alpha and OpenVMS I64 (as well as all Microsoft
		// Windows implementations)
		// all support and all use the little-endian byte ordering. BUT in our
		// case - found by trial and error -
		// we have to use the big-endian byte order.

		for (int index = 0; index < 8; index++) {

			// if ...else loop: just necessary to format the byte in the string
			// correct for later use
			if (vmsQuadWord[index] > 0xf) {
				hexString += Integer.toHexString(vmsQuadWord[index]);
			} else {
				hexString += "0" + Integer.toHexString(vmsQuadWord[index]);
			}
		}

		BigInteger bi = new BigInteger(hexString, 16);
		BigInteger epochAsBigInteger = new BigInteger("007C95674BEB4000", 16);
		bi = bi.subtract(epochAsBigInteger);
		BigInteger divisor = BigInteger.valueOf(10000);
		bi = bi.divide(divisor);
		long value = bi.longValue();
		Date date = new Date();
		date.setTime(value);
		return date;
	}

	/**
	 * adds the content of a string to the ImagePlus property which is labeled
	 * "Info" only the content of "Info" is displayed with the
	 * "Show Info"-Command from the menu.
	 */
	private String appendToFileInfo(ImagePlus imp,
			String newinfo) {

		String contentOfImagePlusPropertyInfo = (String) imp
				.getProperty("Info");

		if (contentOfImagePlusPropertyInfo == null) {
			contentOfImagePlusPropertyInfo = newinfo;
		} else {
			contentOfImagePlusPropertyInfo = contentOfImagePlusPropertyInfo
					+ "\n------------------------\n" + newinfo;
		}

		return contentOfImagePlusPropertyInfo;
	}

}