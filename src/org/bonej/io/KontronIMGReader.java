package org.bonej.io;

import ij.IJ;
import ij.ImagePlus;
import ij.io.FileInfo;
import ij.io.OpenDialog;
import ij.plugin.PlugIn;
import ij.process.ByteProcessor;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;

import org.doube.util.UsageReporter;

/**
 * This plugin implements the Import > Kontron IMG command.
 * 
 * @author Michael Doube, RVC, London, UK.
 */
public class KontronIMGReader implements PlugIn {

	/** All IMG start with 01 00 47 12 6D B0  (hex view) */
	private static final byte[] MAGIC = { 1, 0, 71, 18, 109, -80 };
	
	/** 128-byte header */
	private static final int HEADER_LENGTH = 128;

	private int eofErrorCount;

	public void run(String arg) {

		// the IMG-File is selected
		OpenDialog od = new OpenDialog("Open IMG...", arg);
		String directory = od.getDirectory();
		String fileName = od.getFileName();
		String path = directory + fileName;
		if (fileName == null)
			return;
		if (!isKontronIMG(path)) {
			IJ.error("IMG Reader",
					"Not an IMG file. Magic number does not match.");
			return;
		}

		// Open the file
		try {
			ImagePlus imp = openKontronIMG(path);
			imp.show();
			UsageReporter.reportEvent(this).send();
		} catch (IllegalArgumentException e) {
			IJ.error("IMG Reader", e.getMessage());
			return;
		}
	}

	/**
	 * Opens a Kontron IMG file as an ImageJ ImagePlus
	 * 
	 * @param path
	 * @return
	 */
	public ImagePlus openKontronIMG(String path) {

		int[] imageSize = getImageSize(path);
		int width = imageSize[0];
		int height = imageSize[1];	

		// FileInfo
		FileInfo fi = new FileInfo();
		fi.fileName = new File(path).getName();
		fi.directory = new File(path).getParent()
				+ ((IJ.isWindows()) ? "\\" : "/");
		fi.width = width;
		fi.height = height;

		fi.intelByteOrder = true;
		fi.whiteIsZero = false;
		fi.fileType = FileInfo.GRAY8;

		byte[] pixels = new byte[width * height];
		try {
			FileInputStream is = new FileInputStream(path);
			pixels = readPixels(is, HEADER_LENGTH, width, height);
			is.close();
		} catch (Exception e) {
			IJ.log("" + e);
		} catch (OutOfMemoryError e) {
			IJ.outOfMemory(fi.fileName);
		}
		ByteProcessor bp = new ByteProcessor(width, height, pixels);
		ImagePlus imp = new ImagePlus(fi.fileName, bp);
		return imp;
	}

	/** *********************************************************************** **/
	/**
	 * from ImageReader.java: Skips the specified number of bytes, then reads an
	 * image and returns the pixel array (byte, short, int or float). Returns
	 * null if there was an IO exception. Does not close the InputStream.
	 */
	private byte[] readPixels(FileInputStream in, long skipCount, int width,
			int height) {
		byte[] pixels = readPixels(in, width, height);
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
	private byte[] readPixels(FileInputStream in, int width, int height) {
		try {
			byte[] pixels = new byte[width * height]; 
			in.read(pixels, HEADER_LENGTH, pixels.length);
			return pixels;
			
		} catch (IOException e) {
			IJ.log("" + e);
			return null;
		}
	}

	/**
	 * Check the magic number to determine if a file is a Kontron IMG
	 * 
	 * @param path
	 * @return true if the file is a Scanco ISQ
	 */
	public boolean isKontronIMG(String path) {
		return getMagic(path).equals(MAGIC);
	}

	/**
	 * Get the magic number from a file
	 * 
	 * @param path
	 * @return an array containing the first six bytes of the file
	 */
	public byte[] getMagic(String path) {
		return readBytes(path, 0, 6);
	}

	public int getPatientIndex(String path) {
		return readInt(path, 28);
	}

	/**
	 * Get the pixel size of the image
	 * 
	 * @param path
	 * @return {x, y} image size in pixels
	 */
	public int[] getImageSize(String path) {
		int[] sizes = { readInt(path, 6), readInt(path, 8) };
		return sizes;
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


	private byte[] readBytes(String path, int firstByte, int length) {
		if (path == null)
			throw new IllegalArgumentException();
		try {
			File iFile = new File(path);
			FileInputStream p = new FileInputStream(iFile);
			p.skip(firstByte);
			byte[] bytes = new byte[length];
			for (int kh = 0; kh < length; kh++) {
				bytes[kh] = (byte) p.read();
			}
			p.close();
			return bytes;
		} catch (IOException e) {
			IJ.handleException(e);
		}
		return null;
	}

}