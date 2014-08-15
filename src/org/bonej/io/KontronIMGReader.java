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

	/** All IMG start with 01 00 47 12 6D B0 (hex view) */
	private static final int MAGIC = 306642945;

	/** 128-byte header */
	private static final int HEADER_LENGTH = 128;

	public void run(String arg) {

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

		ByteProcessor bp = new ByteProcessor(width, height, readBytes(path,
				HEADER_LENGTH, width * height));
		ImagePlus imp = new ImagePlus(fi.fileName, bp);
		return imp;
	}

	/**
	 * Check the magic number to determine if a file is a Kontron IMG
	 * 
	 * @param path
	 * @return true if the file is a Scanco ISQ
	 */
	public boolean isKontronIMG(String path) {
		return getMagic(path) == MAGIC;
	}

	/**
	 * Get the magic number from a file
	 * 
	 * Although the first six bytes appear to be invariant, the first 4 should
	 * suffice, and pack neatly into a single int
	 * 
	 * @param path
	 * @return an int made of the first 4 byes of the file
	 */
	public int getMagic(String path) {
		return readInt(path, 0);
	}

	/**
	 * Get the pixel size of the image
	 * 
	 * @param path
	 * @return {x, y} image size in pixels
	 */
	public int[] getImageSize(String path) {
		int[] sizes = { readShort(path, 6), readShort(path, 8) };
		return sizes;
	}

	private byte[] readBytes(String path, int firstByte, int length) {
		if (path == null)
			throw new IllegalArgumentException();
		try {
			File iFile = new File(path);
			FileInputStream p = new FileInputStream(iFile);
			p.skip(firstByte);
			byte[] bytes = new byte[length];
			p.read(bytes, 0, length);
			p.close();
			return bytes;
		} catch (IOException e) {
			IJ.handleException(e);
		}
		return null;
	}
	
	private int readShort(String path, int firstByte) {
		if (path == null)
			throw new IllegalArgumentException();
		try {
			File iFile = new File(path);
			FileInputStream p = new FileInputStream(iFile);
			p.skip(firstByte);
			int value = (p.read() + p.read() * 256);
			p.close();
			return value;
		} catch (IOException e) {
			IJ.handleException(e);
		}
		return -1;
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
}