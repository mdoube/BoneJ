package org.bonej.io;

/*
 * ImageJ comes already with some useful routines to read the Shimadzu µCT files. However, the files will not contain the necessary information on the resolution and the voxel dimensions. To attach these information to an ImageJ file you will need the plugin “KHKs Shimadzu µCT HeaderfileReader”. This plugin will get a subset of information from the header files which are associated with the Shimadzu µCT datafiles. The data files have the extension *.tif, while the header files have the extension *.inf. Please take care to leave the *.inf files in the same directory as the data files.
 * To import the Shimadzu µCT files please use the Files > Import > Image Sequence
 * And enter in the following window:
 *
 * -	the number of images you would like to open
 * -	the starting image
 * -	and “tif” in the “File Name Contains” field.
 *
 * After you click “OK” you will import the tif-files into ImageJ and you will get what is called an image stack. This stack displays the images, but it does not contain the scaling information which is necessary for the subsequent evaluation.
 *
 * To import the scaling information you have to open Plugins > KHKs: Shimadzu µCT HeaderfileReader
 *
 * The action of this plugin is rather unspectacular. If you have the *.inf files in the same directory as the *.tif files, then you will see the following “success message”.
 *
 *
 * This file reads the header of the Shimadzu MicroCT / Univ. Tokyo.
 * this special version of the file does not import the package LEDataStream
 * it just neads the file: 
 * 
 * LEDataInputStream.java in the same directory as this file is stored.
 * 
 * anything else is comparable.
 * 
 * little endian import routines: copyright (c) 1999-2008 Roedy Green, Canadian Mind Products, http://mindprod.com";
 * 
 * Immer wieder das gleiche Problem.
 * In Java sind alle Datentypen ausser Char signed.
 * In Java ist Datentyp Char 2 Byte lang, in C nur 1 Byte, 
 * Datentype Char in C entspricht Byte in Java. Allerdings muss noch das Vorzeichen beachtet werden.
 * In C Char von 0 - 256 in Java von -127 bis + 127.
 * In JAVA, a byte always considered as signed when converted to another type. We must mask the sign bit to JAVA, cast to an integer and process the masked bit if needed. The following method implements this idea :
 * 
 * Intel = Little Endian, kleines Ende zuerst
 * Ich nehme einfach an, dass wir immer Little Endian haben
 * Dummerweise arbeitet Java immer mit Big Endian... also umsortieren

 */

// TODO:
// Idee: in 1 Makro 3 Schritte - 
// FolderOpener.java
// hier prüfen, ob man aus Titel und Slice den Pfad für das inf generieren kann
// ImageProperities.java 
// hier Abfrage mit weitermachen
// Process Math Subtract - Air-Value to scale.
// elegant wäre später mit randomaccessstream die Werte für den header auszulesen
// 

/**
 *
 * 
 * @author Prof. Dr. Karl-Heinz Kunzelmann
 * @author www.kunzelmann.de
 * 
 *
 *  little endian import routines: copyright (c) 1999-2008 Roedy Green, Canadian Mind Products, http://mindprod.com";
 * 
 */

/***********************************************************************************
 * File header of Shimadzu microCT files
 *
 * The data are saved in little endian byte format.
 * This makes it really complicated to import them into java
 * Especially the double varibles.
 * Java used big endian file format.
 *
 * In this case ...
 *
 * char = 1 Byte
 * Long = 4 Byte - in contract: long in Java is 8 Byte
 * Double = 8 Byte
 *
 * ... long.
 *
 * 2004.11.29 SHIMADZU corporation NDIBU
 CT data *.inf file structure
 No. 	Type 	Name	Application
 1 	char 	cFileCreateYMD[10] 	Created date and time
 2 	char 	c2xxxx[6]
 3 	char 	cFolderName[260] 	Folder name
 4 	char 	cFileName[256] 	File name
 5 	long 	lDispSizeX 	Lengthwise size of the image
 6 	long 	lDispSizeY 	Crosswise size of the image
 7 	long 	lXrayVol 	X-ray tube voltage
 8 	long 	lXrayCur 	Tube current
 9 	double 	dIISize 	Sensor size
 10 	double 	dSID 	??? value
 11 	double 	dSOD 	??? value
 12 	double 	d12xxxx
 13 	double 	d13xxxx
 14 	long 	l14xxxx
 15 	long 	l15xxxx
 16 	long 	l16xxxx
 17 	char 	c17xxxx[10]
 ** 	char 	c_dummy[6]
 18 	long 	l18xxxx
 19 	long 	l19xxxx
 20 	long 	l20xxxx
 21 	double 	dSliceWidth 	Thickness of the slice
 22 	long 	lViews 	view number
 23 	long 	lAveCnt 	Average number
 24 	double 	dImageConvC1 	Scaling coefficient
 25 	char 	c25xxxx[32]
 26 	long 	lCTMode 	CT mode 1 0:2DCT 1:corn CT
 27 	long 	lOffsetScan 	CT mode 2 0:Normal 1:Offset
 28 	BOOL 	bHalfScan 	CT mode 3 0:Full Scan 1:Half scan
 29 	double 	dFovXY 	FOV (XY)
 30 	double 	dFovZ 	FOV (Z)
 31 	double 	dPitchX 	Pixel equivalent length X
 32 	double 	dPitchY 	Pixel equivalent length Y
 33 	double 	dSlicePitch 	Slice pitch
 34 	long 	l34xxxx
 35 	BOOL 	bReverseImage 	0:BottomView 1:TopView
 36 	long 	l36xxxx
 37 	double 	d37xxxx
 38 	double 	d38xxxx
 39 	double 	d39xxxx
 40 	short 	n40xxxx
 41 	double 	d41xxxx
 42 	double 	d42xxxx
 43 	double 	d43xxxx
 44 	short 	n44xxxx
 45 	double 	d45xxxx
 46 	double 	d46xxxx
 47 	char 	c47xxxx[400];
 ***********************************************************************************
 */

import ij.plugin.filter.PlugInFilter;
import ij.*;
import ij.io.*;
import ij.gui.*;
import ij.process.*;
import ij.measure.*;
import java.io.*;
import java.util.*;
import java.text.DecimalFormat;

public class ShimadzuReader implements PlugInFilter {

	private String path;
	static ImagePlus orgImg;
	FileInfo fi;
	Calibration cal;
	Calibration calOrg;
	boolean infFileFound = true;
	String label;
	String defaultDir;

	public int setup(String arg, ImagePlus imp) {
		orgImg = imp;
		fi = imp.getFileInfo();
		cal = imp.getCalibration();
		calOrg = cal.copy();

		// hiermit bekomme ich den vollständigen Pfad inkl. \ am Ende
		defaultDir = OpenDialog.getDefaultDirectory();
		System.out.println("defDir: " + defaultDir);

		ImageStack stack = imp.getStack();
		label = stack.getSliceLabel(1);
		System.out.println("slicelabel: " + label);

		fi.directory = defaultDir;
		fi.fileName = label;

		path = defaultDir + label;

		return DOES_8G + DOES_16 + DOES_32;
	}

	public void run(ImageProcessor processor) {

		checkFile(path);

		System.out.println("boolean: " + infFileFound);

		if (infFileFound == false) {
			negativeFeedback();
			return;
		}

		String infFileName = labelTif2Inf(label);

		fi.unit = "mm";
		cal.setUnit("mm");

		path = defaultDir + infFileName;
		System.out.println("path: " + path);

		// O P E N inf-File

		FileInputStream fis;
		try {
			fis = new FileInputStream(path);
		} catch (FileNotFoundException e) {
			IJ.log("Inf-File not found" + e);
			System.out.println("Inf-File not found - fis exception");
			return;
		}

		LEDataInputStream ledis = new LEDataInputStream(fis);
		try {
			for (int i = 0; i < 704; i++) {
				// erst wird das signed Byte in ein Int gecastet, dann wird das
				// int in Char gecastet.
				ledis.readByte();
			}

			double d1 = ledis.readDouble();
			System.out.println("dFovXY: " + d1);

			double d2 = ledis.readDouble();
			System.out.println("dFovZ: " + d2);

			double d3 = ledis.readDouble();
			System.out.println("dPitchX: " + d3);

			double d4 = ledis.readDouble();
			System.out.println("dPitchY: " + d4);

			double d5 = ledis.readDouble();
			System.out.println("dSlicePitch: " + d5);

			ledis.close();

			cal.pixelWidth = d3;
			cal.pixelHeight = d4;
			cal.pixelDepth = d5;

		} catch (IOException e) {
			IJ.log("Error opening LEDataInputStream" + e);
			System.out
					.println("Unexpected IOException opening LEDataInputStream");
			return;
		}

		orgImg.setCalibration(cal);
		positiveFeedback();
	}

	void negativeFeedback() {

		GenericDialog gd = new GenericDialog(
				"Prof. Kunzelmann's: Import Shimadzu MicroCT Files");

		gd.addMessage("Header-file could not be found");
		gd.addMessage("Please check that it is in the same directory");
		gd.addMessage("as the Tif-File sequence");

		gd.showDialog();
		if (gd.wasCanceled()) {
		}
	}

	void positiveFeedback() {

		GenericDialog gd = new GenericDialog(
				"Prof. Kunzelmann's: Import Shimadzu MicroCT Files");

		String text1 = "pixelSize in x: " + formatDouble(cal.pixelWidth) + " "
				+ fi.unit;
		String text2 = "pixelSize in y: " + formatDouble(cal.pixelHeight) + " "
				+ fi.unit;
		String text3 = "pixelSize in z: " + formatDouble(cal.pixelDepth) + " "
				+ fi.unit;

		gd.addMessage("Header-file successfully read!");
		gd.addMessage("Click 'OK', if the following data are correct!");
		gd.addMessage(text1);
		gd.addMessage(text2);
		gd.addMessage(text3);

		gd.showDialog();
		if (gd.wasCanceled()) {

			cal = calOrg.copy();
			orgImg.setCalibration(cal);
			System.out.println("Original calibration info restored");

		}
	}

	boolean checkFile(String path) {
		File infFile = new File(path);
		if (infFile.exists() == true) {
			System.out.println("inf-File found at " + path);
			return infFileFound = true;
		} else {
			System.out.println("no inf-File found at " + fi.directory);
			return infFileFound = false;
		}
	}

	String labelTif2Inf(String orgLabel) {
		String oldLabel = orgLabel;

		String newLabel = oldLabel.replace(".tif", ".inf");
		System.out.println("Methode: labelTif2Inf: " + newLabel);
		return newLabel;
	}

	public String formatDouble(double in) {

		Locale.setDefault(Locale.US);

		DecimalFormat format = new DecimalFormat("#,##0.0000");
		return format.format(in);

	}
}
