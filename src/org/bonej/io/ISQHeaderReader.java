package org.bonej.io;

/*
 * KHKs_Scanco_ISQ_HeaderReader opens a Scanco ISQ File and reads the header info
 *
 * We use it to display the information saved as "patient name".
 * The variable patient name is used to describe all our projects.
 *
 */

/**
 *
 * @author Prof. Dr. Karl-Heinz Kunzelmann
 * @author karl-heinz@kunzelmann.de
 */

import ij.plugin.PlugIn;
import java.io.*;
import ij.*;
import ij.io.*;

/** This plugin implements the an header reader for Scanco ISQ files command. */
public class ISQHeaderReader implements PlugIn {

	int i, offset, offset1, offset2, offset3, xdimension, ydimension,
			zdimension, mu_scaling;

	int tmpInt;
	String nameStringInHeader = "";
	float el_size_mm_x, el_size_mm_y, el_size_mm_z;

	public void run(String arg) {

		OpenDialog od = new OpenDialog("Open ISQ Header...", arg);
		String directory = od.getDirectory();
		String fileName = od.getFileName();
		if (fileName == null)
			return;

		try {
			File iFile = new File(directory + fileName);
			FileInputStream p = new FileInputStream(iFile);

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

			p.skip(20);
			// 82 int mu_scaling;
			tmpInt = (p.read() + p.read() * 256 + p.read() * 65536 + p.read() * 256 * 65536);
			mu_scaling = tmpInt;

			p.skip(36); // war 60, wegen mu_scaling jetzt 36

			for (int kh = 0; kh < 40; kh++) {
				char ch = (char) p.read();
				nameStringInHeader = nameStringInHeader + ch;
				// System.out.println(nameStringInHeader);
			}

			p.skip(340);

			offset = (p.read() + p.read() * 256 + p.read() * 65536 + 1) * 512;
			// System.out.println("offset from header: "+offset);

			IJ.showMessage("Header Data:\n " + "\ninfo: " + nameStringInHeader
					+ "\nmu_scaling: " + mu_scaling + "\n..................."
					+ "\nheader byte offset: " + offset + "\nx dimension: "
					+ xdimension + "\ny dimension: " + ydimension
					+ "\nz dimension: " + zdimension + "\nel_size x (in mm): "
					+ el_size_mm_x + "\nel_size y (in mm): " + el_size_mm_y
					+ "\nel_size z (in mm): " + el_size_mm_z);

			p.close();

		} catch (IOException e) {
			System.out.println("IOException error!" + e.getMessage());
		}

	}

}

/*
 * Scanco ISQ Header Information:
 * 
 * typedef struct { /*--------------------------------------------- 00 char
 * check[16]; // Char is in Java 2 Byte 16 int data_type; // Int = 4 Byte 20 int
 * nr_of_bytes; /* either one of them 24 int nr_of_blocks; /* or both, but min.
 * of 1 28 int patient_index; /* 1 block = 512 bytes 32 int scanner_id; 36 int
 * creation_date[2]; /*--------------------------------------------- 40 int
 * dimx_p; 44 int dimy_p; 48 int dimz_p; 52 int dimx_um; 56 int dimy_um; 60 int
 * dimz_um; 64 int slice_thickness_um; 68 int slice_increment_um; 72 int
 * slice_1_pos_um; 76 int min_data_value; 78 int max_data_value; 82 int
 * mu_scaling; /* p(x,y,z)/mu_scaling = value [1/cm] 86 int nr_of_samples; 90
 * int nr_of_projections; 94 int scandist_um; 98 int scanner_type; 102 int
 * sampletime_us; 106 int index_measurement; 110 int site; /* Coded value 114
 * int reference_line_um; 120 int recon_alg; /* Coded value 124 char name[40];
 * // KHK char evtl. nur 1 Byte... siehe unten ?????? int energy; /*V int
 * intensity; /* uA int fill[83];
 * /*--------------------------------------------- int data_offset; /* in
 * 512-byte-blocks } ima_data_type, *ima_data_typeP;
 * 
 * So the first 16 bytes are a string 'CTDATA-HEADER_V1', used to identify the
 * type of data. The 'int' are all 4-byte integers.
 * 
 * dimx_p is the dimension in pixels, dimx_um the dimension in microns.
 * 
 * So dimx_p is at byte-offset 40, then dimy_p at 44, dimz_p (=number of slices)
 * at 48.
 * 
 * The microCT calculates so called 'x-ray linear attenuation' values. These
 * (float) values are scaled with 'mu_scaling' (see header, e.g. 4096) to get to
 * the signed 2-byte integers values that we save in the .isq file.
 * 
 * e.g. Pixel value 8192 corresponds to lin. att. coeff. of 2.0 [1/cm]
 * (8192/4096)
 * 
 * Following to the headers is the data part. It is in 2-byte short integers
 * (signed) and starts from the top-left pixel of slice 1 to the left, then the
 * next line follows, until the last pixel of the last sclice in the lower
 * right.
 */