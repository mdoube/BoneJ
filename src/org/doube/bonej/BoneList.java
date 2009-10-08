package org.doube.bonej;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

import ij.ImagePlus;

public class BoneList {
	/**
	 * Return an array of bone names
	 * 
	 * @return array of bone names
	 */
	public String[] getBoneList(){
		String[] boneList = {
			"unknown",
			"scapula",
			"humerus",
			"radius",
			"ulna",
			"metacarpal",
			"pelvis",
			"femur",
			"tibia",
			"fibula",
			"metatarsal"
			};
		return boneList;
	}
	
	/**
	 * Guess from the image title which bone is in the image
	 * 
	 * @param imp
	 * @return integer code relating to the position of the bone's name in the bone list
	 */
	public int guessBone(ImagePlus imp){
		String[] bones = getBoneList();
		String title = imp.getTitle();
		int boneID = 0;
		for (int n = 0; n < bones.length; n++) {
			Pattern p = Pattern.compile(bones[n], Pattern.CASE_INSENSITIVE);
			Matcher m = p.matcher(title);
			if (m.find()) {
				boneID = n;
				continue;
			}
		}
		return boneID;
	}
}
