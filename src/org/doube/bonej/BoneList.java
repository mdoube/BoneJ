package org.doube.bonej;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

import ij.ImagePlus;

/**
 * Utility methods for making lists of bones and guessing which bone is in an image 
 * 
 * @author Michael Doube
 *
 */
public class BoneList {
	/**
	 * Return an array of bone names
	 * 
	 * @return array of bone names
	 */
	public final static String[] getBoneList(){
		//ONLY ADD TO THE END OF THIS ARRAY
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
			"metatarsal",
			"calcaneus"
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
		String boneString = imp.getTitle();
		return guessBone(boneString);
	}
	
	/**
	 * Return the boneID of a bone in boneList that matches the input string
	 * 
	 * @param boneString
	 * @return
	 */
	public int guessBone(String boneString){
		String[] bones = getBoneList();
		int boneID = 0;
		for (int n = 0; n < bones.length; n++) {
			Pattern p = Pattern.compile(bones[n], Pattern.CASE_INSENSITIVE);
			Matcher m = p.matcher(boneString);
			if (m.find()) {
				boneID = n;
				continue;
			}
		}
		return boneID;
	}
}
