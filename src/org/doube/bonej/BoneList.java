package org.doube.bonej;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

import ij.ImagePlus;

/**
 * Utility methods for making lists of bones and guessing which bone is in an
 * image
 *
 * @author Michael Doube
 *
 */
public class BoneList {

	/**
	 * List of bone names
	 */
	// Only add new bone names to the END of this list
	private final static String[] boneList = { "unknown", "scapula", "humerus", "radius", "ulna", "metacarpal",
			"pelvis", "femur", "tibia", "fibula", "metatarsal", "calcaneus", "tibiotarsus", "tarsometatarsal",
			"sacrum" };

	/**
	 * Return the array of bone names
	 *
	 * @return array of bone names
	 */
	public final static String[] get() {
		return boneList;
	}

	/**
	 * Guess from the image title which bone is in the image
	 *
	 * @param imp
	 * @return integer code relating to the position of the bone's name in the
	 *         bone list
	 */
	public static int guessBone(final ImagePlus imp) {
		final String boneString = imp.getTitle();
		return guessBone(boneString);
	}

	/**
	 * Return the boneID of a bone in boneList that matches the input string
	 *
	 * @param boneString
	 * @return
	 */
	public static int guessBone(final String boneString) {
		int boneID = 0;
		for (int n = 0; n < boneList.length; n++) {
			final Pattern p = Pattern.compile(boneList[n], Pattern.CASE_INSENSITIVE);
			final Matcher m = p.matcher(boneString);
			if (m.find()) {
				boneID = n;
				continue;
			}
		}
		return boneID;
	}
}
