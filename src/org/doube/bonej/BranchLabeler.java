package org.doube.bonej;

/**
 * BranchLabeler plugin for ImageJ
 * Copyright 2010 Michael Doube 
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

import org.doube.util.ImageCheck;

import ij.IJ;
import ij.ImagePlus;
import ij.plugin.PlugIn;

/**
 * 
 * @author Michael Doube
 * 
 */
public class BranchLabeler implements PlugIn {

	public void run(String arg) {
		if (!ImageCheck.checkIJVersion())
			return;
		ImagePlus imp = IJ.getImage();
		ImageCheck ic = new ImageCheck();
		if (!ic.isBinary(imp)) {
			IJ.error("Branch Labeler requires a binary image.");
			return;
		}
		Skeletonize3D skel = new Skeletonize3D();
		ImagePlus skeleton = skel.getSkeleton(imp);
		
		//Label individual branches of skeleton...
		
		//create int working stack from labelled branches
		
		//iteratively grow neighbourhoods from labelled branches
		//until all the space in imp is filled
		
		//return an imp containing the labelled branches

	}
}
