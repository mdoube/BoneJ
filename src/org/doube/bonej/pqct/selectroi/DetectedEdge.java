/*
	This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

	N.B.  the above text was copied from http://www.gnu.org/licenses/gpl.html
	unmodified. I have not attached a copy of the GNU license to the source...

    Copyright (C) 2011 Timo Rantalainen
*/

package org.doube.bonej.pqct.selectroi;
import java.util.*;	//Vector, Collections
import java.lang.Math; //atan2
import java.awt.*;			//Polygon, Rectangle
import org.doube.bonej.pqct.io.*;	//image data
import ij.*;		//ImagePlus
import ij.gui.*;	//ImagePlus ROI
import ij.text.*; 	//Debugging ...
import ij.process.*;	//Debugging
@SuppressWarnings(value ={"serial","unchecked"}) //Unchecked for obtaining Vector<Object> as a returnvalue

public class DetectedEdge implements Comparable<DetectedEdge>{
	public Vector<Integer> iit;		//indexes for x-coordinates
	public Vector<Integer> jiit;	//indexes for y-coordinates
	public int area;
	public int length;
	
	public DetectedEdge(Vector<Integer> iit,Vector<Integer> jiit,int area){
		this.iit = iit;
		this.jiit = jiit;
		this.length = iit.size();
		this.area = area;
	}
	
	public int compareTo(DetectedEdge o){
		int returnValue = 0;
		if (o == null || this == null) {throw new NullPointerException();}
		if (this.area == o.area) {return 0;}
		return this.area < o.area ? -1 : 1;		
	}

	
}
