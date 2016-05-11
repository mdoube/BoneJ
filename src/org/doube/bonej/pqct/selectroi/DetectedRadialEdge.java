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

public class DetectedRadialEdge implements Comparable<DetectedRadialEdge>{
	public int ii;		//indexes for x-coordinates
	public int jj;	//indexes for y-coordinates
	public double theta;
	public double radius;

	public DetectedRadialEdge(int ii,int jj, double theta, double radius){
		this.ii = ii;
		this.jj = jj;
		this.theta = theta;
		this.radius = radius;
	}
	
	@Override
	public int compareTo(DetectedRadialEdge o){
		int returnValue = 0;
		if (o == null || this == null) {throw new NullPointerException();}
		if (this.theta == o.theta) {return 0;}
		return this.theta < o.theta ? -1 : 1;		
	}	
}
