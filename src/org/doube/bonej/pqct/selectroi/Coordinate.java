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
public class Coordinate{
	public int ii;
	public int jj;
	/*Constructors*/
	public Coordinate(){
		this(-1,-1);
	}
	
	public Coordinate(int ii,int jj){
		this.ii = ii;
		this.jj = jj;
	}
	
	public Coordinate subtract(Coordinate a){
		return new Coordinate(this.ii-a.ii,this.jj-a.ii);
	}
	

	
	public int maxVal(){
		return max(Math.abs(ii),Math.abs(jj));
	}
	
	private int max(int a,int b){
		return a >= b ? a:b;
	}
	
	
}