package org.doube.geometry;

import java.util.ArrayList;
import java.util.List;

import javax.vecmath.Point3f;

public class TestDataMaker {

	static public List<Point3f> octahedron(){
		final float a = 2.0f;
		final float b = 1.198039f;
		final float c = 2.801961f;
		List<Point3f> points = new ArrayList<Point3f>();
		/*
		 solid untitled.stl
		 facet normal 5.773502E-01 5.773502E-01 5.773502E-01
		  outer loop
		   vertex 2.000000E+00 2.000000E+00 1.198039E+00
		   vertex 1.198039E+00 2.000000E+00 2.000000E+00
		   vertex 2.000000E+00 1.198039E+00 2.000000E+00
		  endloop
		 endfacet */
		 points.add(new Point3f(a, a, b));
		 points.add(new Point3f(b, a, a));
		 points.add(new Point3f(a, b, a));
		 /*
		 facet normal 5.773503E-01 -5.773503E-01 5.773503E-01
		  outer loop
		   vertex 2.000000E+00 2.000000E+00 1.198039E+00
		   vertex 2.000000E+00 2.801961E+00 2.000000E+00
		   vertex 1.198039E+00 2.000000E+00 2.000000E+00
		  endloop
		 endfacet
		 */
		 points.add(new Point3f(a, a, b));
		 points.add(new Point3f(a, a, a));
		 points.add(new Point3f(b, a, a));
		 /*
		 facet normal -5.773503E-01 5.773503E-01 5.773503E-01
		  outer loop
		   vertex 2.000000E+00 1.198039E+00 2.000000E+00
		   vertex 2.801961E+00 2.000000E+00 2.000000E+00
		   vertex 2.000000E+00 2.000000E+00 1.198039E+00
		  endloop
		 endfacet
		 */
		 points.add(new Point3f(a, b, a));
		 points.add(new Point3f(c, a, a));
		 points.add(new Point3f(a, a, b));
		 /*
		 facet normal -5.773503E-01 -5.773503E-01 5.773502E-01
		  outer loop
		   vertex 2.801961E+00 2.000000E+00 2.000000E+00
		   vertex 2.000000E+00 2.801961E+00 2.000000E+00
		   vertex 2.000000E+00 2.000000E+00 1.198039E+00
		  endloop
		 endfacet
		 */
		 points.add(new Point3f(c, a, a));
		 points.add(new Point3f(a, c, a));
		 points.add(new Point3f(a, a, b));
		 /*
		 facet normal 5.773503E-01 5.773503E-01 -5.773503E-01
		  outer loop
		   vertex 2.000000E+00 1.198039E+00 2.000000E+00
		   vertex 1.198039E+00 2.000000E+00 2.000000E+00
		   vertex 2.000000E+00 2.000000E+00 2.801961E+00
		  endloop
		 endfacet
		 */
		 points.add(new Point3f(a, b, a));
		 points.add(new Point3f(b, a, a));
		 points.add(new Point3f(a, a, c));
		 /*
		 facet normal 5.773502E-01 -5.773503E-01 -5.773503E-01
		  outer loop
		   vertex 1.198039E+00 2.000000E+00 2.000000E+00
		   vertex 2.000000E+00 2.801961E+00 2.000000E+00
		   vertex 2.000000E+00 2.000000E+00 2.801961E+00
		  endloop
		 endfacet
		 */
		 points.add(new Point3f(b, a, a));
		 points.add(new Point3f(a, c, a));
		 points.add(new Point3f(a, a, c));
		 /*
		 facet normal -5.773503E-01 5.773502E-01 -5.773503E-01
		  outer loop
		   vertex 2.000000E+00 1.198039E+00 2.000000E+00
		   vertex 2.000000E+00 2.000000E+00 2.801961E+00
		   vertex 2.801961E+00 2.000000E+00 2.000000E+00
		  endloop
		 endfacet
		 */
		 points.add(new Point3f(a, c, a));
		 points.add(new Point3f(a, a, c));
		 points.add(new Point3f(c, a, a));
		 /*
		 facet normal -5.773503E-01 -5.773503E-01 -5.773503E-01
		  outer loop
		   vertex 2.801961E+00 2.000000E+00 2.000000E+00
		   vertex 2.000000E+00 2.000000E+00 2.801961E+00
		   vertex 2.000000E+00 2.801961E+00 2.000000E+00
		  endloop
		 endfacet
		 */
		 points.add(new Point3f(c, a, a));
		 points.add(new Point3f(a, a, c));
		 points.add(new Point3f(a, c, a));
		 /*
		  endsolid untitled.stl
		  */
		return null;
		
	}
	
}
