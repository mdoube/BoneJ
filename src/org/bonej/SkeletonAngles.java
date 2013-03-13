package org.bonej;

/**
 * SkeletonAngles class for ImageJ
 * Copyright 2012 Michael Doube 
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

import java.util.ArrayList;

import org.doube.geometry.Centroid;
import org.doube.geometry.Trig;
import org.doube.skeleton.AnalyzeSkeleton;
import org.doube.skeleton.Edge;
import org.doube.skeleton.Graph;
import org.doube.skeleton.Point;
import org.doube.skeleton.Vertex;
import org.doube.util.ResultInserter;
import org.doube.util.UsageReporter;

import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;

public class SkeletonAngles implements PlugIn {

	/** Measure angles between vertices */
	public static final int VERTEX_TO_VERTEX = -1;

	public void run(String arg) {
		ImagePlus imp = WindowManager.getCurrentImage();
		if (imp == null) {
			IJ.noImage();
			return;
		}

		GenericDialog gd = new GenericDialog("Options");
		String[] options = { "Opposite vertex", "nth edge pixel" };
		gd.addChoice("Use:", options, options[0]);
		gd.addNumericField("nth pixel", 0, 0);
		gd.addHelp("http://bonej.org/triplepointangles");

		gd.showDialog();

		if (gd.wasCanceled())
			return;

		int choice = gd.getNextChoiceIndex();

		int nthPixel = (int) Math.floor(gd.getNextNumber());

		if (choice == 0) {
			nthPixel = VERTEX_TO_VERTEX;
		}

		ResultInserter ri = ResultInserter.getInstance();

		double[][][] result = calculateTriplePointAngles(imp, nthPixel);

		for (int g = 0; g < result.length; g++) {
			for (int v = 0; v < result[g].length; v++) {
				if (result[g][v] == null)
					continue;
				ri.setResultInRow(imp, "Skeleton", g);
				ri.setResultInRow(imp, "Vertex", v);
				ri.setResultInRow(imp, "Theta 0", result[g][v][0]);
				ri.setResultInRow(imp, "Theta 1", result[g][v][1]);
				ri.setResultInRow(imp, "Theta 2", result[g][v][2]);
			}
		}
		ri.updateTable();
		UsageReporter.reportEvent(this).send();
	}

	/**
	 * Calculate the three angles formed by the branches at every triple point
	 * in the skeleton. Assumes the image has been skeletonized. Angles can be
	 * calculated between opposite vertices or an arbitrary number of points
	 * from the triple point.
	 * 
	 * @param imp
	 *            Skeletonized image with the skeleton in the foreground
	 * @param nthPixel
	 *            Number of points along the edge away from the triple point to
	 *            use for angle calculation. Set to
	 *            SkeletonAngles.VERTEX_TO_VERTEX to use opposite vertices
	 * @return 3D array containing 3 angles (in radians) for each triple point
	 *         in each skeleton in the image
	 */
	public double[][][] calculateTriplePointAngles(ImagePlus imp, int nthPixel) {
		AnalyzeSkeleton skeletonAnalyzer = new AnalyzeSkeleton();
		skeletonAnalyzer.setup("", imp);
		skeletonAnalyzer.run();
		Graph[] graphs = skeletonAnalyzer.getGraphs();
		double[][][] angleList = new double[graphs.length][][];
		int g = 0;
		for (Graph graph : graphs) {
			ArrayList<Vertex> vertices = graph.getVertices();
			double[][] verts = new double[vertices.size()][3];
			int v = 0;
			for (Vertex vertex : vertices) {
				// triple point
				if (vertex.getBranches().size() == 3) {
					ArrayList<Edge> edges = vertex.getBranches();
					if (edges.size() != 3) {
						try {
							throw new Exception("Unexpected number of edges: "
									+ edges.size());
						} catch (Exception e) {
							e.printStackTrace();
						}
					}
					Edge edge0 = edges.get(0);
					Edge edge1 = edges.get(1);
					Edge edge2 = edges.get(2);

					double theta0 = vertexAngle(vertex, edge0, edge1, nthPixel);
					double theta1 = vertexAngle(vertex, edge0, edge2, nthPixel);
					double theta2 = vertexAngle(vertex, edge1, edge2, nthPixel);

					double[] thetas = { theta0, theta1, theta2 };

					verts[v] = thetas;
				} else {
					verts[v] = null;
				}
				v++;
			}
			angleList[g] = verts;
			g++;
		}
		return angleList;
	}

	private double vertexAngle(Vertex vertex, Edge edge0, Edge edge1) {
		Vertex v0 = edge0.getOppositeVertex(vertex);
		Vertex v1 = edge1.getOppositeVertex(vertex);
		ArrayList<Point> pointsv = vertex.getPoints();
		ArrayList<Point> points0 = v0.getPoints();
		ArrayList<Point> points1 = v1.getPoints();
		double[] cv = Centroid.getCentroid(pointsv);
		double[] c0 = Centroid.getCentroid(points0);
		double[] c1 = Centroid.getCentroid(points1);

		return Trig.angle3D(c0[0], c0[1], c0[2], c1[0], c1[1], c1[2], cv[0],
				cv[1], cv[2]);
	}

	private double vertexAngle(Vertex vertex, Edge edge0, Edge edge1,
			int nthPoint) {
		if (nthPoint == VERTEX_TO_VERTEX)
			return vertexAngle(vertex, edge0, edge1);
		Point p0 = getNthPoint(vertex, edge0, nthPoint);
		Point p1 = getNthPoint(vertex, edge1, nthPoint);
		double[] cv = Centroid.getCentroid(vertex.getPoints());
		return Trig.angle3D(p0.x, p0.y, p0.z, p1.x, p1.y, p1.z, cv[0], cv[1],
				cv[2]);
	}

	private Point getNthPoint(Vertex vertex, Edge edge, int nthPoint) {
		ArrayList<Point> vertexPoints = vertex.getPoints();
		ArrayList<Point> edgePoints = edge.getSlabs();
		boolean startAtZero = false;
		outerloop: for (Point v : vertexPoints) {
			Point p0 = edgePoints.get(0);
			for (int x = v.x - 1; x <= v.x + 1; x++) {
				for (int y = v.y - 1; y <= v.y + 1; y++) {
					for (int z = v.z - 1; z <= v.z + 1; z++) {
						if (x == 0 && y == 0 && z == 0)
							continue;
						if (x == p0.x && y == p0.y && z == p0.z) {
							startAtZero = true;
							break outerloop;
						}
					}
				}
			}
		}

		if (startAtZero) {
			if (nthPoint < edgePoints.size())
				return edgePoints.get(nthPoint);
			else
				return edgePoints.get(edgePoints.size() - 1);
		} else {
			if (nthPoint < edgePoints.size())
				return edgePoints.get(edgePoints.size() - nthPoint - 1);
			else
				return edgePoints.get(0);
		}
	}
}
