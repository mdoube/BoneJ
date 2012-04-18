package org.bonej;

import java.util.ArrayList;

import org.doube.geometry.Trig;
import org.doube.skeleton.AnalyzeSkeleton;
import org.doube.skeleton.Edge;
import org.doube.skeleton.Graph;
import org.doube.skeleton.Point;
import org.doube.skeleton.Vertex;
import org.doube.util.ResultInserter;

import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.plugin.PlugIn;

public class SkeletonAngles implements PlugIn {

	public void run(String arg) {
		ImagePlus imp = WindowManager.getCurrentImage();
		if (imp == null) {
			IJ.noImage();
			return;
		}
		
		AnalyzeSkeleton skeletonAnalyzer = new AnalyzeSkeleton();
		skeletonAnalyzer.setup("", imp);
		skeletonAnalyzer.run();
		
		ResultInserter ri = ResultInserter.getInstance();
		
		Graph[] graphs = skeletonAnalyzer.getGraphs();
		int g = 0; 
		int v = 0;
		for (Graph graph : graphs) {
			ArrayList<Vertex> vertices = graph.getVertices();
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

					double theta01 = vertexAngle(vertex, edge0, edge1);
					double theta02 = vertexAngle(vertex, edge0, edge2);
					double theta12 = vertexAngle(vertex, edge1, edge2);
					
					ri.setResultInRow(imp, "Skeleton", g);
					ri.setResultInRow(imp, "Vertex", v);
					ri.setResultInRow(imp, "Theta 0", theta01);
					ri.setResultInRow(imp, "Theta 1", theta02);
					ri.setResultInRow(imp, "Theta 2", theta12);	
				}
				v++;
			}
			g++;
		}
		ri.updateTable();
	}

	private double vertexAngle(Vertex vertex, Edge edge0, Edge edge1) {
		Vertex v0 = edge0.getOppositeVertex(vertex);
		Vertex v1 = edge1.getOppositeVertex(vertex);
		ArrayList<Point> pointsv = vertex.getPoints();
		ArrayList<Point> points0 = v0.getPoints();
		ArrayList<Point> points1 = v1.getPoints();
		double[] centroidv = calculateCentroid(pointsv);
		double[] centroid0 = calculateCentroid(points0);
		double[] centroid1 = calculateCentroid(points1);

		double x0 = centroid0[0] - centroidv[0];
		double y0 = centroid0[1] - centroidv[1];
		double z0 = centroid0[2] - centroidv[2];
		double x1 = centroid1[0] - centroidv[0];
		double y1 = centroid1[1] - centroidv[1];
		double z1 = centroid1[2] - centroidv[2];

		double dot = x0 * x1 + y0 * y1 + z0 * z1;
		double d0 = Trig.distance3D(centroidv, centroid0);
		double d1 = Trig.distance3D(centroidv, centroid1);

		double cosTheta = dot / (d0 * d1);
		
		return Math.acos(cosTheta);
	}

	private double[] calculateCentroid(ArrayList<Point> points) {
		double xsum = 0;
		double ysum = 0;
		double zsum = 0;
		double n = points.size();

		for (Point p : points) {
			xsum += p.x;
			ysum += p.y;
			zsum += p.z;
		}
		double[] centroid = { xsum / n, ysum / n, zsum / n };
		return centroid;
	}
}
