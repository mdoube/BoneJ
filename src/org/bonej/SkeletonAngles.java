package org.bonej;

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

	private static final int VERTEX_TO_VERTEX = -1;

	public void run(String arg) {
		ImagePlus imp = WindowManager.getCurrentImage();
		if (imp == null) {
			IJ.noImage();
			return;
		}

		GenericDialog gd = new GenericDialog("Options");
		String[] options = { "Opposite edges", "Vertex and nth edge pixel" };
		gd.addChoice("Measure angle between...", options, options[0]);
		gd.addNumericField("nth pixel", 0, 0);

		gd.showDialog();

		if (gd.wasCanceled())
			return;

		int choice = gd.getNextChoiceIndex();

		int nthPixel = (int) Math.floor(gd.getNextNumber());

		if (choice == 0) {
			nthPixel = VERTEX_TO_VERTEX;
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

					double theta01 = vertexAngle(vertex, edge0, edge1, nthPixel);
					double theta02 = vertexAngle(vertex, edge0, edge2, nthPixel);
					double theta12 = vertexAngle(vertex, edge1, edge2, nthPixel);

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
		UsageReporter.reportEvent(this).send();
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

		// double x0 = centroid0[0] - centroidv[0];
		// double y0 = centroid0[1] - centroidv[1];
		// double z0 = centroid0[2] - centroidv[2];
		// double x1 = centroid1[0] - centroidv[0];
		// double y1 = centroid1[1] - centroidv[1];
		// double z1 = centroid1[2] - centroidv[2];
		//
		// double dot = x0 * x1 + y0 * y1 + z0 * z1;
		// double d0 = Trig.distance3D(centroidv, centroid0);
		// double d1 = Trig.distance3D(centroidv, centroid1);
		//
		// double cosTheta = dot / (d0 * d1);
		//
		// return Math.acos(cosTheta);
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
