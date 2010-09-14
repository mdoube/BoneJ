package org.doube.geometry;

import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics;

import ij.ImagePlus;
import ij.gui.Line;
import ij.gui.Overlay;
import ij.gui.Roi;
import ij.gui.TextRoi;

/**
 * Indicator to show anatomic directions such as medial, proximal, cranial
 * 
 * Instantiate 2 line ROI's, lock their lengths, lock their midpoint to always
 * be the same position.
 * 
 * @author Michael Doube
 * 
 */
@SuppressWarnings("serial")
public class AnatomicAxes {

	private static final String[][] axisLabels = { { "medial", "lateral" },
			{ "cranial", "caudal" }, { "rostral", "caudal" },
			{ "dorsal", "ventral" }, { "anterior", "posterior" },
			{ "superior", "inferior" }, { "proximal", "distal" },
			{ "dorsal", "palmar" }, { "dorsal", "plantar" },
			{ "dorsal", "volar" }, { "axial", "abaxial" } };

	/** Common midpoint */
	private int x, y;

	/** Axis length */
	private int length;
	private ImagePlus imp;
	private Line line1, line2;

	/** TextRoi labels for axes */
	private TextRoi tr1, tr2, tr3, tr4;

	private Graphics g;
	private Overlay overlay;

	public AnatomicAxes(ImagePlus imp) {
		this.imp = imp;
		int w = imp.getWidth();
		int h = imp.getHeight();
		this.length = Math.min(w, h) / 4;
		this.x = w / 2;
		this.y = h / 2;
		this.line1 = new Line(x, y - length, x, y + length);
		this.line2 = new Line(x - length, y, x + length, y);
		TextRoi.setFont("SansSerif", 10, Font.PLAIN, true);
		this.tr1 = new TextRoi(x, y - length - 10, axisLabels[0][0]);
		this.tr2 = new TextRoi(x, y + length + 10, axisLabels[0][1]);
		this.tr3 = new TextRoi(x - length - 10, y, axisLabels[1][0]);
		this.tr4 = new TextRoi(x + length + 10, y, axisLabels[1][1]);
		this.overlay = new Overlay();
		overlay.add(line1);
		overlay.add(line2);
		overlay.add(tr1);
		overlay.add(tr2);
		overlay.add(tr3);
		overlay.add(tr4);
		for (Roi r : overlay.toArray()){
			r.setImage(this.imp);
		}
		overlay.setStrokeColor(new Color(0, 0, 255));
		this.imp.setOverlay(overlay);
		this.g = this.imp.getCanvas().getGraphics();
		create();
	}

	public void create() {
		for (Roi r : overlay.toArray())
			r.draw(this.g);
	}

	public void reflect() {
		// TODO method stub
		// get whichever line is selected
		// swap the labels
	}

}
