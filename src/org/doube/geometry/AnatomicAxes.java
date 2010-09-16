package org.doube.geometry;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Frame;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.Scrollbar;
import java.awt.Shape;
import java.awt.event.AdjustmentEvent;
import java.awt.event.AdjustmentListener;
import java.awt.event.WindowEvent;
import java.awt.geom.AffineTransform;
import java.awt.geom.GeneralPath;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.GUI;
import ij.gui.Overlay;
import ij.gui.Roi;
import ij.gui.ShapeRoi;
import ij.gui.TextRoi;
import ij.plugin.frame.PlugInFrame;

/**
 * Indicator to show anatomic directions such as medial, proximal, cranial
 * 
 * Instantiate 2 line ROI's, lock their lengths, lock their midpoint to always
 * be the same position.
 * 
 * @author Michael Doube
 * @author Wayne Rasband
 * 
 */
@SuppressWarnings("serial")
public class AnatomicAxes extends PlugInFrame implements AdjustmentListener {

	public static final String LOC_KEY = "aa.loc";
	private static final String[][] axisLabelsFull = { { "medial", "lateral" },
			{ "cranial", "caudal" }, { "rostral", "caudal" },
			{ "dorsal", "ventral" }, { "anterior", "posterior" },
			{ "superior", "inferior" }, { "proximal", "distal" },
			{ "dorsal", "palmar" }, { "dorsal", "plantar" },
			{ "dorsal", "volar" }, { "axial", "abaxial" },
			{ "north", "south" }, { "east", "west" }, { "up", "down" },
			{ "right", "left" } };

	private static final String[][] axisLabels = { { "M", "L" },
			{ "Cr", "Ca" }, { "Ro", "Ca" }, { "D", "V" }, { "A", "P" },
			{ "Sup", "Inf" }, { "Pr", "Di" }, { "Do", "Pa" }, { "Do", "Pl" },
			{ "Do", "Vo" }, { "Ax", "Ab" }, { "N", "S" }, { "E", "W" },
			{ "Up", "D" }, { "R", "L" } };

	private String[] directions = { axisLabels[1][0], axisLabels[1][1],
			axisLabels[0][0], axisLabels[0][1] };

	private static Frame instance;
	ImageJ ij;
	GridBagLayout gridbag;
	GridBagConstraints c;

	/** Common midpoint */
	private int x, y;

	/** Current orientation in radians */
	private double theta = 0;

	/** Axis length */
	private int length;
	private ImagePlus imp;
	private Point p;
	private GeneralPath path;
	private BasicStroke stroke;
	private Scrollbar slider;

	private Overlay overlay = new Overlay();
	private int fontSize = 12;

	public AnatomicAxes() {
		super("Orientation");
	}

	public void run(String arg) {
		if (instance != null) {
			if (!instance.getTitle().equals(getTitle())) {
				AnatomicAxes aa = (AnatomicAxes) instance;
				Prefs.saveLocation(LOC_KEY, aa.getLocation());
				aa.close();
			} else {
				instance.toFront();
				return;
			}
		}
		instance = this;
		IJ.register(AnatomicAxes.class);
		WindowManager.addWindow(this);

		ij = IJ.getInstance();
		gridbag = new GridBagLayout();
		c = new GridBagConstraints();
		setLayout(gridbag);

		slider = new Scrollbar(Scrollbar.HORIZONTAL,
				(int) (theta * 180 / Math.PI), 1, 0, 360);
		c.gridy = y++;
		c.insets = new Insets(2, 10, 0, 10);
		gridbag.setConstraints(slider, c);
		add(slider);
		slider.addAdjustmentListener(this);
		slider.addKeyListener(ij);
		slider.setUnitIncrement(1);
		slider.setFocusable(false); // prevents blinking on Windows
		slider.setPreferredSize(new Dimension(180, 16));
		pack();
		Point loc = Prefs.getLocation(LOC_KEY);
		if (loc != null)
			setLocation(loc);
		else
			GUI.center(this);
		if (IJ.isMacOSX())
			setResizable(false);
		setVisible(true);

		// show the crosshairs
		this.imp = WindowManager.getCurrentImage();
		if (imp == null)
			return;
		int w = imp.getWidth();
		int h = imp.getHeight();
		this.length = Math.min(w, h) / 4;
		this.x = w / 2;
		this.y = h / 2;
		this.p = new Point(x, y);
		path = new GeneralPath();
		drawCross(p, path);
		stroke = new BasicStroke(1, BasicStroke.CAP_ROUND,
				BasicStroke.JOIN_MITER);
		this.imp.setOverlay(path, Color.BLUE, stroke);
		rotateTo(0);
	}

	public void reflect() {
		// TODO method stub
		// get whichever line is selected
		// swap the labels
	}

	/**
	 * Retrieve the orientation of the named axis. If the axis name is not
	 * found, the principal axis orientation is returned.
	 * 
	 * @param direction
	 *            Label you wish to match
	 * @return orientation of the axis in radians clockwise from 12 o'clock
	 */
	public double getOrientation(String direction) {
		double orientation = this.theta;
		int quadrant = 0;
		for (int i = 0; i < 4; i++) {
			if (directions[i].equals(direction)) {
				quadrant = i;
				break;
			}
		}

		switch (quadrant) {
		case 0:
			return orientation;
		case 1:
			orientation += Math.PI; break;
		case 2:
			orientation += Math.PI / 2; break;
		case 3:
			orientation += 3 * Math.PI / 2; break;
		}

		if (orientation > 2 * Math.PI)
			return orientation - 2 * Math.PI;
		else
			return orientation;
	}

	public void rotate(double deltaTheta) {
		AffineTransform transform = new AffineTransform();
		transform.rotate(deltaTheta, p.x, p.y);
		this.theta += deltaTheta;
		path.transform(transform);
		overlay.clear();
		addPath(path, Color.BLUE, stroke);
		addLabels();
		imp.setOverlay(overlay);
	}

	public void rotateTo(double newTheta) {
		rotate(newTheta - this.theta);
		this.theta = newTheta;
	}

	private void addLabels() {
		Font font = new Font("SansSerif", Font.PLAIN, fontSize);
		final double lsinTheta = (length + fontSize) * Math.sin(theta);
		final double lcosTheta = (length + fontSize) * Math.cos(theta);
		addString(directions[0], (int) (p.x + lsinTheta),
				(int) (p.y - lcosTheta), Color.RED, font);
		addString(directions[1], (int) (p.x - lsinTheta),
				(int) (p.y + lcosTheta), Color.BLUE, font);
		addString(directions[2], (int) (p.x + lcosTheta),
				(int) (p.y + lsinTheta), Color.BLUE, font);
		addString(directions[3], (int) (p.x - lcosTheta),
				(int) (p.y - lsinTheta), Color.BLUE, font);
	}

	void addPath(Shape shape, Color color, BasicStroke stroke) {
		Roi roi = new ShapeRoi(shape);
		roi.setStrokeColor(color);
		roi.setStroke(stroke);
		overlay.add(roi);
	}

	void addString(String text, int x, int y, Color color, Font font) {
		TextRoi roi = new TextRoi(x, y, text, font);
		roi.setLocation(x - text.length() * fontSize / 4, y - fontSize / 2);
		roi.setStrokeColor(color);
		overlay.add(roi);
	}

	/** draws the crosses in the images */
	private void drawCross(Point p, GeneralPath path) {
		float x = p.x;
		float y = p.y;
		path.moveTo(x - length, y);
		path.lineTo(x + length, y);
		path.moveTo(x, y - length);
		path.lineTo(x, y + length);
	}

	public void close() {
		super.close();
		instance = null;
		imp.setOverlay(null);
	}

	public void windowClosing(WindowEvent e) {
		close();
		Prefs.saveLocation(LOC_KEY, getLocation());
	}

	@Override
	public void adjustmentValueChanged(AdjustmentEvent e) {
		if (e.getSource().equals(slider)) {
			rotateTo(slider.getValue() * Math.PI / 180);
		}
	}
}
