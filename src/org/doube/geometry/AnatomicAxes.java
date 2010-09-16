package org.doube.geometry;

import java.awt.BasicStroke;
import java.awt.Checkbox;
import java.awt.Choice;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Frame;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.Panel;
import java.awt.Point;
import java.awt.Scrollbar;
import java.awt.Shape;
import java.awt.event.AdjustmentEvent;
import java.awt.event.AdjustmentListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
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
public class AnatomicAxes extends PlugInFrame implements AdjustmentListener,
		ItemListener, KeyListener {

	public static final String LOC_KEY = "aa.loc";
	private static final String[][] axisLabels = {
			{ "medial", "lateral", "M", "L" },
			{ "cranial", "caudal", "Cr", "Ca" },
			{ "rostral", "caudal", "Ro", "Ca" },
			{ "dorsal", "ventral", "D", "V" },
			{ "anterior", "posterior", "A", "P" },
			{ "superior", "inferior", "Sup", "Inf" },
			{ "proximal", "distal", "Pr", "Di" },
			{ "dorsal", "palmar", "Do", "Pa" },
			{ "dorsal", "plantar", "Do", "Pl" },
			{ "dorsal", "volar", "Do", "Vo" },
			{ "axial", "abaxial", "Ax", "Ab" }, { "north", "south", "N", "S" },
			{ "east", "west", "E", "W" }, { "up", "down", "Up", "D" },
			{ "right", "left", "R", "L" } };

	private int axis0 = 1, axis1 = 0;
	private String[] directions = { axisLabels[axis0][2], axisLabels[axis0][3],
			axisLabels[axis1][2], axisLabels[axis1][3] };
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
	private Panel panel;
	private Choice axis0Choice;
	private Choice axis1Choice;
	private Checkbox reflect;
	private boolean isReflected = false;

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
		slider.addKeyListener(this);
		slider.setUnitIncrement(1);
		slider.setFocusable(false); // prevents blinking on Windows
		slider.setPreferredSize(new Dimension(360, 16));

		panel = new Panel();
		axis0Choice = new Choice();
		for (int i = 0; i < axisLabels.length; i++)
			axis0Choice.addItem(axisLabels[i][0] + " - " + axisLabels[i][1]);
		axis0Choice.select("cranial - caudal");
		axis0Choice.addItemListener(this);
		// methodChoice.addKeyListener(ij);
		panel.add(axis0Choice);

		reflect = new Checkbox("Reflect");
		reflect.setState(isReflected);
		reflect.addItemListener(this);
		c.gridx = 0;
		c.gridy = y++;
		c.gridwidth = 2;
		c.insets = new Insets(5, 35, 0, 5);
		add(reflect, c);

		axis1Choice = new Choice();
		for (int i = 0; i < axisLabels.length; i++)
			axis1Choice.addItem(axisLabels[i][0] + " - " + axisLabels[i][1]);
		axis1Choice.select("medial - lateral");
		axis1Choice.addItemListener(this);
		// modeChoice.addKeyListener(ij);
		panel.add(axis1Choice);
		c.gridx = 0;
		c.gridy = y++;
		c.gridwidth = 2;
		c.insets = new Insets(5, 5, 0, 5);
		c.anchor = GridBagConstraints.CENTER;
		c.fill = GridBagConstraints.NONE;
		add(panel, c);

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
			orientation += Math.PI;
			break;
		case 2:
			orientation += Math.PI / 2;
			break;
		case 3:
			orientation += 3 * Math.PI / 2;
			break;
		}

		if (orientation > 2 * Math.PI)
			return orientation - 2 * Math.PI;
		else
			return orientation;
	}

	/**
	 * Overloaded version of getOrientation that takes no argument
	 * 
	 * @return orientation of the principal direction in radians clockwise from
	 *         12 o'clock
	 */
	public double getOrientation() {
		double orientation = this.theta;
		return orientation;
	}

	/**
	 * Rotate the direction indicator by a given angle
	 * 
	 * @param deltaTheta
	 *            number of radians to rotate by (+ve is clockwise, -ve is
	 *            anti-clockwise)
	 */
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

	/**
	 * Rotate the principal direction to a new angle
	 * 
	 * @param newTheta
	 *            desired orientation in radians clockwise from 12 o'clock of
	 *            the principal direction
	 */
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

	private void updateDirections() {
		if (!isReflected) {
			directions[0] = axisLabels[axis0][2];
			directions[1] = axisLabels[axis0][3];
		} else {
			directions[0] = axisLabels[axis0][3];
			directions[1] = axisLabels[axis0][2];
		}
		directions[2] = axisLabels[axis1][2];
		directions[3] = axisLabels[axis1][3];
		rotate(0);
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

	public void windowActivated(WindowEvent e) {
		super.windowActivated(e);
		setup();
		WindowManager.setWindow(this);
	}

	private void setup() {
		// TODO Auto-generated method stub
		// Handle situation where instance of orientation exists and a new
		// image is opened
	}

	public void windowClosing(WindowEvent e) {
		close();
		Prefs.saveLocation(LOC_KEY, getLocation());
	}

	public void adjustmentValueChanged(AdjustmentEvent e) {
		if (e.getSource().equals(slider)) {
			rotateTo(slider.getValue() * Math.PI / 180);
		}
	}

	@Override
	public void itemStateChanged(ItemEvent e) {
		if (e.getSource().equals(axis0Choice)) {
			int i = axis0Choice.getSelectedIndex();
			if (i == axis1Choice.getSelectedIndex()) {
				axis0Choice.select(axis0);
				IJ.error("Both axes cannot indicate the same direction");
				return;
			}
			axis0 = i;
			updateDirections();
		} else if (e.getSource().equals(axis1Choice)) {
			int i = axis1Choice.getSelectedIndex();
			if (i == axis0Choice.getSelectedIndex()) {
				axis1Choice.select(axis1);
				IJ.error("Both axes cannot indicate the same direction");
				return;
			}
			axis1 = i;
			updateDirections();
		} else if (e.getSource().equals(reflect)) {
			isReflected = reflect.getState();
			updateDirections();
		}

	}

	@Override
	public void keyPressed(KeyEvent e) {
		// TODO Auto-generated method stub

	}

	@Override
	public void keyReleased(KeyEvent e) {
		// TODO Auto-generated method stub

	}

	@Override
	public void keyTyped(KeyEvent e) {
		// TODO Auto-generated method stub

	}
}
