package org.doube.geometry;

import java.awt.BasicStroke;
import java.awt.Checkbox;
import java.awt.CheckboxGroup;
import java.awt.Choice;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.Label;
import java.awt.Panel;
import java.awt.Point;
import java.awt.Scrollbar;
import java.awt.Shape;
import java.awt.TextField;
import java.awt.event.AdjustmentEvent;
import java.awt.event.AdjustmentListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.awt.event.TextEvent;
import java.awt.event.TextListener;
import java.awt.event.WindowEvent;
import java.awt.geom.AffineTransform;
import java.awt.geom.GeneralPath;
import java.util.Hashtable;

import org.doube.util.UsageReporter;

import ij.IJ;
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
 * Indicator to show directions such as medial, proximal, cranial, north, left
 * 
 * @author Michael Doube
 * @author Wayne Rasband
 * 
 */
@SuppressWarnings("serial")
public class Orienteer extends PlugInFrame implements AdjustmentListener,
		ItemListener, TextListener, MouseWheelListener {

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

	/** Current principal direction choice */
	private int axis0 = 1;

	/** Current secondary direction choice */
	private int axis1 = 0;

	/** Direction labels */
	private String[] directions = { axisLabels[axis0][2], axisLabels[axis0][3],
			axisLabels[axis1][2], axisLabels[axis1][3] };
	private static Orienteer instance;

	/** Current orientation in radians */
	private double theta = 0;

	/** Axis length */
	private int length;

	/** Compass centre coordinates */
	private Point p;

	private Integer activeImpID;
	private Hashtable<Integer, Double> thetaHash = new Hashtable<Integer, Double>();
	private Hashtable<Integer, Integer> lengthHash = new Hashtable<Integer, Integer>();
	private Hashtable<Integer, Point> centreHash = new Hashtable<Integer, Point>();
	private Hashtable<Integer, GeneralPath> pathHash = new Hashtable<Integer, GeneralPath>();
	private Hashtable<Integer, int[]> axisHash = new Hashtable<Integer, int[]>();
	private Hashtable<Integer, boolean[]> reflectHash = new Hashtable<Integer, boolean[]>();
	private Hashtable<Integer, boolean[]> unitHash = new Hashtable<Integer, boolean[]>();

	private Overlay overlay = new Overlay();
	private int fontSize = 12;
	private double scale = 1;
	private GeneralPath path;
	private BasicStroke stroke;

	GridBagLayout gridbag;
	GridBagConstraints c;
	private Scrollbar slider;
	private Panel panel0;
	private Panel panel1;
	private Choice axis0Choice;
	private Choice axis1Choice;
	private Checkbox reflect0;
	private Checkbox reflect1;
	private boolean isReflected0 = false;
	private boolean isReflected1 = false;
	private TextField text;
	private Panel degRadPanel;
	private Checkbox deg;
	private Checkbox rad;

	public Orienteer() {
		super("Orientation");
		if (instance != null) {
			if (!instance.getTitle().equals(getTitle())) {
				Orienteer aa = (Orienteer) instance;
				Prefs.saveLocation(LOC_KEY, aa.getLocation());
				aa.close();
			} else {
				instance.toFront();
				return;
			}
		}
		instance = this;
		IJ.register(Orienteer.class);
		WindowManager.addWindow(this);

		gridbag = new GridBagLayout();
		c = new GridBagConstraints();
		setLayout(gridbag);

		slider = new Scrollbar(Scrollbar.HORIZONTAL,
				(int) (theta * 180 / Math.PI), 1, 0, 360);
		int y = 0;
		c.gridy = y++;
		c.insets = new Insets(2, 10, 0, 10);
		gridbag.setConstraints(slider, c);
		add(slider);
		slider.addAdjustmentListener(this);
		slider.setUnitIncrement(1);
		slider.setFocusable(false); // prevents blinking on Windows
		slider.setPreferredSize(new Dimension(360, 16));
		slider.addMouseWheelListener(this);

		degRadPanel = new Panel();
		Label degRadLabel = new Label("Orientation");
		degRadPanel.add(degRadLabel);
		text = new TextField(IJ.d2s(theta * 180 / Math.PI, 3), 7);
		degRadPanel.add(text);
		text.addTextListener(this);

		CheckboxGroup degRad = new CheckboxGroup();
		deg = new Checkbox("deg", degRad, true);
		rad = new Checkbox("rad", degRad, false);
		degRadPanel.add(deg);
		degRadPanel.add(rad);
		deg.addItemListener(this);
		rad.addItemListener(this);

		panel0 = new Panel();

		Label label0 = new Label("Principal direction");
		panel0.add(label0);

		axis0Choice = new Choice();
		for (int i = 0; i < axisLabels.length; i++)
			axis0Choice.addItem(axisLabels[i][0] + " - " + axisLabels[i][1]);
		axis0Choice.select(axis0);
		axis0Choice.addItemListener(this);
		panel0.add(axis0Choice);

		reflect0 = new Checkbox("Reflect");
		reflect0.setState(isReflected0);
		reflect0.addItemListener(this);
		panel0.add(reflect0);

		panel1 = new Panel();
		Label label1 = new Label("Secondary direction");
		panel1.add(label1);

		axis1Choice = new Choice();
		for (int i = 0; i < axisLabels.length; i++)
			axis1Choice.addItem(axisLabels[i][0] + " - " + axisLabels[i][1]);
		axis1Choice.select(axis1);
		axis1Choice.addItemListener(this);
		panel1.add(axis1Choice);

		reflect1 = new Checkbox("Reflect");
		reflect1.setState(isReflected1);
		reflect1.addItemListener(this);
		panel1.add(reflect1);

		c.gridx = 0;
		c.gridy = y++;
		c.gridwidth = 2;
		c.insets = new Insets(5, 5, 0, 5);
		c.anchor = GridBagConstraints.EAST;
		c.fill = GridBagConstraints.NONE;
		add(degRadPanel, c);

		c.gridx = 0;
		c.gridy = y++;
		c.gridwidth = 2;
		c.insets = new Insets(5, 5, 0, 5);
		c.anchor = GridBagConstraints.EAST;
		c.fill = GridBagConstraints.NONE;
		add(panel0, c);

		c.gridx = 0;
		c.gridy = y++;
		c.gridwidth = 2;
		c.insets = new Insets(5, 5, 0, 5);
		c.anchor = GridBagConstraints.EAST;
		c.fill = GridBagConstraints.NONE;
		add(panel1, c);

		pack();
		Point loc = Prefs.getLocation(LOC_KEY);
		if (loc != null)
			setLocation(loc);
		else
			GUI.center(this);
		if (IJ.isMacOSX())
			setResizable(false);
		setVisible(true);
		ImagePlus imp = WindowManager.getCurrentImage();
		if (imp != null)
			setup(imp);
	}

	private void setup(ImagePlus imp) {
		if (imp == null)
			return;
		if (checkHash(imp)) {
			IJ.log("Image has already been set up");
			return;
		}
		Integer id = new Integer(imp.getID());
		activeImpID = id;
		int w = imp.getWidth();
		int h = imp.getHeight();
		this.theta = 0;
		slider.setValue((int) (theta * 180 / Math.PI));
		this.isReflected0 = false;
		reflect0.setState(isReflected0);
		this.isReflected1 = false;
		reflect1.setState(isReflected1);
		this.length = Math.min(w, h) / 4;
		this.p = new Point(w / 2, h / 2);
		path = new GeneralPath();
		path.moveTo(p.x - length, p.y);
		path.lineTo(p.x + length, p.y);
		path.moveTo(p.x, p.y - length);
		path.lineTo(p.x, p.y + length);
		stroke = new BasicStroke(1, BasicStroke.CAP_ROUND,
				BasicStroke.JOIN_MITER);
		imp.setOverlay(path, Color.BLUE, stroke);
		rotateTo(theta);
		centreHash.put(id, new Point(p));
		thetaHash.put(id, new Double(theta));
		pathHash.put(id, new GeneralPath(path));
		int[] axes = { axis0, axis1 };
		axisHash.put(id, axes.clone());
		lengthHash.put(id, new Integer(this.length));
		boolean[] reflectors = { isReflected0, isReflected1 };
		reflectHash.put(id, reflectors.clone());
		deg.setState(true);
		rad.setState(false);
		boolean[] units = { deg.getState(), rad.getState() };
		unitHash.put(id, units);
		updateTextbox();
		UsageReporter.reportEvent(this).send();
	}

	private void update() {
		ImagePlus imp = WindowManager.getCurrentImage();
		activeImpID = new Integer(imp.getID());
		instance.setTitle("Orientation - " + imp.getTitle());
		if (!checkHash(imp)) {
			setup(imp);
			return;
		}
		this.p = centreHash.get(activeImpID);
		this.theta = thetaHash.get(activeImpID);
		int[] axes = axisHash.get(activeImpID);
		this.axis0 = axes[0];
		this.axis1 = axes[1];
		this.path = pathHash.get(activeImpID);
		this.length = lengthHash.get(activeImpID);
		boolean[] reflectors = reflectHash.get(activeImpID);
		this.isReflected0 = reflectors[0];
		this.isReflected1 = reflectors[1];
		axis0Choice.select(axis0);
		axis1Choice.select(axis1);
		slider.setValue((int) (theta * 180 / Math.PI));
		reflect0.setState(isReflected0);
		reflect1.setState(isReflected1);
		boolean[] units = unitHash.get(activeImpID);
		deg.setState(units[0]);
		rad.setState(units[1]);
		updateTextbox();
		updateDirections();
	}

	/**
	 * Check if an image is already handled by Orientation
	 * 
	 * @param imp
	 * @return true if this image is already handled by Orientation
	 */
	private boolean checkHash(ImagePlus imp) {
		Integer i = new Integer(imp.getID());
		return thetaHash.containsKey(i);
	}

	public static Orienteer getInstance() {
		return instance;
	}

	/**
	 * Gets the labels associated with an ImagePlus. Returns null if the image
	 * isn't being tracked by Orienteer.
	 * 
	 * @param imp
	 * @return an array of axis labels, with the principal direction in the
	 *         zeroth position, the principal tail in the 1st position, the
	 *         secondary head in the 2nd position and the secondary tail in the
	 *         3rd position.
	 */
	public String[] getDirections(ImagePlus imp) {
		if (!checkHash(imp))
			return null;
		Integer id = new Integer(imp.getID());
		int[] axes = axisHash.get(id);
		boolean[] ref = reflectHash.get(id);
		String[] dirs = new String[4];

		if (!ref[0]) {
			dirs[0] = axisLabels[axes[0]][2];
			dirs[1] = axisLabels[axes[0]][3];
		} else {
			dirs[0] = axisLabels[axes[0]][3];
			dirs[1] = axisLabels[axes[0]][2];
		}
		if (!ref[1]) {
			dirs[2] = axisLabels[axes[1]][2];
			dirs[3] = axisLabels[axes[1]][3];
		} else {
			dirs[2] = axisLabels[axes[1]][3];
			dirs[3] = axisLabels[axes[1]][2];
		}
		return dirs;
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
		ImagePlus imp = WindowManager.getCurrentImage();
		return getOrientation(imp, direction);
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
	 * Get the orientation of the principal direction
	 * 
	 * @param imp
	 * @return Orientation in radians clockwise from 12 o'clock
	 * @throws IllegalArgumentException
	 *             if imp has not been activated by Orientation
	 */
	public double getOrientation(ImagePlus imp) {
		Integer id = new Integer(imp.getID());
		Double o = thetaHash.get(id);
		if (o == null) {
			throw new IllegalArgumentException();
		} else {
			return o.doubleValue();
		}
	}

	public double getOrientation(ImagePlus imp, String direction) {
		double orientation = getOrientation(imp);
		String[] dir = getDirections(imp);

		int quadrant = 0;
		for (int i = 0; i < 4; i++) {
			if (dir[i].equals(direction)) {
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
	 * Given a set of (x,y) coordinates, find the caliper diameters across the
	 * axes
	 * 
	 * @param points
	 *            in double[n][2] format
	 * @return caliper diameters across the principal and secondary axes (zeroth
	 *         and first elements respectively)
	 */
	public double[] getDiameters(double[][] points) {
		double xMin = Double.POSITIVE_INFINITY;
		double xMax = Double.NEGATIVE_INFINITY;
		double yMin = Double.POSITIVE_INFINITY;
		double yMax = Double.NEGATIVE_INFINITY;
		final int nPoints = points.length;
		for (int i = 0; i < nPoints; i++) {
			double x = points[i][0];
			double y = points[i][1];
			double xr = x * Math.cos(-theta) - y * Math.sin(-theta);
			double yr = x * Math.sin(-theta) + y * Math.cos(-theta);
			xMin = Math.min(xMin, xr);
			xMax = Math.max(xMax, xr);
			yMin = Math.min(yMin, yr);
			yMax = Math.max(yMax, yr);
		}
		double[] result = { yMax - yMin, xMax - xMin };
		return result;
	}

	/**
	 * Rotate the direction indicator by a given angle
	 * 
	 * @param deltaTheta
	 *            number of radians to rotate by (+ve is clockwise, -ve is
	 *            anti-clockwise)
	 */
	public void rotate(double deltaTheta) {
		ImagePlus imp = WindowManager.getCurrentImage();
		this.activeImpID = new Integer(imp.getID());
		AffineTransform transform = new AffineTransform();
		transform.rotate(deltaTheta, p.x, p.y);
		this.theta += deltaTheta;
		path.transform(transform);
		overlay.clear();
		addPath(path, Color.BLUE, stroke);
		addLabels();
		imp.setOverlay(overlay);
		thetaHash.put(activeImpID, new Double(theta));
		pathHash.put(activeImpID, new GeneralPath(path));
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
		thetaHash.put(activeImpID, newTheta);
	}

	private void addLabels() {
		ImagePlus imp = WindowManager.getImage(activeImpID);
		scale = imp.getCanvas().getMagnification();
		Font font = new Font("SansSerif", Font.PLAIN, (int) (fontSize / scale));
		final double lsinTheta = (length + (fontSize / scale))
				* Math.sin(theta);
		final double lcosTheta = (length + (fontSize / scale))
				* Math.cos(theta);
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
		directions = getDirections(WindowManager.getImage(activeImpID));
		rotate(0);
	}

	private void updateTextbox() {
		if (deg.getState())
			text.setText(IJ.d2s(this.theta * 180 / Math.PI, 3));
		else
			text.setText(IJ.d2s(this.theta, 5));
	}

	void addPath(Shape shape, Color color, BasicStroke stroke) {
		Roi roi = new ShapeRoi(shape);
		roi.setStrokeColor(color);
		roi.setStroke(stroke);
		roi.setStrokeWidth(roi.getStrokeWidth() / (float) scale);
		overlay.add(roi);
	}

	void addString(String text, int x, int y, Color color, Font font) {
		TextRoi roi = new TextRoi(x, y, text, font);
		roi.setLocation(x - text.length() * (int) (fontSize / scale) / 4, y
				- (int) (fontSize / scale) / 2);
		roi.setStrokeColor(color);
		overlay.add(roi);
	}

	public void close() {
		super.close();
		instance = null;
		for (Integer i : thetaHash.keySet())
			WindowManager.getImage(i.intValue()).setOverlay(null);
	}

	public void windowActivated(WindowEvent e) {
		super.windowActivated(e);
		update();
		WindowManager.setWindow(this);
	}

	public void windowClosing(WindowEvent e) {
		close();
		Prefs.saveLocation(LOC_KEY, getLocation());
	}

	public void adjustmentValueChanged(AdjustmentEvent e) {
		if (e.getSource().equals(slider)) {
			rotateTo(slider.getValue() * Math.PI / 180);
			updateTextbox();
		}
	}

	public void itemStateChanged(ItemEvent e) {
		Object source = e.getSource();
		if (source.equals(axis0Choice)) {
			int i = axis0Choice.getSelectedIndex();
			if (i == axis1Choice.getSelectedIndex()) {
				axis0Choice.select(axis0);
				IJ.error("Both axes cannot indicate the same direction");
				return;
			}
			axis0 = i;
			int[] axes = { axis0, axis1 };
			axisHash.put(activeImpID, axes.clone());
			updateDirections();
		} else if (source.equals(axis1Choice)) {
			int i = axis1Choice.getSelectedIndex();
			if (i == axis0Choice.getSelectedIndex()) {
				axis1Choice.select(axis1);
				IJ.error("Both axes cannot indicate the same direction");
				return;
			}
			axis1 = i;
			int[] axes = { axis0, axis1 };
			axisHash.put(activeImpID, axes.clone());
			updateDirections();
		} else if (source.equals(reflect0)) {
			isReflected0 = reflect0.getState();
			boolean[] reflectors = { isReflected0, isReflected1 };
			reflectHash.put(activeImpID, reflectors.clone());
			updateDirections();
		} else if (source.equals(reflect1)) {
			isReflected1 = reflect1.getState();
			boolean[] reflectors = { isReflected0, isReflected1 };
			reflectHash.put(activeImpID, reflectors.clone());
			updateDirections();
		} else if (source.equals(deg) || source.equals(rad)) {
			boolean[] units = { deg.getState(), rad.getState() };
			unitHash.put(activeImpID, units);
			updateTextbox();
		}
	}

	public void textValueChanged(TextEvent e) {
		TextField field = (TextField) e.getSource();
		double value = Double.parseDouble(field.getText());
		if (deg.getState())
			value *= Math.PI / 180;
		value = value % (2 * Math.PI);
		if (value < 0)
			value += 2 * Math.PI;
		slider.setValue((int) Math.round((value * 180 / Math.PI)));
		rotateTo(value);
	}

	public void mouseWheelMoved(MouseWheelEvent e) {
		final int oldPos = slider.getValue();
		int newPos = oldPos + e.getWheelRotation();
		if (newPos < 0)
			newPos += 360;
		else if (newPos >= 360)
			newPos -= 360;
		rotateTo(newPos * Math.PI / 180);
		updateTextbox();
		slider.setValue(newPos);
	}
}
