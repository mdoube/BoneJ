package org.doube.geometry;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Frame;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.Point;
import java.awt.Scrollbar;
import java.awt.event.AdjustmentEvent;
import java.awt.event.AdjustmentListener;
import java.awt.geom.AffineTransform;
import java.awt.geom.GeneralPath;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.GUI;
import ij.gui.ImageCanvas;
import ij.plugin.frame.ContrastAdjuster;
import ij.plugin.frame.PlugInFrame;

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
public class AnatomicAxes extends PlugInFrame implements AdjustmentListener {

	public static final String LOC_KEY = "aa.loc";
	private static final String[][] axisLabels = { { "medial", "lateral" },
			{ "cranial", "caudal" }, { "rostral", "caudal" },
			{ "dorsal", "ventral" }, { "anterior", "posterior" },
			{ "superior", "inferior" }, { "proximal", "distal" },
			{ "dorsal", "palmar" }, { "dorsal", "plantar" },
			{ "dorsal", "volar" }, { "axial", "abaxial" } };

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
	private ImageCanvas canvas;
	private Point p;
	private GeneralPath path;
	private BasicStroke stroke;
	private Scrollbar slider;

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
		IJ.register(ContrastAdjuster.class);
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
		this.canvas = imp.getCanvas();
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
		this.imp.setOverlay(path, Color.RED, stroke);
		canvas.setCustomRoi(true);
	}

	public void reflect() {
		// TODO method stub
		// get whichever line is selected
		// swap the labels
	}

	public double getOrientation(String direction) {
		double principal = this.theta;
		return principal;
	}

	public void rotate(double deltaTheta) {
		AffineTransform transform = new AffineTransform();
		transform.rotate(deltaTheta, p.x, p.y);
		this.theta += deltaTheta;
		path.transform(transform);
		this.imp.setOverlay(path, Color.RED, stroke);
	}

	public void rotateTo(double newTheta) {
		rotate(newTheta - this.theta);
		this.theta = newTheta;
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

	@Override
	public void adjustmentValueChanged(AdjustmentEvent e) {
		if (e.getSource().equals(slider))
			rotateTo(slider.getValue() * Math.PI / 180);
	}
}
