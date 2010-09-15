package org.doube.geometry;

import java.awt.AWTEvent;
import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Component;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.Frame;
import java.awt.Graphics;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.Label;
import java.awt.Panel;
import java.awt.Point;
import java.awt.Scrollbar;
import java.awt.event.AdjustmentEvent;
import java.awt.event.AdjustmentListener;
import java.awt.geom.AffineTransform;
import java.awt.geom.GeneralPath;
import java.util.Vector;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.DialogListener;
import ij.gui.GUI;
import ij.gui.GenericDialog;
import ij.gui.ImageCanvas;
import ij.gui.Line;
import ij.gui.NonBlockingGenericDialog;
import ij.gui.Overlay;
import ij.gui.Roi;
import ij.gui.TextRoi;
import ij.plugin.frame.ContrastAdjuster;
import ij.plugin.frame.PlugInFrame;
//import ij.plugin.frame.TrimmedLabel;
import ij.process.ColorProcessor;
import ij.process.ImageProcessor;

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
	private double theta;

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
		this.imp = IJ.getImage();
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
		NonBlockingGenericDialog nbgd = new NonBlockingGenericDialog(
				"Orientation");
		nbgd.addSlider("Orientation", 0, 360, 0);
		for (Component c : nbgd.getComponents()) {
			if (c instanceof Panel) {
				Panel p = (Panel) c;
				for (Component pc : p.getComponents()) {
					if (pc instanceof Scrollbar) {
						slider = (Scrollbar) pc;
					}
				}
			}
		}
		slider.addAdjustmentListener(this);
		nbgd.showDialog();
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

		slider = new Scrollbar(Scrollbar.HORIZONTAL, 0, 1, 0, 360);
		c.gridy = y++;
		c.insets = new Insets(2, 10, 0, 10);
		gridbag.setConstraints(slider, c);
		add(slider);
		slider.addAdjustmentListener(this);
		slider.addKeyListener(ij);
		slider.setUnitIncrement(1);
		slider.setFocusable(false); // prevents blinking on Windows
		pack();
		Point loc = Prefs.getLocation(LOC_KEY);
		if (loc!=null)
			setLocation(loc);
		else
			GUI.center(this);
		if (IJ.isMacOSX()) setResizable(false);
		setVisible(true);
//		addLabel("Minimum", null);
	}

//	private void addLabel(String text, Label label2) {
//		if (label2 == null && IJ.isMacOSX())
//			text += "    ";
//		panel = new Panel();
//		c.gridy = y++;
//		int bottomInset = IJ.isMacOSX() ? 4 : 0;
//		c.insets = new Insets(0, 10, bottomInset, 0);
//		gridbag.setConstraints(panel, c);
//		panel.setLayout(new FlowLayout(label2 == null ? FlowLayout.CENTER
//				: FlowLayout.LEFT, 0, 0));
//		Label label = new TrimmedLabel(text);
//		label.setFont(sanFont);
//		panel.add(label);
//		if (label2 != null) {
//			label2.setFont(monoFont);
//			label2.setAlignment(Label.LEFT);
//			panel.add(label2);
//		}
//		add(panel);
//	}

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
