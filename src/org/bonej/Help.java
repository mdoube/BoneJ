package org.bonej;

import java.awt.BorderLayout;
import java.io.IOException;

import javax.swing.JEditorPane;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.WindowConstants;
import javax.swing.event.HyperlinkEvent;
import javax.swing.event.HyperlinkListener;

import org.doube.util.UsageReporter;

import ij.plugin.BrowserLauncher;
import ij.plugin.PlugIn;

public class Help implements PlugIn {

	/**
	 * branch
	 */
	public static final String branch = "-testing";

	/**
	 * BoneJ version
	 */
	public static final String bonejVersion = "1.4.1" + branch;

	@Override
	public void run(final String arg) {
		if (arg.equals("about")) {
			showAbout();
			// testing
			UsageReporter.reportEvent(this).send();
			return;
		}
	}

	private void showAbout() {
		final JEditorPane htmlPane = new JEditorPane("text/html",
				"<html>\n" + "  <body>\n" + "<p><b>BoneJ version " + bonejVersion + "</b>	</p>"
						+ "<p>BoneJ is an ImageJ plugin designed for (but not limited to) bone image analysis.</p>"
						+ "<p>User and developer documentation can be found at <a href=http://bonej.org/>bonej.org</a></p>"
						+ "\n" + "  </body>\n" + "</html>");
		htmlPane.setEditable(false);
		htmlPane.setOpaque(false);
		htmlPane.addHyperlinkListener(new HyperlinkListener() {
			@Override
			public void hyperlinkUpdate(final HyperlinkEvent e) {
				if (HyperlinkEvent.EventType.ACTIVATED.equals(e.getEventType()))
					try {
						BrowserLauncher.openURL(e.getURL().toString());
					} catch (final IOException exception) {
						// ignore
					}
			}
		});

		final JPanel panel = new JPanel(new BorderLayout());
		panel.add(htmlPane, BorderLayout.CENTER);

		final JFrame frame = new JFrame("About BoneJ...");
		frame.setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
		frame.getContentPane().add(panel, BorderLayout.CENTER);
		frame.pack();
		frame.setVisible(true);
	}

}
