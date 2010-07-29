import ij.plugin.PlugIn;
import ij.plugin.BrowserLauncher;

import java.awt.BorderLayout;

import java.io.IOException;

import javax.swing.JEditorPane;
import javax.swing.JFrame;
import javax.swing.JPanel;

import javax.swing.event.HyperlinkEvent;
import javax.swing.event.HyperlinkListener;

public class Help_ implements PlugIn {

	/**
	 * BoneJ version
	 */
	public static final String bonejVersion = "1.1.10";

	public void run(String arg) {
		if (arg.equals("about")) {
			showAbout();
			return;
		}
	}

	private void showAbout() {
		JEditorPane htmlPane = new JEditorPane(
				"text/html",
				"<html>\n"
						+ "  <body>\n"
						+ "<p><b>BoneJ version "
						+ bonejVersion
						+ "</b>	</p>"
						+ "<p>BoneJ is an ImageJ plugin designed for (but not limited to) bone image analysis.</p>"
						+ "<p>User and developer documentation can be found at <a href=http://bonej.org/>bonej.org</a></p>"
						+ "\n" + "  </body>\n" + "</html>");
		htmlPane.setEditable(false);
		htmlPane.setOpaque(false);
		htmlPane.addHyperlinkListener(new HyperlinkListener() {
			public void hyperlinkUpdate(HyperlinkEvent e) {
				if (HyperlinkEvent.EventType.ACTIVATED.equals(e.getEventType()))
					try {
						BrowserLauncher.openURL(e.getURL().toString());
					} catch (IOException exception) {
						// ignore
					}
			}
		});

		JPanel panel = new JPanel(new BorderLayout());
		panel.add(htmlPane, BorderLayout.CENTER);

		final JFrame frame = new JFrame("About BoneJ...");
		frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		frame.getContentPane().add(panel, BorderLayout.CENTER);
		frame.pack();
		frame.setVisible(true);
	}

}
