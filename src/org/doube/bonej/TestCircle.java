package org.doube.bonej;

import ij.measure.ResultsTable;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
import org.doube.bonej.FitCircle;

/**
 * Check out how different circle fitting algorithms handle data.
 * 
 * @author Michael Doube
 *
 */
public class TestCircle implements PlugIn {
    double x;

    double y;

    double r;

    double sA;

    double eA;

    int n;

    double noise;

    boolean doKasa;

    boolean doPrattNTN;

    boolean doPrattSVD;

    boolean doHyperSimple;

    boolean doHyperStable;

    boolean doTaubinNTN;

    boolean doTaubinSVD;

    boolean doLevenMarqFull;

    boolean doLevenMarqRed;

    ResultsTable rt;

    public void run(String arg) {
	if (!showDialog()) {
	    return;
	}
	FitCircle fc = new FitCircle();
	double[][] testCircle = fc.getTestCircle(x, y, r, n, sA, eA, noise);

	rt = ResultsTable.getResultsTable();

	if (doKasa) {
	    double[] circle = fc.kasaFit(testCircle);
	    addResult(circle, "Kåsa");
	}
	if (doHyperStable) {
	    double[] circle = fc.hyperStable(testCircle);
	    addResult(circle, "Hyper Stable");
	}

	if (doHyperSimple) {
	    double[] circle = fc.hyperSimple(testCircle);
	    addResult(circle, "Hyper Simple");
	}

	if (doPrattNTN) {
	    double[] circle = fc.prattNewton(testCircle);
	    addResult(circle, "Pratt-Newton");
	}

	if (doPrattSVD) {
	    double[] circle = fc.prattSVD(testCircle);
	    addResult(circle, "Pratt-SVD");
	}
	if (doTaubinNTN) {
	    double[] circle = fc.taubinNewton(testCircle);
	    addResult(circle, "Taubin-Newton");
	}
	if (doTaubinSVD) {
	    double[] circle = fc.taubinSVD(testCircle);
	    addResult(circle, "Taubin-SVD");
	}
	if (doLevenMarqFull) {
	    double[] circle = fc.levenMarqFull(testCircle);
	    addResult(circle, "Levenburg-Marquardt (full)");
	}
	if (doLevenMarqRed) {
	    double[] circle = fc.levenMarqRed(testCircle);
	    addResult(circle, "Levenburg-Marquardt (reduced)");
	}
	return;
    }

    private void addResult(double[] circle, String string) {
	rt.incrementCounter();
	rt.addLabel(string);
	rt.addValue("X", circle[0]);
	rt.addValue("Y", circle[1]);
	rt.addValue("R", circle[2]);
	rt.show("Results");
	return;
    }

    private boolean showDialog() {
	GenericDialog gd = new GenericDialog("Options");
	gd.addMessage("Circle parameters");
	gd.addNumericField("Centre X", 0, 2);
	gd.addNumericField("Centre Y", 0, 2);
	gd.addNumericField("Radius", 10, 2);
	gd.addNumericField("Start angle", 0, 4);
	gd.addNumericField("End angle", 2 * Math.PI, 4);
	gd.addNumericField("N points", 50, 0);
	gd.addNumericField("Noise", 0.01, 2);

	gd.addMessage("Algorithms");
	gd.addCheckbox("Kåsa", true);
	gd.addCheckbox("Pratt-Newton", true);
	gd.addCheckbox("Pratt-SVD", true);
	gd.addCheckbox("Hyper (Simple)", true);
	gd.addCheckbox("Hyper (Stable)", true);
	gd.addCheckbox("Taubin-Newton", true);
	gd.addCheckbox("Taubin-SVD", true);
	gd.addCheckbox("Levenburg-Marquardt (Full)", true);
	gd.addCheckbox("Levenburg-Marquardt (Reduced)", true);

	gd.showDialog();
	if (gd.wasCanceled()) {
	    return false;
	} else {
	    x = gd.getNextNumber();
	    y = gd.getNextNumber();
	    r = gd.getNextNumber();
	    sA = gd.getNextNumber();
	    eA = gd.getNextNumber();
	    n = (int) gd.getNextNumber();
	    noise = gd.getNextNumber();

	    doKasa = gd.getNextBoolean();
	    doPrattNTN = gd.getNextBoolean();
	    doPrattSVD = gd.getNextBoolean();
	    doHyperSimple = gd.getNextBoolean();
	    doHyperStable = gd.getNextBoolean();
	    doTaubinNTN = gd.getNextBoolean();
	    doTaubinSVD = gd.getNextBoolean();
	    doLevenMarqFull = gd.getNextBoolean();
	    doLevenMarqRed = gd.getNextBoolean();
	    return true;
	}
    }
}
