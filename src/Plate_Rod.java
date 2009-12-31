
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;
import ij.plugin.filter.PlugInFilter;
import ij.gui.GenericDialog;
import ij.macro.Interpreter;
import ij.measure.Calibration;

import org.doube.jama.*;
import org.doube.util.ImageCheck;
import org.doube.util.ResultInserter;

/**
 * <p>
 * <b>Plate_Rod</b>
 * </p>
 * <p>
 * ImageJ plugin to describe the local geometry of a binary image in an
 * oblate/prolate spheroid space. Uses Skeletonize3D to generate a 3D skeleton,
 * the points of which are used as centres for star volumes. Local geometry is
 * determined by the ratio between the first and second eigenvalues and first
 * and third eigenvalues of each star volume.
 * </p>
 * 
 * @author Michael Doube
 * 
 */
/*
 * TODO Summarise plateness / rodness, e.g. ratio of sums of middle and biggest
 * Eigenvalues plot ev2/ev1 on 3d skeleton
 */

public class Plate_Rod implements PlugInFilter {
    private ImagePlus imp;

    protected int nVectors = 1000;

    protected ImageStack stack;

    protected double vW, vH, vD, samplingIncrement;

    protected String units;

    public int setup(String arg, ImagePlus imp) {
	if (imp == null || imp.getNSlices() < 2) {
	    IJ.showMessage("A stack must be open");
	    return DONE;
	}
	stack = imp.getStack();
	this.imp = imp;
	Calibration cal = imp.getCalibration();
	vW = cal.pixelWidth;
	vH = cal.pixelHeight;
	vD = cal.pixelDepth;
	units = cal.getUnits();
	if (this.imp != null
		&& (this.imp.getType() == ImagePlus.GRAY8 || this.imp.getType() == ImagePlus.COLOR_256)) {
	    ImageStatistics stats = this.imp.getStatistics();
	    if (stats.histogram[0] + stats.histogram[255] != stats.pixelCount) {
		IJ.error("8-bit binary (black and white only) image required.");
		return DONE;
	    }
	}
	return DOES_8G + STACK_REQUIRED;
    }

    public void run(ImageProcessor ip) {
    	if (!ImageCheck.checkIJVersion())
			return;
	samplingIncrement = Math.max(vH, Math.max(vW, vD)); // need to get
	// this from a
	// dialog
	showDialog();
	double[][] randomVectors = randomVectors(nVectors);
	double[][] skeletonPoints = skeletonPoints(this.imp);
	double[][] localEigenValues = localEigenValues(stack, randomVectors,
		skeletonPoints, samplingIncrement);

	//	ResultsTable rt = ResultsTable.getResultsTable();
	double sumEv1 = 0, sumEv2 = 0, sumEv3 = 0;
	int NaNs = 0;
	for (int l = 0; l < localEigenValues.length; l++) {
	    //	  i.e. so long as it's not NaN
	    if (localEigenValues[l][0] >= 0 || localEigenValues[l][0] < 0) { 
		sumEv1 += localEigenValues[l][0];
		sumEv2 += localEigenValues[l][1];
		sumEv3 += localEigenValues[l][2];
	    } else {
		NaNs++;
	    }
	}
	IJ.log(NaNs + " tests hit the sides");

	ResultInserter ri = ResultInserter.getInstance();
	ri.setResultInRow(this.imp, "ΣeV1", sumEv1);
	ri.setResultInRow(this.imp, "ΣeV2", sumEv2);
	ri.setResultInRow(this.imp, "ΣeV3", sumEv3);
	ri.setResultInRow(this.imp, "eV2/eV1", sumEv2 / sumEv1);
	ri.setResultInRow(this.imp, "eV3/eV1", sumEv3 / sumEv1);
	ri.updateTable();
    }

    /* ----------------------------------------------------------------------- */
    /**
     * Generate an array of randomly-oriented 3D unit vectors
     * 
     * @param nVectors
     *                number of vectors to generate
     * @return 2D array (nVectors x 3) containing unit vectors
     */
    public double[][] randomVectors(int nVectors) {
	double[][] randomVectors = new double[nVectors][3];
	for (int n = 0; n < nVectors; n++) {
	    randomVectors[n][2] = 2 * Math.random() - 1;
	    double rho = Math.sqrt(1 - randomVectors[n][2]
	                                                * randomVectors[n][2]);
	    double phi = Math.PI * (2 * Math.random() - 1);
	    randomVectors[n][0] = rho * Math.cos(phi);
	    randomVectors[n][1] = rho * Math.sin(phi);
	}
	return randomVectors;
    } /* end randomVectors */

    public double[][] skeletonPoints(ImagePlus imp) {
	ImageStack skeletonStack = new ImageStack(this.imp.getWidth(), this.imp
		.getHeight(), this.imp.getNSlices());
	for (int s = 1; s <= stack.getSize(); s++) {
	    byte[] pixels = (byte[]) stack.getPixels(s);
	    byte[] pixelsCopy = pixels.clone();
	    skeletonStack.setPixels(pixelsCopy, s);
	}

	ImagePlus skeletonImp = new ImagePlus("Skeleton", skeletonStack);
	skeletonImp.setCalibration(this.imp.getCalibration());
	skeletonImp.show();
	IJ.run("Invert LUT");
	IJ.run("Skeletonise 3D");
	skeletonStack = skeletonImp.getStack();
	int d = skeletonImp.getStackSize();
	int h = skeletonImp.getHeight();
	int w = skeletonImp.getWidth();
	double vW = this.imp.getCalibration().pixelWidth;
	double vH = this.imp.getCalibration().pixelHeight;
	double vD = this.imp.getCalibration().pixelDepth;
	int count = 0;
	for (int z = 1; z <= d; z++) {
	    byte[] slicePixels = (byte[]) skeletonStack.getPixels(z);
	    for (int y = 0; y < h; y++) {
		int offset = y * w;
		for (int x = 0; x < w; x++) {
		    if (slicePixels[offset + x] < 0) {
			count++;
		    }
		}
	    }
	}
	IJ.log("Counted " + count + " skeleton points");
	double[][] skeletonPoints = new double[count][3];
	int p = 0;
	for (int z = 0; z < d; z++) {
	    byte[] slicePixels = (byte[]) skeletonStack.getPixels(z + 1);
	    for (int y = 0; y < h; y++) {
		int offset = y * w;
		for (int x = 0; x < w; x++) {
		    if (slicePixels[offset + x] < 0) {
			skeletonPoints[p][0] = (double) x * vW;
			skeletonPoints[p][1] = (double) y * vH;
			skeletonPoints[p][2] = (double) z * vD;
			p++;
		    }
		}
	    }
	}
	return skeletonPoints;
    }

    public double[][] localEigenValues(ImageStack stack,
	    double[][] randomVectors, double[][] skeletonPoints,
	    double samplingIncrement) {
	double[][] localEigenValues = new double[skeletonPoints.length][3];
	int w = stack.getWidth();
	int h = stack.getHeight();
	int d = stack.getSize();
	int nP = skeletonPoints.length;
	int nV = randomVectors.length;
	// make a work array containing the stack
	byte[] workArray = new byte[w * h * d];
	int pixPerSlice = w * h;
	for (int s = 0; s < d; s++) {
	    byte[] slicePixels = (byte[]) stack.getPixels(s + 1); // slice
	    // number
	    // starts
	    // at 1,
	    // array
	    // starts
	    // at 0
	    System.arraycopy(slicePixels, 0, workArray, s * pixPerSlice,
		    pixPerSlice);
	}

	for (int p = 0; p < nP; p++) {
	    IJ.showStatus("Calculating local eigenvalues");
	    IJ.showProgress(p, nP);
	    boolean hitSide = false;
	    double[][] localStar = new double[nV][3];
	    for (int v = 0; v < nV; v++) {
		if (hitSide)
		    break;
		byte pixelValue = -1;
		double vectorLength = 0;
		while (pixelValue < 0 && !hitSide) {
		    int testPixelX = (int) Math
		    .floor((skeletonPoints[p][0] + randomVectors[v][0]
		                                                    * vectorLength)
		                                                    / vW);
		    int testPixelY = (int) Math
		    .floor((skeletonPoints[p][1] + randomVectors[v][1]
		                                                    * vectorLength)
		                                                    / vH);
		    int testPixelZ = (int) Math
		    .floor((skeletonPoints[p][2] + randomVectors[v][2]
		                                                    * vectorLength)
		                                                    / vD);
		    if (testPixelX < 0 || testPixelX >= w || testPixelY < 0
			    || testPixelY >= h || testPixelZ < 0
			    || testPixelZ >= d) {
			hitSide = true;
			break;
		    } else {
			pixelValue = workArray[testPixelZ * pixPerSlice
			                       + testPixelY * w + testPixelX];
			vectorLength += samplingIncrement;
		    }
		}
		localStar[v][0] = vectorLength * randomVectors[v][0];
		localStar[v][1] = vectorLength * randomVectors[v][1];
		localStar[v][2] = vectorLength * randomVectors[v][2];
	    }
	    if (!hitSide) {
		EigenvalueDecomposition E = PrincipalComponents(localStar);
		localEigenValues[p][0] = E.getD().get(2, 2);
		localEigenValues[p][1] = E.getD().get(1, 1);
		localEigenValues[p][2] = E.getD().get(0, 0);
	    } else {
		localEigenValues[p][0] = 0.0 / 0.0;
		localEigenValues[p][1] = 0.0 / 0.0;
		localEigenValues[p][2] = 0.0 / 0.0;
	    }
	}
	return localEigenValues;
    }

    /*---------------------------------------------------------*/
    /**
     * Calculate the eigenvectors and eigenvalues of a set of points by the
     * covariance method and eigendecomposition.
     * 
     * @param coOrdinates
     *                n x 3 array centred on (0,0,0)
     * @return EigenvalueDecomposition containing eigenvectors and
     *         eigenvalues
     * 
     */
    public EigenvalueDecomposition PrincipalComponents(double[][] coOrdinates) {
	double sumX = 0, sumY = 0, sumZ = 0;
	for (int n = 0; n < coOrdinates.length; n++) {
	    sumX += coOrdinates[n][0];
	    sumY += coOrdinates[n][1];
	    sumZ += coOrdinates[n][2];
	}
	double centX = sumX / coOrdinates.length;
	double centY = sumY / coOrdinates.length;
	double centZ = sumZ / coOrdinates.length;

	double[][] C = new double[3][3];
	double count = 0;
	for (int n = 0; n < coOrdinates.length; n++) {
	    double x = coOrdinates[n][0] - centX;
	    double y = coOrdinates[n][1] - centY;
	    double z = coOrdinates[n][2] - centZ;
	    C[0][0] += x * x;
	    C[1][1] += y * y;
	    C[2][2] += z * z;
	    C[0][1] += x * y;
	    C[0][2] += x * z;
	    C[1][0] += x * y;
	    C[1][2] += y * z;
	    C[2][0] += x * z;
	    C[2][1] += y * z;
	    count += 1;
	}
	double invCount = 1 / count;
	Matrix covarianceMatrix = new Matrix(C).times(invCount);
	EigenvalueDecomposition E = new EigenvalueDecomposition(
		covarianceMatrix);
	return E;
    }/* end PrincipalComponents */

    public boolean showDialog() {
	GenericDialog gd = new GenericDialog("Setup");
	gd
	.addNumericField("Sampling increment", samplingIncrement, 3, 8,
		units);
	gd.addNumericField("Vectors", nVectors, 0, 8, "");
	gd.showDialog();
	if (gd.wasCanceled()) {
	    return false;
	} else {
	    if (!Interpreter.isBatchMode()){
		samplingIncrement = gd.getNextNumber();
		nVectors = (int) Math.round(gd.getNextNumber());
		return true;
	    } else {
		return true;
	    }
	}
    }
}
