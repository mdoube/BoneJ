/**
 *Neck Shaft Angle ImageJ plugin
 *Copyright 2008 2009 Michael Doube 
 *
 *This program is free software: you can redistribute it and/or modify
 *it under the terms of the GNU General Public License as published by
 *the Free Software Foundation, either version 3 of the License, or
 *(at your option) any later version.
 *
 *This program is distributed in the hope that it will be useful,
 *but WITHOUT ANY WARRANTY; without even the implied warranty of
 *MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *GNU General Public License for more details.
 *
 *You should have received a copy of the GNU General Public License
 *along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ImageProcessor;
import ij.plugin.filter.PlugInFilter;
import ij.plugin.frame.RoiManager;
import ij.measure.Calibration;
import ij.gui.*;

import java.awt.Rectangle;
import java.awt.event.*;

import org.doube.bonej.FitSphere;
import org.doube.bonej.ResultInserter;
import org.doube.bonej.FitCircle;
import org.doube.jama.*;

/*
 * TODO incorporate curvature
 * curvature = radius of circle fit to diaphyseal slice centroids / bone length
 */
/**
 *<p>
 * Neck Shaft Angle<br />
 * 
 * Tool to calculate the neck shaft angle of 3D images of femora.
 * </p>
 * <p>
 * Neck shaft angle is the angle formed at the intersection between the coplanar
 * lines <b>N</b> and <b>S</b>. <br>
 * <b>S</b> is the singular value decomposition (orthogonal distance regression)
 * vector that passes through the centroid (<b>B</b>) of the bone and describes
 * the long axis of the bone.<br/>
 * <b>C</b> is the centre of a sphere fit to the femoral head.<br/>
 * 
 * <b>P</b> is the plane that contains <b>S</b> and <b>C</b>. <br />
 * <b>N</b> is the projection onto <b>P</b> of a vector originating at <b>C</b>
 * and passing through the 'middle' of the femoral neck.<br />
 * 
 * Singular value decomposition performed with the <a
 * href="http://math.nist.gov/javanumerics/jama/">Jama</a> package
 * </p>
 * 
 *@author Michael Doube
 *@version 0.1
 */
public class Neck_Shaft_Angle implements PlugInFilter, MouseListener {
    ImagePlus imp;

    ImageCanvas canvas;

    protected ImageStack stack;

    public float[] CTable;

    public double[] coeff = { 0, 1 }, neckPoint = { 78.75, 100.55, 80 },
	    headCentre, centroid;

    public double[][] shaftVector;

    public int minT = 0, maxT = 4000; // min and maximum bone value in HU

    public int startSlice = 1, endSlice;

    public boolean doCurvature;

    public String title, units, valueUnit;

    public Calibration cal;

    public int setup(String arg, ImagePlus imp) {
	if (imp == null || imp.getNSlices() < 2) {
	    IJ.showMessage("A stack must be open");
	    return DONE;
	}
	this.imp = imp;
	this.stack = imp.getStack();
	this.endSlice = stack.getSize();
	this.cal = imp.getCalibration();
	this.coeff = cal.getCoefficients();
	this.title = imp.getShortTitle();
	return DOES_16 + STACK_REQUIRED;
    }

    public void run(ImageProcessor ip) {
	// set up pixel calibration
	if (!this.cal.isSigned16Bit() && !this.cal.calibrated()) {
	    IJ.run("Threshold...");
	    new WaitForUserDialog(
		    "This image is not density calibrated.\nSet the threshold, then click OK.")
		    .show();
	    this.minT = (int) ip.getMinThreshold();
	    this.maxT = (int) ip.getMaxThreshold();
	    showDialog(valueUnit);
	    IJ.log("Image is uncalibrated: using user-determined threshold "
		    + minT + " to " + maxT);
	} else if (this.coeff[0] == -1000 && this.coeff[1] == 1.0) {
	    // looks like an HU calibrated image
	    // convert HU limits to pixel values
	    showDialog("Hounsfield units");
	    this.minT = (int) Math.round(cal.getRawValue(minT));
	    this.maxT = (int) Math.round(cal.getRawValue(maxT));
	    IJ.log("Image looks like it is HU calibrated. Using " + minT
		    + " and " + maxT + " " + valueUnit + " as bone cutoffs");
	} else if (cal.isSigned16Bit() && !cal.calibrated()) {
	    new WaitForUserDialog(
		    "This image is not density calibrated.\nSet the threshold, then click OK.")
		    .show();
	    this.minT = (int) ip.getMinThreshold();
	    this.maxT = (int) ip.getMaxThreshold();
	    showDialog(valueUnit);
	    IJ.log("Image is uncalibrated: using user-determined threshold "
		    + minT + " to " + maxT);
	} else {
	    IJ.error("Unrecognised file type");
	    return;
	}
	// get coordinates from the ROI manager and fit a sphere
	RoiManager roiMan = RoiManager.getInstance();
	if (roiMan == null && imp != null) {
	    IJ.run("ROI Manager...");
	    IJ.error("Please populate ROI Manager with point ROIs\n"
		    + "placed on the boundary of the femoral head");
	    return;
	} else {
	    FitSphere fs = new FitSphere();
	    double[][] points = fs.getRoiManPoints(imp, roiMan);
	    this.headCentre = fs.fitSphere(points);
	}
	ImageWindow win = this.imp.getWindow();
	this.canvas = win.getCanvas();

	// work out the centroid and regression vector of the bone
	this.centroid = findCentroid3D(this.stack, startSlice, endSlice);
	if (this.centroid[0] < 0) {
	    IJ
		    .error(
			    "Empty Stack",
			    "No voxels are a"
				    + "vailable for calculation.\nCheck your ROI and threshold.");
	    return;
	}

	this.shaftVector = regression3D(this.stack, this.centroid);

	if (doCurvature)
	    calculateCurvature(this.stack, this.shaftVector, this.headCentre,
		    this.centroid);

	// remove stale MouseListeners
	MouseListener[] l = this.canvas.getMouseListeners();
	for (int n = 0; n < l.length; n++) {
	    this.canvas.removeMouseListener(l[n]);
	}
	// add a new MouseListener
	this.canvas.addMouseListener(this);

	new WaitForUserDialog("Click on the middle of the femoral neck.\n"
		+ "Neck-shaft angle and out-of-plane skew\n"
		+ "will be recorded until you hit \'OK\'").show();
	this.canvas.removeMouseListener(this);
	return;
    }

    /**
     * Find the centroid of the voxels of interest
     * 
     * @param stack
     * @param startSlice
     * @param endSlice
     * @param rt
     * @return centroid double[3]
     */
    public double[] findCentroid3D(ImageStack stack, int startSlice,
	    int endSlice) {
	Rectangle r = stack.getRoi();
	int offset, i;
	int w = stack.getWidth();
	int h = stack.getHeight();
	int sliceSize = w * h;
	double sumx = 0;
	double sumy = 0;
	double sumz = 0;
	double count = 0;
	for (int z = startSlice; z <= endSlice; z++) {
	    IJ.showStatus("Calculating centroid...");
	    IJ.showProgress(z, endSlice);
	    short[] slicePixels = new short[sliceSize];
	    slicePixels = (short[]) stack.getPixels(z);
	    for (int y = r.y; y < (r.y + r.height); y++) {
		offset = y * w;
		for (int x = r.x; x < (r.x + r.width); x++) {
		    i = offset + x;
		    int testPixel = slicePixels[i] & 0xffff;
		    if (testPixel >= minT && testPixel <= maxT) {
			sumx += (double) x;
			sumy += (double) y;
			sumz += (double) z;
			count++;
		    }
		}
	    }
	}
	// centroid in pixels
	double centX = sumx / count;
	double centY = sumy / count;
	double centZ = sumz / count;
	if (sumx == 0) {
	    centX = -1;
	    centY = -1;
	    centZ = -1;
	}
	// centroid in real units
	double[] centroid = { centX * cal.pixelWidth, centY * cal.pixelHeight,
		centZ * cal.pixelDepth };
	return centroid;
    }/* end findCentroid3D */

    /**
     * Calculate the vector associated with the projection plane from the
     * regression vector and the vector connecting the centroid and the femoral
     * head centre
     * 
     * @param shaftVector
     *            double[][]
     * @param headCentre
     *            double[][]
     * @param centroid
     *            double[][]
     * @return double[][] projectionPlane
     */
    public double[][] projectionPlane(double[][] shaftVector,
	    double[] headCentre, double[] centroid) {
	// have to calculate distance between points so that we find a unit
	// vector
	double d = Math.sqrt((headCentre[0] - centroid[0])
		* (headCentre[0] - centroid[0]) + (headCentre[1] - centroid[1])
		* (headCentre[1] - centroid[1]) + (headCentre[2] - centroid[2])
		* (headCentre[2] - centroid[2]));
	double[][] cHVec = new double[3][1];
	cHVec[0][0] = (headCentre[0] - centroid[0]) / d;
	cHVec[1][0] = (headCentre[1] - centroid[1]) / d;
	cHVec[2][0] = (headCentre[2] - centroid[2]) / d;

	Matrix cH = new Matrix(cHVec);
	printMatrix(cH, "cHVec");

	// projectionPlane is the cross product of cHVec and shaftVector
	double[][] projectionPlane = new double[3][1];
	projectionPlane[0][0] = cHVec[1][0] * shaftVector[2][0] - cHVec[2][0]
		* shaftVector[1][0];
	projectionPlane[1][0] = cHVec[2][0] * shaftVector[0][0] - cHVec[0][0]
		* shaftVector[2][0];
	projectionPlane[2][0] = cHVec[0][0] * shaftVector[1][0] - cHVec[1][0]
		* shaftVector[0][0];

	d = Math.sqrt(projectionPlane[0][0] * projectionPlane[0][0]
		+ projectionPlane[1][0] * projectionPlane[1][0]
		+ projectionPlane[2][0] * projectionPlane[2][0]);
	projectionPlane[0][0] /= d;
	projectionPlane[1][0] /= d;
	projectionPlane[2][0] /= d;

	return projectionPlane;
    }

    /**
     * Calculate the vector associated with the plane formed between neckVector
     * and the normal to projectionPlane from the regression vector and the
     * vector connecting the centroid and the femoral head centre
     * 
     * @param projectionPlane
     *            double[][]
     * @param neckVector
     *            double[][]
     * @return double[][] neckPlane
     */
    public double[][] neckPlane(double[][] projectionPlane,
	    double[][] neckVector) {
	// neckPlane is the cross product of neckVector and projectionPlane
	double[][] neckPlane = new double[3][1];
	neckPlane[0][0] = projectionPlane[1][0] * neckVector[2][0]
		- projectionPlane[2][0] * neckVector[1][0];
	neckPlane[1][0] = projectionPlane[2][0] * neckVector[0][0]
		- projectionPlane[0][0] * neckVector[2][0];
	neckPlane[2][0] = projectionPlane[0][0] * neckVector[1][0]
		- projectionPlane[1][0] * neckVector[0][0];
	return neckPlane;
    }

    /**
     * Find the intersection between neckPlane and projectionPlane
     * 
     * @param projectionPlane
     * @param neckPlane
     * @return double[][] testVector.
     */
    public double[][] testVector(double[][] projectionPlane,
	    double[][] neckPlane) {
	// testVector is the cross product of neckPlane and projectionPlane
	double[][] testVector = new double[3][1];
	testVector[0][0] = projectionPlane[1][0] * neckPlane[2][0]
		- projectionPlane[2][0] * neckPlane[1][0];
	testVector[1][0] = projectionPlane[2][0] * neckPlane[0][0]
		- projectionPlane[0][0] * neckPlane[2][0];
	testVector[2][0] = projectionPlane[0][0] * neckPlane[1][0]
		- projectionPlane[1][0] * neckPlane[0][0];
	return testVector;
    }/* end testVector */

    /**
     * Calculate the orthogonal distance regression plane of a set of points by
     * the covariance method and Singular Value Decomposition
     * 
     * @param stack
     * @param centroid
     * @return SingularValueDecomposition containing eigenvector and eigenvalue
     * 
     * @see <a
     *      href="http://mathforum.org/library/drmath/view/63765.html">Description
     *      on Ask Dr Math</a>
     * 
     */
    public double[][] regression3D(ImageStack stack, double[] centroid) {
	IJ.showStatus("Calculating SVD");
	Rectangle r = stack.getRoi();
	IJ.log("Rectangle r has top left coordinates of (" + r.x + ", " + r.y
		+ ") and size of " + r.width + " x " + r.height);
	int offset, i;
	int w = stack.getWidth();
	int h = stack.getHeight();
	int sliceSize = w * h;
	double vW = cal.pixelWidth;
	double vH = cal.pixelHeight;
	double vD = cal.pixelDepth;
	double[][] C = new double[3][3];
	double count = 0;
	for (int z = startSlice; z <= endSlice; z++) {
	    IJ.showStatus("Calculating centroid...");
	    IJ.showProgress(z, endSlice);
	    short[] slicePixels = new short[sliceSize];
	    slicePixels = (short[]) stack.getPixels(z);
	    for (int y = r.y; y < (r.y + r.height); y++) {
		offset = y * w;
		for (int x = r.x; x < (r.x + r.width); x++) {
		    i = offset + x;
		    int testPixel = slicePixels[i] & 0xffff;
		    if (testPixel >= minT && testPixel <= maxT) {
			double dx = x * vW - centroid[0];
			double dy = y * vH - centroid[1];
			double dz = z * vD - centroid[2];
			C[0][0] += dx * dx;
			C[1][1] += dy * dy;
			C[2][2] += dz * dz;
			C[0][1] += dx * dy;
			C[0][2] += dx * dz;
			C[1][0] += dx * dy;
			C[1][2] += dy * dz;
			C[2][0] += dx * dz;
			C[2][1] += dy * dz;
			count += 1;
		    }
		}
	    }
	}
	double invCount = 1 / count;
	Matrix covarianceMatrix = new Matrix(C).times(invCount);
	printMatrix(covarianceMatrix, "Covariance matrix");
	SingularValueDecomposition S = new SingularValueDecomposition(
		covarianceMatrix);
	Matrix leftVectors = S.getU();
	printMatrix(leftVectors, "Left vectors");
	double[][] orthogonalDistanceRegression = new double[3][1];
	orthogonalDistanceRegression[0][0] = leftVectors.get(0, 0);
	orthogonalDistanceRegression[1][0] = leftVectors.get(1, 0);
	orthogonalDistanceRegression[2][0] = leftVectors.get(2, 0);
	return orthogonalDistanceRegression;
    }/* end Regression3D */

    public double[][] neckVector(double[] headCentre, double[] neckPoint) {
	// have to calculate d to make sure that neckVector is a unit vector
	double d = Math.sqrt((headCentre[0] - neckPoint[0])
		* (headCentre[0] - neckPoint[0])
		+ (headCentre[1] - neckPoint[1])
		* (headCentre[1] - neckPoint[1])
		+ (headCentre[2] - neckPoint[2])
		* (headCentre[2] - neckPoint[2]));

	double[][] neckVector = new double[3][1];
	neckVector[0][0] = (headCentre[0] - neckPoint[0]) / d;
	neckVector[1][0] = (headCentre[1] - neckPoint[1]) / d;
	neckVector[2][0] = (headCentre[2] - neckPoint[2]) / d;
	return neckVector;
    }

    public void calculateAngles() {
	double[][] neckVector = neckVector(headCentre, neckPoint);
	double[][] projectionPlane = projectionPlane(shaftVector, headCentre,
		centroid);
	double[][] neckPlane = neckPlane(neckVector, projectionPlane);
	double[][] testVector = testVector(projectionPlane, neckPlane);
	// P . Q = ||P|| ||Q|| cos(a) so if P and Q are unit vectors, then P.Q =
	// cos(a)
	Matrix pP = new Matrix(projectionPlane);
	printMatrix(pP, "projectionPlane");

	Matrix tV = new Matrix(testVector);
	printMatrix(tV, "testVector");

	Matrix sV = new Matrix(shaftVector);
	printMatrix(sV, "shaftVector");

	Matrix nV = new Matrix(neckVector);
	printMatrix(nV, "neckVector");

	double cosA1 = sV.get(0, 0) * tV.get(0, 0) + sV.get(1, 0)
		* tV.get(1, 0) + sV.get(2, 0) * tV.get(2, 0);
	// printMatrix(cosA1, "cosA1");
	IJ.log("cosA1: " + cosA1);

	double cosA2 = nV.get(0, 0) * tV.get(0, 0) + nV.get(1, 0)
		* tV.get(1, 0) + nV.get(2, 0) * tV.get(2, 0);
	// printMatrix(cosA2, "cosA2");
	IJ.log("cosA2: " + cosA2);

	double neckShaftAngle = Math.acos(cosA1);
	double neckShaftSkew = Math.acos(cosA2);
	ResultInserter ri = ResultInserter.getInstance();
	ri.setResultInRow(this.imp, "Angle (rad)", neckShaftAngle);
	ri.setResultInRow(this.imp, "Skew (rad)", neckShaftSkew);
	ri.updateTable();
    }

    /**
     * <p>
     * Calculate curvature of bone using shaft vector as a reference axis and
     * centre of femoral head to define reference plane
     * </p>
     * 
     * @param stack
     * @param shaftVector
     * @param headCentre
     */
    private void calculateCurvature(ImageStack stack, double[][] shaftVector,
	    double[] headCentre, double[] centroid) {
	// calculate the eigenvector of the reference plane containing
	// the shaftVector and the headCentre

	// get the 2D centroids
	double vW = this.cal.pixelWidth;
	double vH = this.cal.pixelHeight;
	Rectangle r = stack.getRoi();
	// pixel counters
	double cstack = 0;
	int w = stack.getWidth();

	boolean[] emptySlices = new boolean[this.stack.getSize()];
	double[] cortArea = new double[this.stack.getSize()];
	double[][] sliceCentroids = new double[2][this.stack.getSize()];

	double pixelArea = this.cal.pixelWidth * this.cal.pixelHeight;
	for (int s = this.startSlice; s <= this.endSlice; s++) {
	    double sumX = 0;
	    double sumY = 0;
	    double cslice = 0;
	    short[] pixels = (short[]) this.stack.getPixels(s);
	    for (int y = r.y; y < (r.y + r.height); y++) {
		int offset = y * w;
		for (int x = r.x; x < (r.x + r.width); x++) {
		    int i = offset + x;
		    if (pixels[i] >= this.minT && pixels[i] <= this.maxT) {
			cslice++;
			cortArea[s] += pixelArea;
			sumX += x * vW;
			sumY += y * vH;
		    }
		}
	    }
	    if (cslice > 0) {
		sliceCentroids[0][s] = sumX / cslice;
		sliceCentroids[1][s] = sumY / cslice;
		cstack += cslice;
		emptySlices[s] = false;
	    } else {
		emptySlices[s] = true;
	    }
	}

	double[][] projPlane = projectionPlane(shaftVector, headCentre,
		centroid);
	double pPx = projPlane[0][0];
	double pPy = projPlane[1][0];
	double pPz = projPlane[2][0];

	double x1x = this.centroid[0];
	double x1y = this.centroid[1];
	double x1z = this.centroid[2];
	double x2x = x1x + shaftVector[0][0];
	double x2y = x1y + shaftVector[1][0];
	double x2z = x1z + shaftVector[2][0];

	// for each centroid, calculate the vector to the 3D regression line
	// using equation 10 from
	// http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
	// with denominator = 1 because we use a unit vector for |(x2 - x1)|

	double[][] mL = new double[this.endSlice - this.startSlice + 1][2];
	double[][] cC = new double[this.endSlice - this.startSlice + 1][2];
	int i = 0;
	for (int s = this.startSlice; s <= this.endSlice; s++) {
	    if (!emptySlices[s]) {
		double x0x = sliceCentroids[0][s];
		double x0y = sliceCentroids[1][s];
		double x0z = s * this.cal.pixelDepth;

		// distance is magnitude of cross product of (x0 - x1) and (x0 -
		// x2)
		double x0x1x = x0x - x1x;
		double x0x1y = x0y - x1y;
		double x0x1z = x0z - x1z;

		double x0x2x = x0x - x2x;
		double x0x2y = x0y - x2y;
		double x0x2z = x0z - x2z;

		double cpX = x0x1y * x0x2z - x0x1z * x0x2y;
		double cpY = x0x1z * x0x2x - x0x1x * x0x2z;
		double cpZ = x0x1x * x0x2y - x0x1y * x0x2x;

		double distance = Math.sqrt(cpX * cpX + cpY * cpY + cpZ * cpZ);
		// IJ.log("distance to regression line is "+ distance +
		// " "+this.units+" for slice "+s);

		// work out t (number of unit vectors from centroid along
		// regression)
		// as per equation 3

		double x1x0x = x1x - x0x;
		double x1x0y = x1y - x0y;
		double x1x0z = x1z - x0z;

		double x2x1x = x2x - x1x;
		double x2x1y = x2y - x1y;
		double x2x1z = x2z - x1z;

		double t = -1 * (x1x0x * x2x1x + x1x0y * x2x1y + x1x0z * x2x1z);

		// So now the intersection point x3 of the perpendicular is
		// known
		// as centroid + t * unitVector
		// and the vector to the deflection as (x0 - x3)

		double x3x = x1x + t * shaftVector[0][0];
		double x3y = x1y + t * shaftVector[1][0];
		double x3z = x1z + t * shaftVector[2][0];

		double defVectX = x0x - x3x;
		double defVectY = x0y - x3y;
		double defVectZ = x0z - x3z;

		// project deflection vector onto projection plane vector by
		// taking the dot product
		// this is the craniocaudal deflection

		double cranioCaudal = (defVectX * pPx + defVectY * pPy + defVectZ
			* pPz);
		cC[i][0] = cranioCaudal;
		cC[i][1] = t;
		// IJ.log("Craniocaudal deflection at slice "+s+" is "+cranioCaudal);

		// mediolateral deflection is distance in projectionPlane, i.e.
		// deflection projected onto projectionPlane // double cross
		// product
		// B x (A x B), provided that B is a unit vector

		double aBx = defVectY * pPz - defVectZ * pPy;
		double aBy = defVectZ * pPx - defVectX * pPz;
		double aBz = defVectX * pPy - defVectY * pPx;

		double mLx = pPy * aBz - pPz * aBy;
		double mLy = pPz * aBx - pPx * aBz;
		double mLz = pPx * aBy - pPy * aBx;

		double medioLateral = Math.sqrt(mLx * mLx + mLy * mLy + mLz
			* mLz);

		// give the scalar a direction
		double sign = (mLx * mLy * mLz) / Math.abs(mLx * mLy * mLz);

		medioLateral *= sign;

		// IJ.log("Mediolateral deflection at slice "+s+" is "+medioLateral);
		IJ.log(s + "," + t + ", " + distance + ", " + medioLateral
			+ ", " + cranioCaudal);
		mL[i][0] = medioLateral;
		mL[i][1] = t;
		i++;
	    } else {
		// IJ.log("No pixels to calculate centroid in slice "+s);
	    }
	}
	// Calculate circle fitting for mL and cC deflections
	ResultInserter ri = ResultInserter.getInstance();
	String units = this.cal.getUnits();
	
	FitCircle fc = new FitCircle();
	double[] mLabR = fc.hyperStable(mL);
	double[] cCabR = fc.hyperStable(cC);
	ri.setResultInRow(imp, "M-L radius ("+units+")", mLabR[2]);
	ri.setResultInRow(imp, "M-L centre X ("+units+")", mLabR[0]);
	ri.setResultInRow(imp, "M-L centre Y ("+units+")", mLabR[1]);
	ri.setResultInRow(imp, "Cr-Ca radius ("+units+")", cCabR[2]);
	ri.setResultInRow(imp, "Cr-Ca centre X ("+units+")", cCabR[0]);
	ri.setResultInRow(imp, "Cr-Ca centre Y ("+units+")", cCabR[1]);
	ri.updateTable();
	return;
    }

    public void mousePressed(MouseEvent e) {
	int x = canvas.offScreenX(e.getX());
	int y = canvas.offScreenY(e.getY());
	int z = imp.getCurrentSlice();
	neckPoint[0] = x * cal.pixelWidth;
	neckPoint[1] = y * cal.pixelHeight;
	neckPoint[2] = z * cal.pixelDepth;
	IJ.log("neckPoint: (" + neckPoint[0] + "," + neckPoint[1] + ", "
		+ neckPoint[2] + ")");
	calculateAngles();
    }

    public void mouseReleased(MouseEvent e) {
    }

    public void mouseExited(MouseEvent e) {
    }

    public void mouseClicked(MouseEvent e) {
    }

    public void mouseEntered(MouseEvent e) {
    }

    public void mouseMoved(MouseEvent e) {
    }

    public void showDialog(String valueUnit) {
	GenericDialog gd = new GenericDialog("Setup");
	gd.addNumericField("Shaft Start Slice:", startSlice, 0);
	gd.addNumericField("Shaft End Slice:", endSlice, 0);
	gd.addMessage("Only use pixels between clip values:");
	gd.addNumericField("Clip min.", minT, 0, 6, valueUnit);
	gd.addNumericField("Clip max.", maxT, 0, 6, valueUnit);
	gd.addCheckbox("Calculate curvature", true);
	gd.showDialog();
	if (gd.wasCanceled()) {
	    return;
	}
	startSlice = (int) gd.getNextNumber();
	endSlice = (int) gd.getNextNumber();
	minT = (int) gd.getNextNumber();
	maxT = (int) gd.getNextNumber();
	doCurvature = gd.getNextBoolean();
    }

    public void printMatrix(Matrix matrix, String title) {
	IJ.log(title);
	int nCols = matrix.getColumnDimension();
	int nRows = matrix.getRowDimension();
	double[][] eVal = matrix.getArrayCopy();
	for (int r = 0; r < nRows; r++) {
	    String row = "||";
	    for (int c = 0; c < nCols; c++) {
		row = row + eVal[r][c] + "|";
	    }
	    row = row + "|";
	    IJ.log(row);
	}
	return;
    }
}
