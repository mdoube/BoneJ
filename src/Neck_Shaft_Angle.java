

/** Neck Shaft Angle
 *Tool to calculate the neck shaft angle
 *
 * Neck shaft angle is the angle formed at the intersection
 * between the coplanar lines N and S
 * 
 * S is the singular value decomposition (orthogonal distance regression) vector 
 * that passes through the centroid (B) of the bone and describes 
 * the long axis of the bone  
 * 
 * C is the centre of a sphere fit to the femoral head.
 * 
 * P is the plane that contains S and C
 * 
 * N is the projection onto P of a vector originating at C
 * and passing through the 'middle' of the femoral neck.
 *
 *Singular value decomposition performed with the Jama package
 *http://math.nist.gov/javanumerics/jama/
 *Unfortunately JAMA does not have dot or cross product methods!
 *
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
 *
 *@author Michael Doube
 *@version 0.1
 */

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ImageProcessor;
import ij.plugin.filter.PlugInFilter;
import ij.plugin.frame.RoiManager;
import ij.measure.Calibration;
//import ij.measure.ResultsTable;
import ij.gui.*;

import java.awt.List;
import java.awt.Rectangle;
import java.awt.event.*;

import Jama.*;

import org.doube.bonej.ResultInserter;
/*
 * TODO incorporate curvature
 * curvature = maximum deflection of the slice centroid from the 
 * SVD long axis as a ratio to bone length.  Deflection should
 * be resolved to M-L and Cr-Ca planes; M-L plane contains SVD
 * long axis and centre of femoral head; Cr-Ca plane is perpendicular to
 * M-L plane
 */

public class Neck_Shaft_Angle implements PlugInFilter, MouseListener{
    ImagePlus imp;
    ImageCanvas canvas;
    protected ImageStack stack;
    public float[] CTable;
    public double[] coeff = {0,1}, neckPoint = {78.75, 100.55, 80}, headCentre, centroid;
    public double[][] shaftVector;
    public int minT = 0, maxT = 4000; //min and maximum bone value in HU
    public int startSlice = 1, endSlice;
    public String title, units, valueUnit;
    public Calibration cal;
    //	public ResultsTable rt;

    public int setup(String arg, ImagePlus imp){
	if (imp == null || imp.getNSlices() < 2){
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
    public void run(ImageProcessor ip){		
	//set up pixel calibration
	if (!this.cal.isSigned16Bit() && !this.cal.calibrated()){
	    IJ.run("Threshold...");
	    new WaitForUserDialog("This image is not density calibrated.\nSet the threshold, then click OK.").show();
	    this.minT = (int)ip.getMinThreshold();
	    this.maxT = (int)ip.getMaxThreshold();
	    showDialog(valueUnit);
	    IJ.log("Image is uncalibrated: using user-determined threshold "+minT+" to "+maxT);
	} else if (this.coeff[0] == -1000 && this.coeff[1] == 1.0) {
	    //looks like an HU calibrated image
	    //convert HU limits to pixel values
	    showDialog("Hounsfield units");
	    this.minT = (int)Math.round(cal.getRawValue(minT));
	    this.maxT = (int)Math.round(cal.getRawValue(maxT));
	    IJ.log("Image looks like it is HU calibrated. Using "+minT+" and "+maxT+" "+valueUnit+" as bone cutoffs");
	} else if (cal.isSigned16Bit() && !cal.calibrated()){
	    new WaitForUserDialog("This image is not density calibrated.\nSet the threshold, then click OK.").show();
	    this.minT = (int)ip.getMinThreshold();
	    this.maxT = (int)ip.getMaxThreshold();
	    showDialog(valueUnit);
	    IJ.log("Image is uncalibrated: using user-determined threshold "+minT+" to "+maxT);
	} else {
	    IJ.error("Unrecognised file type");
	    return;
	}
	//get coordinates from the ROI manager and fit a sphere
	RoiManager roiMan = RoiManager.getInstance();
	if (roiMan == null && imp != null){
	    IJ.run("ROI Manager...");
	    IJ.error("Please populate ROI Manager with point ROIs\n" +
	    "placed on the boundary of the femoral head");
	    return;
	} else {
	    headCentre = fitSphere(this.imp, roiMan);
	}
	ImageWindow win = this.imp.getWindow();
	this.canvas = win.getCanvas();

	//work out the centroid and regression vector of the bone
	this.centroid = findCentroid3D(this.stack, startSlice, endSlice);
	if (this.centroid[0] < 0){
	    IJ.error("Empty Stack","No voxels are a" +
	    "vailable for calculation.\nCheck your ROI and threshold.");
	    return;
	}

	this.shaftVector = Regression3D(this.stack, this.centroid);

	//remove stale MouseListeners
	MouseListener[] l = this.canvas.getMouseListeners();
	for (int n = 0; n < l.length; n++){
	    this.canvas.removeMouseListener(l[n]);
	}
	//add a new MouseListener
	this.canvas.addMouseListener(this);

	new WaitForUserDialog("Click on the middle of the femoral neck.\n" +
		"Neck-shaft angle and out-of-plane skew\n" +
	"will be recorded until you hit \'OK\'").show();
	this.canvas.removeMouseListener(this);
	return;
    }

    /**
     * Fit a sphere to points in the ROI manager
     * 
     * @param imp
     * @param roiMan
     * @return double[4] containing (x,y,z) coordinates of centre and radius of sphere
     */
    public double[] fitSphere(ImagePlus imp, RoiManager roiMan){
	double[] voxDim = {cal.pixelWidth, cal.pixelHeight, cal.pixelDepth};
	double[] sphereDim = new double[4];
	int no_p = roiMan.getCount();
	List listRoi = roiMan.getList();
	double[] datapoints = new double[3*no_p];
	double xSum = 0, ySum = 0, zSum = 0;
	Roi[] roiList = roiMan.getRoisAsArray();
	int j = 0;
	for (int i = 0; i<roiMan.getCount(); i++){
	    Roi roi = roiList[i];
	    if (roi.getType() == 10){
		String label = listRoi.getItem(i);
		Rectangle xy = roi.getBoundingRect();
		datapoints[3*j] = xy.getX()*voxDim[0];
		datapoints[3*j+1] = xy.getY()*voxDim[1];
		datapoints[3*j+2] = roiMan.getSliceNumber(label)*voxDim[2];
		xSum += datapoints[3*j];
		ySum += datapoints[3*j+1];
		zSum += datapoints[3*j+2];
		j++;
	    }
	}
	no_p = j;
	if (no_p < 5) IJ.error("ROI Manager contains < 5 point ROIs./n" +
	"Add some point ROIs and try again.");
	sphereDim[0] = xSum/no_p;
	sphereDim[1] = ySum/no_p;
	sphereDim[2] = zSum/no_p;
	double[] r = new double[no_p];
	double g_new = 100.0;
	double g_old = 1.0;
	sphereDim[3] = 0;
	for (int i = 0; i < no_p; i++)
	    sphereDim[3] += Math.sqrt((datapoints[3*i]-sphereDim[0])*(datapoints[3*i]-sphereDim[0])
		    + (datapoints[3*i+1]-sphereDim[1])*(datapoints[3*i+1]-sphereDim[1])
		    + (datapoints[3*i+2]-sphereDim[2])*(datapoints[3*i+2]-sphereDim[2]));
	sphereDim[3] /= no_p;
	while (Math.abs(g_new-g_old) > 1e-10){
	    Matrix J = new Matrix(no_p, 4);
	    double[][] Jp = J.getArray();
	    Matrix d = new Matrix(no_p, 1);
	    double[][] dp = d.getArray(); //dp is a pointer to d's values
	    g_old = g_new;
	    for (int i = 0; i < no_p; i++){
		r[i] = Math.sqrt((datapoints[3*i]-sphereDim[0])*(datapoints[3*i]-sphereDim[0])
			+ (datapoints[3*i+1]-sphereDim[1])*(datapoints[3*i+1]-sphereDim[1])
			+ (datapoints[3*i+2]-sphereDim[2])*(datapoints[3*i+2]-sphereDim[2]));
		dp[i][0] = r[i] - sphereDim[3];
		Jp[i][0] = -(datapoints[3*i]-sphereDim[0])/r[i];
		Jp[i][1] = -(datapoints[3*i+1]-sphereDim[1])/r[i];
		Jp[i][2] = -(datapoints[3*i+2]-sphereDim[2])/r[i];
		Jp[i][3] = -1;
	    }
	    d = d.times(-1);
	    Matrix J1 = J;
	    J = J.transpose();
	    Matrix J2 = J.times(J1);
	    Matrix Jd = J.times(d);
	    Matrix x = J2.inverse().times(Jd);
	    double[][] xp = x.getArray();
	    sphereDim[0] += xp[0][0];
	    sphereDim[1] += xp[1][0];
	    sphereDim[2] += xp[2][0];
	    sphereDim[3] += xp[3][0];
	    d = d.times(-1);
	    Matrix G = J.times(d);
	    double[][] Gp = G.getArray();
	    g_new = 0.0;
	    for (int i = 0; i < 4; i++)
		g_new += Gp[i][0];
	}
	IJ.log("headCentre = ("+sphereDim[0]+", "+sphereDim[1]+", "+sphereDim[2]+"): headRadius = "+sphereDim[3]+" "+cal.getUnits());
	return sphereDim;
    }/* end fitSPhere */

    /**
     * Find the centroid of the voxels of interest
     * 
     * @param stack
     * @param startSlice
     * @param endSlice
     * @param rt
     * @return centroid double[3]
     */
    public double[] findCentroid3D(ImageStack stack, int startSlice, int endSlice){
	Rectangle r = stack.getRoi();
	int offset, i;
	int w = stack.getWidth();
	int h = stack.getHeight();
	int sliceSize = w*h;    	
	double sumx = 0; double sumy = 0; double sumz = 0; double count = 0;
	for (int z = startSlice; z <= endSlice; z++){
	    IJ.showStatus("Calculating centroid...");
	    IJ.showProgress(z, endSlice);
	    short[] slicePixels = new short[sliceSize];
	    slicePixels = (short[])stack.getPixels(z);
	    for (int y = r.y; y < (r.y+r.height); y++) {
		offset = y*w;
		for (int x=r.x; x<(r.x+r.width); x++) {
		    i = offset + x;
		    int testPixel  = slicePixels[i]&0xffff;
		    if (testPixel >= minT && testPixel <= maxT){
			sumx += (double)x;
			sumy += (double)y;
			sumz += (double)z;
			count++;
		    }
		}
	    }
	}
	//centroid in pixels
	double centX = sumx / count;
	double centY = sumy / count;
	double centZ = sumz / count;
	if (sumx == 0){
	    centX = -1; centY = -1; centZ = -1;
	}
	//centroid in real units
	double[] centroid = {centX * cal.pixelWidth, centY * cal.pixelHeight, centZ * cal.pixelDepth};
	return centroid;
    }/*end findCentroid3D */

    /**
     * Calculate the vector associated with the projection plane
     * from the regression vector and the vector connecting the
     * centroid and the femoral head centre
     * 
     * @param shaftVector double[][]
     * @param headCentre double[][]
     * @param centroid double[][]
     * @return double[][] projectionPlane
     */
    public double[][] projectionPlane(double[][] shaftVector, double[] headCentre, double[] centroid){
	//have to calculate distance between points so that we find a unit vector
	double d = Math.sqrt(
		(headCentre[0] - centroid[0]) * (headCentre[0] - centroid[0]) +
		(headCentre[1] - centroid[1]) * (headCentre[1] - centroid[1]) +
		(headCentre[2] - centroid[2]) * (headCentre[2] - centroid[2])
	);
	double[][] cHVec = new double[3][1];
	cHVec[0][0] = (headCentre[0] - centroid[0])/d;
	cHVec[1][0] = (headCentre[1] - centroid[1])/d;
	cHVec[2][0] = (headCentre[2] - centroid[2])/d;

	Matrix cH = new Matrix(cHVec);
	printMatrix(cH, "cHVec");

	//projectionPlane is the cross product of cHVec and shaftVector
	double[][] projectionPlane = new double[3][1];
	projectionPlane[0][0] = cHVec[1][0] * shaftVector[2][0] - cHVec[2][0] * shaftVector[1][0];
	projectionPlane[1][0] = cHVec[2][0] * shaftVector[0][0] - cHVec[0][0] * shaftVector[2][0];
	projectionPlane[2][0] = cHVec[0][0] * shaftVector[1][0] - cHVec[1][0] * shaftVector[0][0];

	d = Math.sqrt(
		projectionPlane[0][0] * projectionPlane[0][0] +
		projectionPlane[1][0] * projectionPlane[1][0] +
		projectionPlane[2][0] * projectionPlane[2][0]
	);
	projectionPlane[0][0] /= d;
	projectionPlane[1][0] /= d;
	projectionPlane[2][0] /= d;

	return projectionPlane;
    }

    /**
     * Calculate the vector associated with the plane formed
     * between neckVector and the normal to projectionPlane
     * from the regression vector and the vector connecting the
     * centroid and the femoral head centre
     * 
     * @param projectionPlane double[][]
     * @param neckVector double[][]
     * @return double[][] neckPlane
     */
    public double[][] neckPlane(double[][] projectionPlane, double[][] neckVector){
	//neckPlane is the cross product of neckVector and projectionPlane
	double[][] neckPlane = new double[3][1];
	neckPlane[0][0] = projectionPlane[1][0] * neckVector[2][0] - projectionPlane[2][0] * neckVector[1][0];
	neckPlane[1][0] = projectionPlane[2][0] * neckVector[0][0] - projectionPlane[0][0] * neckVector[2][0];
	neckPlane[2][0] = projectionPlane[0][0] * neckVector[1][0] - projectionPlane[1][0] * neckVector[0][0];
	return neckPlane;
    }

    /**
     * Find the intersection between neckPlane and projectionPlane 
     * @param projectionPlane
     * @param neckPlane
     * @return double[][] testVector.  
     */
    public double[][] testVector(double[][] projectionPlane, double[][] neckPlane){
	//testVector is the cross product of neckPlane and projectionPlane
	double[][] testVector = new double[3][1];
	testVector[0][0] = projectionPlane[1][0] * neckPlane[2][0] - projectionPlane[2][0] * neckPlane[1][0];
	testVector[1][0] = projectionPlane[2][0] * neckPlane[0][0] - projectionPlane[0][0] * neckPlane[2][0];
	testVector[2][0] = projectionPlane[0][0] * neckPlane[1][0] - projectionPlane[1][0] * neckPlane[0][0];
	return testVector;
    }/* end testVector */

    /**
     * Calculate the orthogonal distance regression plane 
     * of a set of points by
     * the covariance method and Singular Value Decomposition
     * 
     * @param stack
     * @param centroid
     * @return SingularValueDecomposition containing eigenvector and eigenvalue
     * 
     * @see <a href="http://mathforum.org/library/drmath/view/63765.html">Description on Ask Dr Math</a>
     * 
     */
    public double[][] Regression3D(ImageStack stack, double[] centroid){
	IJ.showStatus("Calculating SVD");
	Rectangle r = stack.getRoi();
	IJ.log("Rectangle r has top left coordinates of ("+r.x+", "+r.y+") and size of "+r.width+" x "+r.height);
	int offset, i;
	int w = stack.getWidth();
	int h = stack.getHeight();
	int sliceSize = w*h;
	double vW = cal.pixelWidth;
	double vH = cal.pixelHeight;
	double vD = cal.pixelDepth;
	double[][] C = new double[3][3];
	double count = 0;
	for (int z = startSlice; z <= endSlice; z++){
	    IJ.showStatus("Calculating centroid...");
	    IJ.showProgress(z, endSlice);
	    short[] slicePixels = new short[sliceSize];
	    slicePixels = (short[])stack.getPixels(z);
	    for (int y = r.y; y < (r.y+r.height); y++) {
		offset = y*w;
		for (int x=r.x; x<(r.x+r.width); x++) {
		    i = offset + x;
		    int testPixel = slicePixels[i]&0xffff;
		    if (testPixel >= minT && testPixel <= maxT){
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
	double invCount = 1/count;
	Matrix covarianceMatrix = new Matrix(C).times(invCount);
	printMatrix(covarianceMatrix, "Covariance matrix");
	SingularValueDecomposition S = new SingularValueDecomposition(covarianceMatrix);
	Matrix leftVectors = S.getU();
	printMatrix(leftVectors, "Left vectors");
	double[][] orthogonalDistanceRegression = new double[3][1];
	orthogonalDistanceRegression[0][0] = leftVectors.get(0, 0);
	orthogonalDistanceRegression[1][0] = leftVectors.get(1, 0);
	orthogonalDistanceRegression[2][0] = leftVectors.get(2, 0);
	return orthogonalDistanceRegression;
    }/* end Regression3D */

    public double[][] neckVector(double[] headCentre, double[] neckPoint){		
	//have to calculate d to make sure that neckVector is a unit vector
	double d = Math.sqrt(
		(headCentre[0] - neckPoint[0]) * (headCentre[0] - neckPoint[0]) +
		(headCentre[1] - neckPoint[1]) * (headCentre[1] - neckPoint[1]) +
		(headCentre[2] - neckPoint[2]) * (headCentre[2] - neckPoint[2])
	);

	double[][] neckVector = new double[3][1];
	neckVector[0][0] = (headCentre[0] - neckPoint[0])/d;
	neckVector[1][0] = (headCentre[1] - neckPoint[1])/d;
	neckVector[2][0] = (headCentre[2] - neckPoint[2])/d;
	return neckVector;
    }

    public void calculateAngles(){
	double[][] neckVector = neckVector(headCentre, neckPoint); //contains a 3D coordinate that defines the middle of the neck; perhaps get it with a mouse click and dynamically update neck-shaft angle?
	double[][] projectionPlane = projectionPlane(shaftVector, headCentre, centroid);
	double[][] neckPlane = neckPlane(neckVector, projectionPlane);
	double[][] testVector = testVector(projectionPlane, neckPlane);
	//	P . Q = ||P|| ||Q|| cos(a)  so if P and Q are unit vectors, then P.Q = cos(a)
	Matrix pP = new Matrix(projectionPlane);
	printMatrix(pP, "projectionPlane");

	Matrix tV = new Matrix(testVector);
	printMatrix(tV, "testVector");

	Matrix sV = new Matrix(shaftVector);
	printMatrix(sV, "shaftVector");

	Matrix nV = new Matrix(neckVector);
	printMatrix(nV, "neckVector");

	double cosA1 = sV.get(0, 0)*tV.get(0, 0) + sV.get(1, 0)*tV.get(1, 0) + sV.get(2, 0)*tV.get(2, 0);
	//printMatrix(cosA1, "cosA1");
	IJ.log("cosA1: "+cosA1);

	double cosA2 = nV.get(0, 0)*tV.get(0, 0) + nV.get(1, 0)*tV.get(1, 0) + nV.get(2, 0)*tV.get(2, 0);
	//printMatrix(cosA2, "cosA2");
	IJ.log("cosA2: "+cosA2);

	double neckShaftAngle = Math.acos(cosA1);
	double neckShaftSkew = Math.acos(cosA2);
/*	rt = ResultsTable.getResultsTable();
	rt.incrementCounter();
	rt.addLabel("Label", imp.getTitle());
	rt.addValue("Angle (rad)", neckShaftAngle);  //angle between shaft and neck in plane of head and shaft
	rt.addValue("Skew (rad)", neckShaftSkew);     //angle bewteen neck and plane of head and shaft
	rt.show("Results");*/
	ResultInserter ri = new ResultInserter();
	ri.setResultInRow(this.imp, "Angle (rad)", neckShaftAngle);
	ri.setResultInRow(this.imp, "Skew (rad)", neckShaftSkew);
    }

    public void mousePressed(MouseEvent e) {
	int x = canvas.offScreenX(e.getX());
	int y = canvas.offScreenY(e.getY());
	int z = imp.getCurrentSlice();
	neckPoint[0] = x * cal.pixelWidth;
	neckPoint[1] = y * cal.pixelHeight;
	neckPoint[2] = z * cal.pixelDepth;
	IJ.log("neckPoint: ("+neckPoint[0]+","+neckPoint[1]+", "+neckPoint[2]+")");
	calculateAngles();
    }

    public void mouseReleased(MouseEvent e) {}	
    public void mouseExited(MouseEvent e) {}
    public void mouseClicked(MouseEvent e) {}	
    public void mouseEntered(MouseEvent e) {}
    public void mouseMoved(MouseEvent e) {}


    public void showDialog(String valueUnit){
	GenericDialog gd = new GenericDialog("Setup");
	gd.addNumericField("Start Slice:",startSlice,0);
	gd.addNumericField("End Slice:",endSlice,0);
	gd.addMessage("Only use pixels between clip values:");
	gd.addNumericField("Clip min.", minT, 0, 6, valueUnit);
	gd.addNumericField("Clip max.", maxT, 0, 6, valueUnit);
	gd.showDialog();
	if (gd.wasCanceled()) {
	    return;
	}
	startSlice = (int)gd.getNextNumber();
	endSlice = (int)gd.getNextNumber();
	minT = (int)gd.getNextNumber();
	maxT = (int)gd.getNextNumber();
    }

    public void printMatrix(Matrix matrix, String title){
	IJ.log(title);
	int nCols = matrix.getColumnDimension();
	int nRows = matrix.getRowDimension();
	double[][] eVal = matrix.getArrayCopy();
	for (int r = 0; r < nRows; r++){
	    String row = "||";
	    for (int c = 0; c < nCols; c++){
		row = row + eVal[r][c] + "|";
	    }
	    row = row + "|";
	    IJ.log(row);
	}
	return;
    }
}
