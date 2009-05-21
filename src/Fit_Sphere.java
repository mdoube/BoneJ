/**
 * Fit_Sphere plugin for ImageJ
 * Copyright 2008 2009 Michael Doube
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

import Jama.Matrix;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ImageProcessor;
import ij.plugin.filter.PlugInFilter;
import ij.gui.*;
import ij.plugin.frame.*;
import ij.measure.Calibration;
import ij.measure.ResultsTable;

import java.awt.Rectangle;
import java.awt.List;

/**
 *<p>Takes point selections from ROI manager and returns the
 *centroid and radius of a best fit sphere
 *Ported from Angelo Tardugno's C++</p>
 *
 *
 *@author Michael Doube and Angelo Tardugno
 *@version 0.1
 */
public class Fit_Sphere implements PlugInFilter {
    ImagePlus imp;
    ImageProcessor ip;
    RoiManager roiMan = RoiManager.getInstance();
    protected ImageStack sourceStack;
    public boolean doCopy, doInnerCube, doOuterCube;
    public int padding;
    public double cropFactor;

    public int setup(String arg, ImagePlus imp) {
	if (IJ.versionLessThan("1.42h")){
	    IJ.error("Your version of ImageJ is too old/n" +
	    "Please update (Help->Update ImageJ...)");
	    return DONE;
	}
	this.imp = imp;
	if (imp == null || imp.getNSlices() < 2){
	    IJ.showMessage("A stack must be open");
	    return DONE;
	}
	if (roiMan == null && imp != null){
	    IJ.run("ROI Manager...");
	    IJ.error("Please populate ROI Manager with point ROIs");
	    return DONE;
	}
	sourceStack = imp.getStack();
	return DOES_ALL + STACK_REQUIRED;
    }

    public void run(ImageProcessor ip) {		
	double[] voxDim = getVoxDim(imp);
	if (voxDim[2] > voxDim[0]*20){
	    if (!IJ.showMessageWithCancel("Voxel depth problem", "The voxel depth (slice thickness)" +
	    "is unusually large.\nClick OK if you know it is correct.")){
		imp.unlock();
		IJ.run("Properties...");
	    }
	}
	if (!imp.lock()) imp.lock(); //if we have unlocked the image to reset properties, relock it.
	showDialog();
	double[] sphereDim = fitSphere(imp, roiMan);
	if (doCopy) copySphere(imp, ip, padding, cropFactor, sphereDim);
	if (doInnerCube) copyInnerCube(imp, ip, cropFactor, sphereDim);
	if (doOuterCube) copyOuterCube(imp, ip, cropFactor, sphereDim);
    }

    public boolean showDialog(){
	GenericDialog gd = new GenericDialog("Setup");
	gd.addMessage("");
	gd.addCheckbox("Copy Sphere", true);
	gd.addCheckbox("Inner Cube", true);
	gd.addCheckbox("Outer Cube", true);
	gd.addNumericField("Padding", 2, 0, 2, "voxels");
	gd.addNumericField("Crop Factor", 1.0, 2, 4, "");
	gd.showDialog();
	doCopy = gd.getNextBoolean();
	doInnerCube = gd.getNextBoolean();
	doOuterCube = gd.getNextBoolean();
	padding = (int)gd.getNextNumber();
	cropFactor = gd.getNextNumber();
	return true;
    }

    public static double[] getVoxDim(ImagePlus imp){
	Calibration cal = imp.getCalibration();
	double[] voxDim = new double[3];
	voxDim[0] = cal.pixelWidth;
	voxDim[1] = cal.pixelHeight;
	voxDim[2] = cal.pixelDepth;
	return voxDim;
    }
    /**
     * Fit a sphere to a set of 3D point ROI's
     * 
     * @param imp ImagePlus
     * @param roiMan RoiManager containing points placed on sphere
     * @return double[4] containing calibrated (x,y,z,r) coordinates of the centre and radius, r.
     */
    public double[] fitSphere(ImagePlus imp, RoiManager roiMan){
	ResultsTable rt = ResultsTable.getResultsTable();
	int row = rt.getCounter();
	rt.incrementCounter();
	rt.setLabel(imp.getTitle(), row);
	double[] voxDim = getVoxDim(imp);
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
	Calibration cal = imp.getCalibration();
	rt.setValue("X centroid ("+cal.getUnits()+")", row, sphereDim[0]);
	rt.setValue("Y centroid ("+cal.getUnits()+")", row, sphereDim[1]);
	rt.setValue("Z centroid ("+cal.getUnits()+")", row, sphereDim[2]);
	rt.setValue("Radius ("+cal.getUnits()+")", row, sphereDim[3]);
	rt.show("Results");
	return sphereDim;
    }
    //TODO make this go faster by getting slice pixels and iterating through
    //it's array rather than using setSlice and getPixel
    public void copySphere(ImagePlus imp, ImageProcessor ip, int padding, double cropFactor, double[] sphereDim){
	IJ.showStatus("Copying sphere to new stack");
	double[] voxDim = getVoxDim(imp);
	int startX = (int)Math.round((sphereDim[0]-sphereDim[3]*cropFactor)/voxDim[0])-padding;
	int startY = (int)Math.round((sphereDim[1]-sphereDim[3]*cropFactor)/voxDim[1])-padding;
	int startZ = (int)Math.round((sphereDim[2]-sphereDim[3]*cropFactor)/voxDim[2])-padding;
	int roiWidth = (int)Math.round(2*sphereDim[3]*cropFactor/voxDim[0])+2*padding;
	int roiHeight = (int)Math.round(2*sphereDim[3]*cropFactor/voxDim[1])+2*padding;
	int roiDepth = (int)Math.round(2*sphereDim[3]*cropFactor/voxDim[2])+2*padding;
	ImageStack targetStack = new ImageStack(roiWidth,roiHeight);
	for (int z = startZ; z <= startZ+roiDepth; z++){
	    IJ.showProgress(z, startZ+roiDepth);
	    short[] targetSlice = new short[roiWidth*roiHeight];
	    imp.setSlice(z);
	    int nRows = 0;
	    for (int y = startY; y < startY+roiHeight; y++){
		int index = nRows*roiWidth;
		int nCols = 0;
		for (int x = startX; x < startX+roiWidth; x++){
		    double distance = Math.sqrt(
			    (x*voxDim[0] - sphereDim[0])*(x*voxDim[0] - sphereDim[0]) +
			    (y*voxDim[1] - sphereDim[1])*(y*voxDim[1] - sphereDim[1]) +
			    (z*voxDim[2] - sphereDim[2])*(z*voxDim[2] - sphereDim[2])
		    );
		    if (distance < sphereDim[3]*cropFactor){
			targetSlice[index+nCols] = (short)ip.getPixel(x, y);	
		    }
		    else {
			targetSlice[index+nCols] = 0;
		    }
		    nCols++;
		}
		nRows++;
	    }
	    targetStack.addSlice("slice "+z, targetSlice);
	}
	ImagePlus target = new ImagePlus("Sphere", targetStack);
	target.setCalibration(imp.getCalibration());
	target.setDisplayRange(imp.getDisplayRangeMin(), imp.getDisplayRangeMax());
	target.show();		
	return;
    }
    //TODO make this go faster by getting slice pixels and iterating through
    //its array rather than using setSlice and getPixel
    public void copyInnerCube(ImagePlus imp, ImageProcessor ip, double cropFactor, double[] sphereDim){
	Calibration cal = imp.getCalibration();
	IJ.showStatus("Copying largest enclosed cube");
	double[] voxDim = getVoxDim(imp);
	double h = sphereDim[3] * cropFactor / Math.sqrt(3);
	int startX = (int)Math.round((sphereDim[0] - h)/voxDim[0]);
	int startY = (int)Math.round((sphereDim[1] - h)/voxDim[1]);
	int startZ = (int)Math.round((sphereDim[2] - h)/voxDim[2]);
	int roiWidth = (int)Math.round(2 * h / voxDim[0]);
	int roiHeight = (int)Math.round(2 * h / voxDim[1]);
	int roiDepth = (int)Math.round(2 * h / voxDim[2]);
	ImageStack targetStack = new ImageStack(roiWidth,roiHeight);
	for (int z = startZ; z <= startZ+roiDepth; z++){
	    IJ.showProgress(z, startZ+roiDepth);
	    short[] targetSlice = new short[roiWidth*roiHeight];
	    imp.setSlice(z);
	    int nRows = 0;
	    for (int y = startY; y < startY+roiHeight; y++){
		int index = nRows*roiWidth;
		int nCols = 0;
		for (int x = startX; x < startX+roiWidth; x++){
		    targetSlice[index+nCols] = (short)ip.getPixel(x, y);
		    nCols++;
		}
		nRows++;
	    }
	    targetStack.addSlice("slice "+z, targetSlice);
	}
	ImagePlus target = new ImagePlus("Inner Cube", targetStack);
	target.setCalibration(cal);
	target.setDisplayRange(imp.getDisplayRangeMin(), imp.getDisplayRangeMax());
	target.show();		
	return;
    }

    public void copyOuterCube(ImagePlus imp, ImageProcessor ip, double cropFactor, double[] sphereDim){
	Calibration cal = imp.getCalibration();
	IJ.showStatus("Copying smallest enclosing cube");
	double[] voxDim = getVoxDim(imp);
	double h = sphereDim[3] * cropFactor;
	int startX = (int)Math.round((sphereDim[0] - h)/voxDim[0]);
	int startY = (int)Math.round((sphereDim[1] - h)/voxDim[1]);
	int startZ = (int)Math.round((sphereDim[2] - h)/voxDim[2]);
	int roiWidth = (int)Math.round(2 * h / voxDim[0]);
	int roiHeight = (int)Math.round(2 * h / voxDim[1]);
	int roiDepth = (int)Math.round(2 * h / voxDim[2]);
	ImageStack targetStack = new ImageStack(roiWidth,roiHeight);
	for (int z = startZ; z <= startZ+roiDepth; z++){
	    IJ.showProgress(z, startZ+roiDepth);
	    short[] targetSlice = new short[roiWidth*roiHeight];
	    imp.setSlice(z);
	    int nRows = 0;
	    for (int y = startY; y < startY+roiHeight; y++){
		int index = nRows*roiWidth;
		int nCols = 0;
		for (int x = startX; x < startX+roiWidth; x++){
		    targetSlice[index+nCols] = (short)ip.getPixel(x, y);
		    nCols++;
		}
		nRows++;
	    }
	    targetStack.addSlice("slice "+z, targetSlice);
	}
	ImagePlus target = new ImagePlus("Outer Cube", targetStack);
	target.setCalibration(cal);
	target.setDisplayRange(imp.getDisplayRangeMin(), imp.getDisplayRangeMax());
	target.show();	
	return;
    }
}