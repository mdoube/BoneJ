

/**
 * Slice_Geometry plugin for ImageJ
 * Copyright 2009 Michael Doube 
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;
import ij.plugin.filter.PlugInFilter;
import ij.measure.Calibration;
import ij.measure.ResultsTable;
import ij.gui.*;

import java.awt.Rectangle;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * <p>Calculate 2D geometrical parameters</p>
 * 
 * @author Michael Doube
 *
 */

public class Slice_Geometry implements PlugInFilter {
    ImagePlus imp;
    protected ImageStack stack;
    public static final double PI = 3.141592653589793;
    private int boneID, al;
    private double vW, vH, vD, airHU, minBoneHU, maxBoneHU;
    private String units, analyse, calString;
    private boolean doThickness, doCentroids, doCopy, doOutline, doAxes, doStack, isCalibrated;
    private double[] cslice, Sx, Sy, Sxx, Syy, Sxy, Myy, Mxx, Mxy, theta, 
    Imax, Imin, Ipm, R1, R2, maxRadMin, maxRadMax, Zmax, Zmin, ImaxFast, IminFast;
    private boolean[] emptySlices;
    private double[][] sliceCentroids;
    Calibration cal;

    public int setup(String arg, ImagePlus imp) {
	if (imp == null) {
	    IJ.noImage();
	    return DONE;
	}
	this.imp = imp;
	this.cal = imp.getCalibration();
	this.stack = imp.getStack();
	this.al = this.stack.getSize()+1;
	//TODO properly support 8bit images
	return DOES_8G + DOES_16 + SUPPORTS_MASKING;
    }

    public void run(ImageProcessor ip) {
	//show a setup dialog and set calibration
	setHUCalibration();

	if (!showDialog()){
	    return;
	}
		
	calculateCentroids();
	calculateMoments();
	
	showSliceResults();
	
    }

    private void calculateCentroids(){
	//2D centroids
	this.sliceCentroids = new double[2][this.al]; 
	//pixel counters
	double cstack = 0;
	Rectangle r = this.stack.getRoi();
	int w = this.stack.getWidth();
	this.emptySlices = new boolean[this.al];
	this.cslice = new double[this.al];
	double[] cortArea;
	cortArea = new double[this.al];
	for (int s = 1; s <= this.stack.getSize(); s++) {
	    double sumx = 0; double sumy = 0;
	    this.cslice[s] = 0;
	    short[] pixels = (short[])this.stack.getPixels(s);
	    for (int y=r.y; y<(r.y+r.height); y++) {
		int offset = y*w;
		for (int x=r.x; x<(r.x+r.width); x++) {
		    int i = offset + x;
		    if (pixels[i] >= this.minBoneHU && pixels[i] <= this.maxBoneHU){
			this.cslice[s]++;
			cortArea[s] += this.vW * this.vH;
			sumx += x * this.vW;
			sumy += y * this.vH;
		    }
		}
	    }
	    if (this.cslice[s] > 0){
		this.sliceCentroids[0][s] = sumx / this.cslice[s];
		this.sliceCentroids[1][s] = sumy / this.cslice[s];
		cstack += this.cslice[s];
		this.emptySlices[s] = false;
	    } else {
		this.emptySlices[s] = true;
	    }
	}
	if (cstack == 0){
	    IJ.error("Empty Stack","No pixels are available for calculation.\nCheck your ROI and threshold.");
	    return;
	}
    }

    private void calculateMoments(){
	Rectangle r = this.stack.getRoi();
	int w = this.stack.getWidth();
	//START OF Ix AND Iy CALCULATION 
	this.Sx = new double[this.al]; this.Sy = new double[this.al]; this.Sxx = new double[this.al]; this.Syy = new double[this.al]; this.Sxy = new double[this.al];
	this.Myy = new double[this.al]; this.Mxx = new double[this.al]; this.Mxy = new double[this.al]; this.theta = new double[this.al];
	for (int s = 1; s <= this.stack.getSize(); s++) {
	    if (!this.emptySlices[s]){
		short[] pixels = (short[])this.stack.getPixels(s);
		for (int y = r.y; y < (r.y + r.height); y++) {
		    int offset = y * w;
		    for (int x = r.x; x < (r.x+r.width); x++) {
			int i = offset + x;
			if (pixels[i] >= this.minBoneHU && pixels[i] <= this.maxBoneHU){
			    this.Sx[s] += x;			
			    this.Sy[s] += y;			
			    this.Sxx[s] += x*x;
			    this.Syy[s] += y*y;
			    this.Sxy[s] += y*x;
			}
		    }
		}
		this.Myy[s] = this.Sxx[s] - (this.Sx[s] * this.Sx[s] / this.cslice[s]) + this.cslice[s]/12; //this.cslice[]/12 is for each pixel's moment around its own centroid
		this.Mxx[s] = this.Syy[s] - (this.Sy[s] * this.Sy[s] / this.cslice[s]) + this.cslice[s]/12;
		this.Mxy[s] = this.Sxy[s] - (this.Sx[s] * this.Sy[s] / this.cslice[s]) + this.cslice[s]/12;
		if (this.Mxy[s] == 0) this.theta[s] = 0;
		else {
		    this.theta[s] = Math.atan((this.Mxx[s]-this.Myy[s]+Math.sqrt(Math.pow(this.Mxx[s]-this.Myy[s],2) + 4*this.Mxy[s]*this.Mxy[s]))/(2*this.Mxy[s]));
		    //thetaFast gives same result except jumps when hits PI/4 and -PI/4
		    //thetaFast[s] = Math.atan(2*Mxy[s]/(Myy[s]-Mxx[s])) / 2;
		}
	    }
	}
	//END OF Ix and Iy CALCULATION
	//START OF Imax AND Imin CALCULATION
	this.Imax = new double[this.al]; this.Imin = new double[this.al]; this.Ipm = new double[this.al]; 
	this.R1 = new double[this.al]; this.R2 = new double[this.al]; this.maxRadMin = new double[this.al]; this.maxRadMax = new double[this.al];
	this.Zmax = new double[this.al]; this.Zmin = new double[this.al]; this.ImaxFast = new double[this.al]; this.IminFast = new double[this.al]; 
	for (int s=1;s<=this.stack.getSize();s++) {
	    if (!this.emptySlices[s]){
		short[] pixels = (short[])this.stack.getPixels(s);
		this.Sx[s]=0; this.Sy[s]=0; this.Sxx[s]=0; this.Syy[s]=0; this.Sxy[s]=0;
		for (int y=r.y; y<(r.y+r.height); y++) {
		    int offset = y*w;
		    for (int x=r.x; x<(r.x+r.width); x++) {
			int i = offset + x;
			if (pixels[i] >= this.minBoneHU && pixels[i] <= this.maxBoneHU){
			    this.Sx[s] += x*Math.cos(this.theta[s]) + y*Math.sin(this.theta[s]);				//normal distance from parallel axis summed over pixels
			    this.Sy[s] += y*Math.cos(this.theta[s]) - x*Math.sin(this.theta[s]);
			    this.Sxx[s] += (x*Math.cos(this.theta[s]) + y*Math.sin(this.theta[s]))
			    	* (x*Math.cos(this.theta[s]) + y*Math.sin(this.theta[s]));				//squared normal distances from parallel axis (Iz)
			    this.Syy[s] += (y*Math.cos(this.theta[s]) - x*Math.sin(this.theta[s]))
			    	* (y*Math.cos(this.theta[s]) - x*Math.sin(this.theta[s]));
			    this.Sxy[s] += (y*Math.cos(theta[s]) - x*Math.sin(theta[s]))
			    	*(x*Math.cos(theta[s])+y*Math.sin(theta[s]));          				
			    //maximum distance from minimum principal axis (longer)
			    this.maxRadMin[s] = Math.max(this.maxRadMin[s],Math.abs((x-this.sliceCentroids[0][s])*Math.cos(this.theta[s])
				    + (y-this.sliceCentroids[1][s])*Math.sin(theta[s])));
			    //maximum distance from maximum principal axis (shorter)
			    this.maxRadMax[s] = Math.max(this.maxRadMax[s],Math.abs((y-this.sliceCentroids[1][s])*Math.cos(this.theta[s]) - (x-sliceCentroids[0][s])*Math.sin(theta[s])));
			}
		    }
		}
		this.Imax[s] = this.Sxx[s] - (this.Sx[s] * this.Sx[s] / this.cslice[s]) + this.cslice[s]*(Math.pow(Math.cos(this.theta[s]),2)+Math.pow(Math.sin(this.theta[s]),2)) / 12;			//2nd moment of area around minimum principal axis (shorter axis, bigger I) 
		this.Imin[s] = this.Syy[s] - (this.Sy[s] * this.Sy[s] / this.cslice[s]) + this.cslice[s]*(Math.pow(Math.cos(this.theta[s]),2)+Math.pow(Math.sin(this.theta[s]),2)) / 12;			//2nd moment of area around maximum principal axis (longer axis, smaller I)
		this.Ipm[s] = this.Sxy[s] - (this.Sy[s] * this.Sx[s] / this.cslice[s]) + this.cslice[s]*(Math.pow(Math.cos(this.theta[s]),2)+Math.pow(Math.sin(this.theta[s]),2)) / 12;			//product moment of area, should be 0 if theta calculated perfectly
		this.R1[s] = Math.sqrt(this.Imin[s] / this.cslice[s]);    				//length of major axis
		this.R2[s] = Math.sqrt(this.Imax[s] / this.cslice[s]);					//length of minor axis
		this.Zmax[s] = this.Imax[s] / this.maxRadMin[s];							//Section modulus around maximum principal axis
		this.Zmin[s] = this.Imin[s] / this.maxRadMax[s];							//Section modulus around minimum principal axis
		this.ImaxFast[s] = (this.Mxx[s]+this.Myy[s])/2 + Math.sqrt(Math.pow(((this.Mxx[s]-this.Myy[s])/2),2)+this.Mxy[s]*this.Mxy[s]);
		this.IminFast[s] = (this.Mxx[s]+this.Myy[s])/2 - Math.sqrt(Math.pow(((this.Mxx[s]-this.Myy[s])/2),2)+this.Mxy[s]*this.Mxy[s]);
	    }
	}
    }


    /**
     * Calculate the caliper diameter of an Roi
     * 
     * @param roi Roi to get the caliper diameter of
     * @param rt ResultsTable
     * 
     * @see <p><a href="http://cgm.cs.mcgill.ca/~orm/rotcal.html">Rotating Calipers homepage</a></p>
     * 
     */
    public void rotatingCalipers(Roi roi, ResultsTable rt){
	double dmax = 0;
	double dmin = 0;
	IJ.run("Convex Hull");
	int perimx[] = ((PolygonRoi)roi).getXCoordinates();
	int perimy[] = ((PolygonRoi)roi).getYCoordinates();
	int podal = 0;
	int antipodal = 0;
	//		Compute the polygon's extreme points in the y direction
	int ymin = 999999;
	int ymax = 0;
	int nPoints = perimx.length;
	for (int n = 0; n < nPoints; n++){
	    if (perimy[n] < ymin){
		ymin = perimy[n];
		podal = n;    		//n is always 0 in ImageJ, start point of convex hull-ed selection is minimum x value within set with minimum y
	    }
	    if (perimy[n] > ymax){
		ymax = perimy[n];
		antipodal = n;
	    }
	}

	//		create list of angles beween points
	double angles[]; angles = new double[nPoints];
	for (int n = 0; n < nPoints - 1; n++){
	    angles[n] = Math.atan2(perimy[n]-perimy[n+1], perimx[n]-perimx[n+1]);
	}

	//		angle between last point and first point of selection
	angles[nPoints-1] = Math.atan2(perimy[nPoints-1]-perimy[0], perimx[nPoints-1]-perimx[0]);

	//		correct any angles < -pi
	for (int n=0; n < nPoints; n++) {
	    if (angles[n] < -PI) angles[n] = angles[n] + 2*PI;
	}

	//		Construct two horizontal lines of support through ymin and ymax.
	//		Since this is already an anti-podal pair, compute the distance, and keep as maximum.
	dmax = ymax - ymin;
	dmin = dmax;
	int end = antipodal;
	double startTheta = 0; double endTheta = 0; double testAngle = 0; 
	while (antipodal < nPoints  && podal <= end){
	    int podincr = 0;
	    int antipodincr = 0;
	    //			find the start condition for theta
	    if (podal > 0) {
		startTheta = endTheta;
	    } else {
		if (angles[antipodal - 1] > 0) { 
		    testAngle = angles[antipodal - 1] - PI;
		} else {
		    testAngle = angles[antipodal - 1] + PI;
		}
		if (angles[nPoints -1] <= testAngle){
		    startTheta = angles[nPoints - 1];
		} else {
		    startTheta = testAngle;
		}
	    }

	    //			Rotate the lines until one is flush with an edge of the polygon.
	    //			find the end condition for theta
	    //			increment the podal or antipodal point depending on
	    //			which caliper blade touches the next point first

	    if (angles[antipodal] > 0) { 
		testAngle = angles[antipodal] - PI;
	    } else {
		testAngle = angles[antipodal] + PI;
	    }
	    if (antipodal <= nPoints - 1) {
		if (angles[podal] > testAngle && angles[podal] * testAngle > 0) {
		    endTheta = angles[podal];
		    podincr = 1;
		} else if (angles[podal] < testAngle) {
		    endTheta = testAngle;
		    antipodincr = 1;
		} else if (angles[podal] > testAngle && angles[podal] * testAngle < 0) {
		    endTheta = testAngle;
		    antipodincr = 1;
		} else {
		    endTheta = angles[podal];
		    podincr = 1;
		    antipodincr = 1;
		}
	    } else {
		endTheta = 0;
	    }

	    //			A new anti-podal pair is determined. Compute the new distance,
	    //			compare to old maximum, and update if necessary.
	    double thetaH = Math.atan2(perimy[antipodal] - perimy[podal] , perimx[antipodal] - perimx[podal]);
	    double dp = Math.sqrt(Math.pow(perimx[antipodal]-perimx[podal],2) + Math.pow(perimy[antipodal]-perimy[podal],2));
	    if (dp > dmax) dmax = dp;
	    if (startTheta >= endTheta){
		double ds = dp*Math.abs(Math.cos(-thetaH + startTheta - PI/2));
		double de = dp*Math.abs(Math.cos(-thetaH + endTheta - PI/2));
		if (ds < dmin) dmin = ds;
		else if (de < dmin) dmin = de;
	    }	
	    else if (startTheta < endTheta && startTheta * endTheta < 0){
		double ds = dp*Math.abs(Math.cos(-thetaH + startTheta - PI/2));
		double de = dp*Math.abs(Math.cos(-thetaH + - 3*PI/2));
		if (ds < dmin) dmin = ds;
		else if (de < dmin) dmin = de;
		ds = dp*Math.abs(Math.cos(-thetaH + PI/2));
		de = dp*Math.abs(Math.cos(-thetaH + endTheta - PI/2));
		if (ds < dmin) dmin = ds;
		else if (de < dmin) dmin = de;
	    }
	    //			increment either or both podal or antipodal
	    podal = podal + podincr;
	    antipodal = antipodal + antipodincr;
	}
	double pWidth = cal.pixelWidth;
	String units = cal.getUnits();
	rt.addLabel("Label", imp.getShortTitle());
	rt.addValue("RCmax ("+units+")", dmax*pWidth);
	rt.addValue("RCmin ("+units+")", dmin*pWidth);
    }
//TODO fix this, it's a mess.
    /**
     * Set up threshold based on HU if image is calibrated or
     * user-based threshold if it is uncalibrated
     * 
     */
    private void setHUCalibration(){
	this.minBoneHU = 0;		//minimum bone value in HU 
	this.maxBoneHU = 4000;		//maximum bone value in HU
	this.vW = this.cal.pixelWidth;
	this.vH = this.cal.pixelHeight;
	this.vD = this.cal.pixelDepth;
	this.units = this.cal.getUnits();
	double[] coeff = this.cal.getCoefficients();
	if (!this.cal.calibrated() || this.cal == null || (this.cal.getCValue(0) == 0 && this.cal.getCoefficients()[1] == 1)){
	    this.isCalibrated = false;
	    this.calString = "Image is uncalibrated\nEnter air and bone pixel values";
	    ImageStatistics stats = imp.getStatistics();
	    if (stats.min < 50 && stats.min >= 0){
		this.airHU = 0;
	    }
	    else if (stats.min > 31000){
		this.airHU = 31768;
	    }
	    else if (stats.min < -800){
		this.airHU = -1000;
	    } else {
		this.airHU = 0;
	    }
	} else {
	    this.isCalibrated = true;
	    this.calString = "Image is calibrated\nEnter HU below:";
	    this.airHU = -1000;
	}
	this.minBoneHU = this.airHU + 1000;
	this.maxBoneHU = this.airHU + 5000;
	if (this.cal.calibrated()) {
	    //convert HU limits to pixel values
	    IJ.log("Image is calibrated, using "+this.minBoneHU+" and "+this.maxBoneHU+" HU as bone cutoffs");
	    this.minBoneHU = (short)Math.round(this.cal.getRawValue(this.minBoneHU));
	    this.maxBoneHU = (short)Math.round(this.cal.getRawValue(this.maxBoneHU));
	    IJ.log("Vox Width: "+vW+"; Vox Height: "+vH+" "+units);
	    IJ.log("Calibration coefficients:"+coeff[0]+","+coeff[1]);
	    IJ.log("this.minBoneHU = "+this.minBoneHU+", this.maxBoneHU = "+this.maxBoneHU);
	}
	else {
	    IJ.log("Image is not calibrated, using user-determined threshold");
	    IJ.run("Threshold...");
	    new WaitForUserDialog("This image is not density calibrated.\nSet the threshold, then click OK.").show();
	    this.minBoneHU = (short)this.stack.getProcessor(this.imp.getCurrentSlice()).getMinThreshold();
	    this.maxBoneHU = (short)this.stack.getProcessor(this.imp.getCurrentSlice()).getMaxThreshold();
	}
	return;
    }

    private boolean showDialog(){
	GenericDialog gd = new GenericDialog("Options");

	gd.addCheckbox("Cortical Thickness", true);
	gd.addCheckbox("Draw Axes", true);
	gd.addCheckbox("Draw Centroids", true);
	gd.addCheckbox("Draw Outline", false);
	gd.addCheckbox("Annotated Copy", true);
	gd.addCheckbox("Process Stack", false);
	String[] bones = {"unknown", "scapula", "humerus", "radius", "ulna", "metacarpal", "pelvis", "femur", "tibia", "fibula", "metatarsal"};
	//guess bone from image title
	String title = this.imp.getTitle();
	this.boneID = 0;
	for (int n = 0; n < bones.length; n++){
	    Pattern p = Pattern.compile(bones[n], Pattern.CASE_INSENSITIVE);
	    Matcher m = p.matcher(title);
	    if (m.find()){
		this.boneID = n;
		continue;
	    }
	}
	gd.addChoice("Bone: ", bones, bones[this.boneID]);
	String[] analyses = {"Weighted", "Unweighted", "Both"};
	gd.addChoice("Calculate: ", analyses, analyses[1]);
	gd.addNumericField("Voxel Size (x): ", vW, 3, 8, units);
	gd.addNumericField("Voxel Size (y): ", vH, 3, 8, units);
	gd.addNumericField("Voxel Size (z): ", vD, 3, 8, units);
	gd.addMessage(this.calString);
	gd.addNumericField("Air:", this.airHU, 0);
	gd.addNumericField("Bone Min:", this.minBoneHU, 0);
	gd.addNumericField("Bone Max:", this.maxBoneHU, 0);
	gd.showDialog();
	this.doThickness = gd.getNextBoolean();
	this.doAxes = gd.getNextBoolean();
	this.doCentroids = gd.getNextBoolean();
	this.doOutline = gd.getNextBoolean();
	this.doCopy = gd.getNextBoolean();
	this.doStack = gd.getNextBoolean();
	String bone = gd.getNextChoice();
	for (int n = 0; n < bones.length; n++){
	    if(bone.equals(bones[n])) {
		this.boneID = n;
		continue;
	    }
	}
	this.analyse = gd.getNextChoice();
	this.vW = gd.getNextNumber();
	this.vH = gd.getNextNumber();
	this.vD = gd.getNextNumber();
	this.airHU = gd.getNextNumber();
	this.minBoneHU = gd.getNextNumber();
	this.maxBoneHU = gd.getNextNumber();
	if(gd.wasCanceled()){
	    return false;
	} else {
	    return true;
	}
    }
    
    private void showSliceResults(){
	ResultsTable rt = ResultsTable.getResultsTable();
	rt.reset();

	//TODO fix spatial calibration: this assumes isotropic pixels
	double unit4 = Math.pow(vW, 4);
	double unit3 = Math.pow(vW, 3);
	for (int s = 1; s <= this.stack.getSize(); s++) {
	    rt.incrementCounter();
	    rt.addValue("Slice", s);
	    rt.addValue("X cent. ("+units+")", this.sliceCentroids[0][s]);
	    rt.addValue("Y cent. ("+units+")", this.sliceCentroids[1][s]);
	    rt.addValue("Theta (rad)", theta[s]);
//	    rt.addValue("CA ("+units+"^2)", cortArea[s]);
	    rt.addValue("Imin ("+units+"^4)", Imin[s]*unit4);
	    rt.addValue("IminFast ("+units+"^4)", IminFast[s]*unit4);
	    rt.addValue("Imax ("+units+"^4)", Imax[s]*unit4);
	    rt.addValue("ImaxFast ("+units+"^4)", ImaxFast[s]*unit4);
	    rt.addValue("Ipm ("+units+"^4)", Ipm[s]*unit4);
	    rt.addValue("R1 ("+units+")", R1[s]);
	    rt.addValue("R2 ("+units+")", R2[s]);
	    rt.addValue("Zmax ("+units+"^3)", Zmax[s]*unit3);
	    rt.addValue("Zmin ("+units+"^3)", Zmin[s]*unit3);
	}
	rt.show("Results");
    }
}