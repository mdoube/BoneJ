

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
    private int boneID;
    private double vW, vH, vD, airHU, minBoneHU, maxBoneHU;
    private String units, analyse;
    private boolean doThickness, doCentroids, doCopy, doOutline, doAxes, doStack;
    Calibration cal;

    public int setup(String arg, ImagePlus imp) {
	if (imp == null) {
	    IJ.noImage();
	    return DONE;
	}
	this.imp = imp;
	this.cal = imp.getCalibration();
	this.stack = imp.getStack();
	//TODO properly support 8bit images
	return DOES_8G + DOES_16 + SUPPORTS_MASKING;
    }

    public void run(ImageProcessor ip) {
	//TODO split different measurements / calculations into their own methods.

	//show a setup dialog and set calibration
	setHUCalibration();
	if (!showDialog()){
	    return;
	}


	ResultsTable rt = ResultsTable.getResultsTable();
	rt.reset();

	int startSlice = 1;
	int endSlice = stack.getSize();
	int al = stack.getSize()+1;
	//2D centroids
	double[] xc; double[] yc;
	xc = new double[al]; yc = new double[al];
	//pixel counters
	double cstack = 0;
	Rectangle r = ip.getRoi();
	int offset, i;
	int w = stack.getWidth();
	boolean[] emptySlices;
	emptySlices = new boolean[al];
	short[] pixels;
	double[] cslice;
	cslice = new double[al];
	double[] cortArea;
	cortArea = new double[al];
	for (int s = startSlice; s <= endSlice; s++) {
	    double sumx = 0; double sumy = 0;
	    cslice[s] = 0;
	    pixels = (short[])stack.getPixels(s);
	    for (int y=r.y; y<(r.y+r.height); y++) {
		offset = y*w;
		for (int x=r.x; x<(r.x+r.width); x++) {
		    i = offset + x;
		    if (pixels[i] >= this.minBoneHU && pixels[i] <= this.maxBoneHU){
			cslice[s]++;
			cortArea[s] += vW*vH;
			sumx += x*vW;
			sumy += y*vH;
		    }
		}
	    }
	    if (cslice[s] > 0){
		xc[s] = sumx/cslice[s];
		yc[s] = sumy/cslice[s];
		cstack += cslice[s];
		emptySlices[s] = false;
	    } else {
		emptySlices[s] = true;
	    }
	}
	if (cstack == 0){
	    IJ.error("Empty Stack","No pixels are available for calculation.\nCheck your ROI and threshold.");
	    return;
	}
	//END OF CENTROID CALCULATION
	//START OF Ix AND Iy CALCULATION
	//slice centroids are now in the arrays xc[] and yc[] and bone area is held in cslice[]
	double[] Sx; double[] Sy; double[] Sxx; double[] Syy; double[] Sxy; double[] Myy; double[] Mxx; double[] Mxy; double[] theta; 
	Sx = new double[al]; Sy = new double[al]; Sxx = new double[al]; Syy = new double[al]; Sxy = new double[al]; Myy = new double[al]; Mxx = new double[al]; Mxy = new double[al]; theta = new double[al];
	for (int s=startSlice;s<=endSlice;s++) {
	    if (!emptySlices[s]){
		pixels = (short[])stack.getPixels(s);
		for (int y=r.y; y<(r.y+r.height); y++) {
		    offset = y*w;
		    for (int x=r.x; x<(r.x+r.width); x++) {
			i = offset + x;
			if (pixels[i] >= this.minBoneHU && pixels[i] <= this.maxBoneHU){
			    Sx[s] += x;			
			    Sy[s] += y;			
			    Sxx[s] += x*x;
			    Syy[s] += y*y;
			    Sxy[s] += y*x;
			}
		    }
		}
		Myy[s] = Sxx[s] - (Sx[s]*Sx[s] / cslice[s]) + cslice[s]/12; //cslice[]/12 is for each pixel's moment around its own centroid
		Mxx[s] = Syy[s] - (Sy[s]*Sy[s] / cslice[s]) + cslice[s]/12;
		Mxy[s] = Sxy[s] - (Sx[s]*Sy[s] / cslice[s]) + cslice[s]/12;
		if (Mxy[s] == 0){theta[s] = 0;}
		else {
		    theta[s] = Math.atan((Mxx[s]-Myy[s]+Math.sqrt(Math.pow(Mxx[s]-Myy[s],2)+4*Math.pow(Mxy[s],2)))/(2*Mxy[s]));
		    //thetaFast gives same result except jumps when hits PI/4 and -PI/4
		    //thetaFast[s] = Math.atan(2*Mxy[s]/(Myy[s]-Mxx[s])) / 2;
		}
	    }
	}
	//END OF Ix and Iy CALCULATION
	//START OF Imax AND Imin CALCULATION
	double[] Imax; double[] Imin; double[] Ipm; double[] R1; double[] R2; double[] maxRadMin;
	double[] maxRadMax; double[] Zmax; double[] Zmin; double[] ImaxFast; double[] IminFast;
	Imax = new double[al]; Imin = new double[al]; Ipm = new double[al]; R1 = new double[al]; R2 = new double[al]; maxRadMin = new double[al];
	maxRadMax = new double[al]; Zmax = new double[al]; Zmin = new double[al]; ImaxFast = new double[al]; IminFast = new double[al]; 
	for (int s=startSlice;s<=endSlice;s++) {
	    if (!emptySlices[s]){
		pixels = (short[])stack.getPixels(s);
		Sx[s]=0;Sy[s]=0;Sxx[s]=0;Syy[s]=0; Sxy[s]=0;
		for (int y=r.y; y<(r.y+r.height); y++) {
		    offset = y*w;
		    for (int x=r.x; x<(r.x+r.width); x++) {
			i = offset + x;
			if (pixels[i] >= this.minBoneHU && pixels[i] <= this.maxBoneHU){
			    Sx[s] += x*Math.cos(theta[s]) + y*Math.sin(theta[s]);				//normal distance from parallel axis summed over pixels
			    Sy[s] += y*Math.cos(theta[s]) - x*Math.sin(theta[s]);
			    Sxx[s] += Math.pow(x*Math.cos(theta[s]) + y*Math.sin(theta[s]),2);				//squared normal distances from parallel axis (Iz)
			    Syy[s] += Math.pow(y*Math.cos(theta[s]) - x*Math.sin(theta[s]),2);
			    Sxy[s] += (y*Math.cos(theta[s]) - x*Math.sin(theta[s]))*(x*Math.cos(theta[s])+y*Math.sin(theta[s]));          				
			    //maximum distance from minimum principal axis (longer)
			    maxRadMin[s] = Math.max(maxRadMin[s],Math.abs((x-xc[s])*Math.cos(theta[s]) + (y-yc[s])*Math.sin(theta[s])));
			    //maximum distance from maximum principal axis (shorter)
			    maxRadMax[s] = Math.max(maxRadMax[s],Math.abs((y-yc[s])*Math.cos(theta[s]) - (x-xc[s])*Math.sin(theta[s])));
			}
		    }
		}
		Imax[s] = Sxx[s] - (Sx[s] * Sx[s] / cslice[s]) + cslice[s]*(Math.pow(Math.cos(theta[s]),2)+Math.pow(Math.sin(theta[s]),2)) / 12;			//2nd moment of area around minimum principal axis (shorter axis, bigger I) 
		Imin[s] = Syy[s] - (Sy[s] * Sy[s] / cslice[s]) + cslice[s]*(Math.pow(Math.cos(theta[s]),2)+Math.pow(Math.sin(theta[s]),2)) / 12;			//2nd moment of area around maximum principal axis (longer axis, smaller I)
		Ipm[s] = Sxy[s] - (Sy[s] * Sx[s] / cslice[s]) + cslice[s]*(Math.pow(Math.cos(theta[s]),2)+Math.pow(Math.sin(theta[s]),2)) / 12;			//product moment of area, should be 0 if theta calculated perfectly
		R1[s] = Math.sqrt(Imin[s] / cslice[s]);    				//length of major axis
		R2[s] = Math.sqrt(Imax[s] / cslice[s]);					//length of minor axis
		Zmax[s] = Imax[s]/maxRadMin[s];							//Section modulus around maximum principal axis
		Zmin[s] = Imin[s]/maxRadMax[s];							//Section modulus around minimum principal axis
		ImaxFast[s] = (Mxx[s]+Myy[s])/2 + Math.sqrt(Math.pow(((Mxx[s]-Myy[s])/2),2)+Mxy[s]*Mxy[s]);
		IminFast[s] = (Mxx[s]+Myy[s])/2 - Math.sqrt(Math.pow(((Mxx[s]-Myy[s])/2),2)+Mxy[s]*Mxy[s]);
	    }
	}
	//TODO fix spatial calibration: this assumes isotropic pixels
	double unit4 = Math.pow(vW, 4);
	double unit3 = Math.pow(vW, 3);
	for (int s=startSlice; s<=endSlice; s++) {
	    rt.incrementCounter();
	    rt.addValue("Slice", s);
	    rt.addValue("X cent. ("+units+")", xc[s]);
	    rt.addValue("Y cent. ("+units+")", yc[s]);
	    rt.addValue("Theta (rad)", theta[s]);
	    rt.addValue("CA ("+units+"^2)", cortArea[s]);
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
    
    private void setHUCalibration(){
	this.minBoneHU = 0;		//minimum bone value in HU 
	this.maxBoneHU = 4000;		//maximum bone value in HU
	this.vW = this.cal.pixelWidth;
	this.vH = this.cal.pixelHeight;
	this.vD = this.cal.pixelDepth;
	this.units = this.cal.getUnits();
	double[] coeff = this.cal.getCoefficients();
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
	    this.minBoneHU = (short)this.imp.getImageStack().getProcessor(1).getMinThreshold();
	    this.maxBoneHU = (short)this.imp.getImageStack().getProcessor(1).getMaxThreshold();
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
	//	Dialog.addChoice("Units: ", units);
	//	gd.addChoice("Units", units, units);
	//	if (calibrate(0) == 0 && calslope == 1){
	//		var isCalibrated = false;
	//		Dialog.addMessage("Image is uncalibrated\nEnter air and bone pixel values");
	//		getRawStatistics(nPixels, mean, min, max, std, histogram);
	//		if (min < 50 && min >=0) {
	//			airhu = 0;
	//		}
	//		else if (min > 31000) {
	//			airhu = 31768;
	//		}
	//		else if (min < -800){
	//			airhu = -1000;
	//		} else {airhu = 0;}
	//	} else {0.156
	//		var isCalibrated = true;
	//		Dialog.addMessage("Image is calibrated\nEnter HU below:");
	//		airhu = -1000;
	//	}
	this.airHU = -1000;
	this.minBoneHU = this.airHU + 1000;
	this.maxBoneHU = this.airHU + 5000;
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
	//	unit = Dialog.getChoice();
	this.airHU = gd.getNextNumber();
	this.minBoneHU = gd.getNextNumber();
	this.maxBoneHU = gd.getNextNumber();
	//	range = maxhu - airhu;
	if(gd.wasCanceled()){
	    return false;
	} else {
	    return true;
	}
    }
}