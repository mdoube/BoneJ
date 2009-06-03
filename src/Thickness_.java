
import ij.*;
import ij.measure.ResultsTable;
import ij.plugin.filter.PlugInFilter;
import ij.process.*;
import java.util.Arrays;
import org.doube.bonej.ResultInserter;

/* Bob Dougherty 8/10/2007
Perform all of the steps for the local thickness calculaton


 License:
	Copyright (c) 2007, OptiNav, Inc.
	All rights reserved.

	Redistribution and use in source and binary forms, with or without
	modification, are permitted provided that the following conditions
	are met:

		Redistributions of source code must retain the above copyright
	notice, this list of conditions and the following disclaimer.
		Redistributions in binary form must reproduce the above copyright
	notice, this list of conditions and the following disclaimer in the
	documentation and/or other materials provided with the distribution.
		Neither the name of OptiNav, Inc. nor the names of its contributors
	may be used to endorse or promote products derived from this software
	without specific prior written permission.

	THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
	"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
	LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
	A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
	CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
	EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
	PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
	PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
	LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
	NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
	SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */
/**
 * @see <p>Hildebrand T, Rüegsegger P (1997) A new method for the model-independent assessment of thickness in three-dimensional images.
 * J Microsc 185: 67-75. <a href="http://dx.doi.org/10.1046/j.1365-2818.1997.1340694.x">doi:10.1046/j.1365-2818.1997.1340694.x</a></p>
 * 
 * <p>Saito T, Toriwaki J (1994) New algorithms for euclidean distance transformation of an n-dimensional digitized picture with applications.
 * Pattern Recognit 27: 1551-1565. <a href="http://dx.doi.org/10.1016/0031-3203(94)90133-3">doi:10.1016/0031-3203(94)90133-3</a></p>
 * 
 * @author Bob Dougherty
 * @author Michael Doube (refactoring for BoneJ)
 * 
 */
public class Thickness_ implements  PlugInFilter {
    private ImagePlus baseImp;
    public int thresh = 128;
    public boolean inverse;
    public byte[][] data;
    public float[][] sNew;//, s;
    public int w,h,d;
    public double vW, vH, vD;

    public int setup(String arg, ImagePlus imp) {
	if (imp == null){
	    IJ.noImage();
	    return DONE;
	}
	if (imp.getType() == ImagePlus.GRAY8 || imp.getType() == ImagePlus.COLOR_256) {
	    ImageStatistics stats = imp.getStatistics();
	    if (stats.histogram[0]+stats.histogram[255]!=stats.pixelCount) {
		IJ.error("8-bit binary (black and white only) image required.");
		return DONE;
	    }
	} else {
	    IJ.error("8-bit binary (black and white only) image required.");
	    return DONE;
	}
	checkVoxelSize(imp);
	
	this.baseImp = imp;
	return DOES_8G;
    }
    public void run(ImageProcessor ip) {
	String title = stripExtension(baseImp.getTitle());
	baseImp.unlock();
	IJ.freeMemory();
	w = baseImp.getWidth();
	h = baseImp.getHeight();
	d = baseImp.getStackSize();
	vW = baseImp.getCalibration().pixelWidth;
	vH = baseImp.getCalibration().pixelHeight;
	vD = baseImp.getCalibration().pixelDepth;
	//calculate trabecular thickness (Tb.Th)
	float[][] s = GeometrytoDistanceMap(baseImp); //8-bit in, 32-bit out
	DistanceMaptoDistanceRidge(s); //32-bit in, 32-bit out
	DistanceRidgetoLocalThickness(s); //32-bit in, 32-bit out
	ImagePlus impLTC = LocalThicknesstoCleanedUpLocalThickness(s); //32-bit in, 32-bit out
	IJ.freeMemory();

	impLTC.setTitle(title+"_Tr.Th");
	impLTC.setCalibration(baseImp.getCalibration());
	meanStdDev(impLTC);

	impLTC.show();
	IJ.run("Fire");

	// check marrow cavity size (i.e. trabcular separation, Tb.Sp)
	inverse = true;
	s = GeometrytoDistanceMap(baseImp); //8-bit in, 32-bit out
	DistanceMaptoDistanceRidge(s); //32-bit in, 32-bit out
	DistanceRidgetoLocalThickness(s); //32-bit in, 32-bit out
	ImagePlus impLTCi = LocalThicknesstoCleanedUpLocalThickness(s); //32-bit in, 32-bit out
	for (int i = 0; i < s.length; i++){
	    Arrays.fill(s[i], 0);
	}
	IJ.freeMemory();
	impLTCi.setTitle(title+"_Tb.Sp");
	impLTCi.setCalibration(baseImp.getCalibration());
	impLTCi.show();
	IJ.run("Fire");

	meanStdDev(impLTCi);
	IJ.showProgress(1.0);
	IJ.showStatus("Done");
    }
    //Modified from ImageJ code by Wayne Rasband
    String stripExtension(String name) {
	if (name!=null) {
	    int dotIndex = name.lastIndexOf(".");
	    if (dotIndex>=0)
		name = name.substring(0, dotIndex);
	}
	return name;
    }

    /**
     * <p>Saito-Toriwaki algorithm for Euclidian Distance Transformation.
     * Direct application of Algorithm 1.  Bob Dougherty 8/8/2006</p>
     *
     *<ul>
     * <li>Version S1A: lower memory usage.</li>
     * <li>Version S1A.1 A fixed indexing bug for 666-bin data set</li>
     * <li>Version S1A.2 Aug. 9, 2006.  Changed noResult value.</li>
     * <li>Version S1B Aug. 9, 2006.  Faster.</li>
     * <li>Version S1B.1 Sept. 6, 2006.  Changed comments.</li>
     * <li>Version S1C Oct. 1, 2006.  Option for inverse case.
     * <br />Fixed inverse behavior in y and z directions.</li>
     * <li>Version D July 30, 2007.  Multithread processing for step 2.</li>
     * </ul>
     *
     *<p>This version assumes the input stack is already in memory, 8-bit, and
     *outputs to a new 32-bit stack.  Versions that are more stingy with memory
     *may be forthcoming.</p> 
     *
     * @param imp 8-bit (binary) ImagePlus
     * 
     */
    private float[][] GeometrytoDistanceMap(ImagePlus imp){
	int nThreads = Runtime.getRuntime().availableProcessors();

	//Create references to input data
	ImageStack stack = imp.getStack();
	byte[][] data = new byte[d][];
	for (int k = 0; k < d; k++) data[k] = (byte[])stack.getPixels(k+1);
	
	//Create 32 bit floating point stack for output, s.  Will also use it for g in Transformation 1.
	float[][] s = new float[d][];
	for(int k = 0; k < d; k++){
	    ImageProcessor ipk = new FloatProcessor(w,h);
	    s[k] = (float[])ipk.getPixels();
	}
	float[] sk;
	//Transformation 1.  Use s to store g.
	IJ.showStatus("EDT transformation 1/3");
	Step1Thread[] s1t = new Step1Thread[nThreads];
	for(int thread = 0; thread < nThreads; thread++){
	    s1t[thread] = new Step1Thread(thread,nThreads,w,h,d,thresh,s,data);
	    s1t[thread].start();
	}
	try{
	    for(int thread = 0; thread< nThreads; thread++){
		s1t[thread].join();
	    }
	}catch(InterruptedException ie){
	    IJ.error("A thread was interrupted in step 1 .");
	}		
	//Transformation 2.  g (in s) -> h (in s)
	IJ.showStatus("EDT transformation 2/3");
	Step2Thread[] s2t = new Step2Thread[nThreads];
	for(int thread = 0; thread < nThreads; thread++){
	    s2t[thread] = new Step2Thread(thread,nThreads,w,h,d,s);
	    s2t[thread].start();
	}
	try{
	    for(int thread = 0; thread< nThreads; thread++){
		s2t[thread].join();
	    }
	}catch(InterruptedException ie){
	    IJ.error("A thread was interrupted in step 2 .");
	}
	//Transformation 3. h (in s) -> s
	IJ.showStatus("EDT transformation 3/3");
	Step3Thread[] s3t = new Step3Thread[nThreads];
	for(int thread = 0; thread < nThreads; thread++){
	    s3t[thread] = new Step3Thread(thread,nThreads,w,h,d,s,data);
	    s3t[thread].start();
	}
	try{
	    for(int thread = 0; thread< nThreads; thread++){
		s3t[thread].join();
	    }
	}catch(InterruptedException ie){
	    IJ.error("A thread was interrupted in step 3 .");
	}		
	//Find the largest distance for scaling
	//Also fill in the background values.
	float distMax = 0;
	int wh = w*h;
	float dist;
	for(int k = 0; k < d; k++){
	    sk = s[k];
	    for(int ind = 0; ind < wh; ind++){
		if(((data[k][ind]&255) < thresh)^inverse){
		    sk[ind] = 0;
		}else{
		    dist = (float)Math.sqrt(sk[ind]);
		    sk[ind] = dist;
		    distMax = (dist > distMax) ? dist : distMax;
		}
	    }
	}
	IJ.showProgress(1.0);
	IJ.showStatus("Done");
	return s;
    }
    class Step1Thread extends Thread{
	int thread,nThreads,w,h,d,thresh;
	float[][] s; 
	byte[][] data;
	public Step1Thread(int thread, int nThreads, int w, int h, int d, int thresh, float[][] s,  byte[][] data){
	    this.thread = thread;
	    this.nThreads = nThreads;
	    this.w = w;
	    this.h = h;
	    this.d = d;
	    this.thresh = thresh;
	    this.data = data;
	    this.s = s;
	}
	public void run(){
	    float[] sk;
	    byte[] dk;
	    int n = w;
	    if(h > n) n = h;
	    if(d > n) n = d;
	    int noResult = 3*(n+1)*(n+1);
	    boolean[] background = new boolean[n];			
//	    boolean nonempty;
	    int test, min;			
	    for(int k = thread; k < d; k+=nThreads){
		IJ.showProgress(k/(1.*d));
		sk = s[k];
		dk = data[k];
		for(int j = 0; j < h; j++){
		    for (int i = 0; i < w; i++){
			background[i] = ((dk[i+w*j]&255) < thresh)^inverse;
		    }
		    for (int i = 0; i < w; i++){
			min = noResult;
			for (int x = i; x < w; x++){
			    if(background[x]){
				test = i - x;
				test *= test;
				min = test;
				break;
			    }
			}
			for (int x = i-1; x >=0 ; x--){
			    if(background[x]){
				test = i - x;
				test *= test;
				if(test < min)min = test;
				break;
			    }
			}
			sk[i+w*j] = min;
		    }
		}
	    }
	}//run
    }//Step1Thread
    class Step2Thread extends Thread{
	int thread,nThreads,w,h,d;
	float[][] s; 
	public Step2Thread(int thread, int nThreads, int w, int h, int d, float[][] s){
	    this.thread = thread;
	    this.nThreads = nThreads;
	    this.w = w;
	    this.h = h;
	    this.d = d;
	    this.s = s;
	}
	public void run(){
	    float[] sk;
	    int n = w;
	    if(h > n) n = h;
	    if(d > n) n = d;
	    int noResult = 3*(n+1)*(n+1);
	    int[] tempInt = new int[n];
	    int[] tempS = new int[n];
	    boolean nonempty;
	    int test, min, delta;
	    for(int k = thread; k < d; k+=nThreads){
		IJ.showProgress(k/(1.*d));
		sk = s[k];
		for (int i = 0; i < w; i++){
		    nonempty = false;
		    for (int j = 0; j < h; j++){
			tempS[j] = (int)sk[i+w*j];
			if(tempS[j] >0)nonempty = true;
		    }
		    if(nonempty){
			for (int j = 0; j < h; j++){
			    min = noResult;
			    delta = j;
			    for(int y = 0; y < h; y++){
				test = tempS[y] + delta*delta--;
				if(test < min)min = test;
			    }
			    tempInt[j] = min;
			}
			for (int j = 0; j < h; j++){
			    sk[i+w*j] = tempInt[j];
			}
		    }
		}
	    }
	}//run
    }//Step2Thread	
    class Step3Thread extends Thread{
	int thread,nThreads,w,h,d;
	float[][] s; 
	byte[][] data;
	public Step3Thread(int thread, int nThreads, int w, int h, int d, float[][] s, byte[][] data){
	    this.thread = thread;
	    this.nThreads = nThreads;
	    this.w = w;
	    this.h = h;
	    this.d = d;
	    this.s = s;
	    this.data = data;
	}
	public void run(){
	    int zStart,zStop,zBegin,zEnd;
//	    float[] sk;
	    int n = w;
	    if(h > n) n = h;
	    if(d > n) n = d;
	    int noResult = 3*(n+1)*(n+1);
	    int[] tempInt = new int[n];
	    int[] tempS = new int[n];
	    boolean nonempty;
	    int test, min, delta;
	    for(int j = thread; j < h; j+=nThreads){
		IJ.showProgress(j/(1.*h));
		for(int i = 0; i < w; i++){
		    nonempty = false;
		    for(int k = 0; k < d; k++){
			tempS[k] = (int)s[k][i+w*j];
			if(tempS[k] >0)nonempty = true;
		    }
		    if(nonempty){
			zStart = 0;
			while((zStart < (d-1))&&(tempS[zStart] == 0))zStart++;
			if(zStart > 0)zStart--;
			zStop = d-1;
			while((zStop > 0)&&(tempS[zStop] == 0))zStop--;
			if(zStop < (d-1))zStop++;

			for(int k = 0; k < d; k++){
			    //Limit to the non-background to save time,
			    if(((data[k][i+w*j]&255) >= thresh)^inverse){
				min = noResult;
				zBegin = zStart;
				zEnd = zStop;
				if(zBegin > k)zBegin = k;
				if(zEnd < k)zEnd = k;
				delta = k - zBegin;
				for (int z = zBegin; z <= zEnd; z++){
				    test = tempS[z] + delta*delta--;
				    if(test < min)min = test;
				    //min = (test < min) ? test : min;
				}
				tempInt[k] = min;
			    }
			}
			for(int k = 0; k < d; k++){
			    s[k][i+w*j] = tempInt[k];
			}
		    }
		}
	    }
	}
    }

    /**
     * <p>DistanceMaptoDistanceRidge</p>
     * <p>Output: Distance ridge resulting from a local scan of the distance map.  Overwrites the input.</p>
     * <p>Note: Non-background points that are not part of the distance ridge are assiged a VERY_SMALL_VALUE.
     * This is used for subsequent processing by other plugins to find the local thickness. Bob Dougherty August 10, 2006</p>
     * 
     * <ul>
     * <li>Version 1: August 10-11, 2006.  Subtracts 0.5 from the distances.</li>
     * <li>Version 1.01: September 6, 2006.  Corrected some typos in the comments.</li>
     * <li>Version 1.01: Sept. 7, 2006.  More tiny edits.</li>
     * <li>Version 2: Sept. 25, 2006.  Creates a separate image stack for symmetry.
     * <br />Temporary version that is very conservative.
     * <br />Admittedly does not produce much impovement on real images.</li>
     * <li>Version 3: Sept. 30, 2006.  Ball calculations based on grid points.  Should be much more accurate.</li>
     * <li>Version 3.1 Oct. 1, 2006.  Faster scanning of search points.</li>
     * </ul>
     *
     * @param imp 3D Distance map (32-bit stack)
     */
    private void DistanceMaptoDistanceRidge(float[][] s){
	sNew = new float[d][];
	for(int k = 0; k < d; k++){
	    ImageProcessor ipk = new FloatProcessor(w,h);
	    sNew[k] = (float[])ipk.getPixels();
	}

	//Do it
	int k1,j1,i1,dz,dy,dx;
	boolean notRidgePoint;
	float[] sk1;
	float[] sk, skNew;
	int sk0Sq,sk0SqInd,sk1Sq;
	//Find the largest distance in the data
	IJ.showStatus("Distance Ridge: scanning the data");
	float distMax = 0;
	for (int k = 0; k < d; k++){
	    sk = s[k];
	    for (int j = 0; j < h; j++){
		for (int i = 0; i < w; i++){
		    int ind = i + w*j;
		    if(sk[ind] > distMax)distMax = sk[ind];
		}
	    }
	}
	int rSqMax = (int)(distMax*distMax + 0.5f)+1;
	boolean[] occurs = new boolean[rSqMax];
	for(int i = 0; i < rSqMax; i++)occurs[i] = false;
	for (int k = 0; k < d; k++){
	    sk = s[k];
	    for (int j = 0; j < h; j++){
		for (int i = 0; i < w; i++){
		    int ind = i + w*j;
		    occurs[(int)(sk[ind]*sk[ind] + 0.5f)] = true;
		}
	    }
	}
	int numRadii = 0;
	for (int i = 0; i < rSqMax; i++){
	    if(occurs[i])numRadii++;
	}
	//Make an index of the distance-squared values
	int[] distSqIndex = new int[rSqMax];
	int[] distSqValues = new int[numRadii];
	int indDS = 0;
	for (int i = 0; i < rSqMax; i++){
	    if(occurs[i]){
		distSqIndex[i] = indDS;
		distSqValues[indDS++] = i;
	    }
	}
	//Build template
	//The first index of the template is the number of nonzero components
	//in the offest from the test point to the remote point.  The second
	//index is the radii index (of the test point).  The value of the template
	//is the minimum square radius of the remote point required to cover the
	//ball of the test point.
	IJ.showStatus("Distance Ridge: creating search templates");
	int[][] rSqTemplate = createTemplate(distSqValues);
	int numCompZ,numCompY,numCompX,numComp;
	for (int k = 0; k < d; k++){
	    IJ.showStatus("Distance Ridge: processing slice "+(k+1)+"/"+d);
	    //IJ.showProgress(k/(1.*d));
	    sk = s[k];
	    skNew = sNew[k];
	    for (int j = 0; j < h; j++){
		for (int i = 0; i < w; i++){
		    int ind = i + w*j;
		    if(sk[ind] > 0){
			notRidgePoint = false;
			sk0Sq = (int)(sk[ind]*sk[ind] + 0.5f);
			sk0SqInd = distSqIndex[sk0Sq];
			for (dz = -1; dz <= 1; dz++){
			    k1 = k + dz;
			    if((k1 >= 0)&&(k1 < d)){
				sk1 = s[k1];
				if(dz == 0){
				    numCompZ = 0;
				}else{
				    numCompZ = 1;
				}
				for (dy = -1; dy <= 1; dy++){
				    j1 = j + dy;
				    if((j1 >= 0)&&(j1 < h)){
					if(dy == 0){
					    numCompY = 0;
					}else{
					    numCompY = 1;
					}
					for (dx = -1; dx <= 1; dx++){
					    i1 = i + dx;
					    if((i1 >= 0)&&(i1 < w)){
						if(dx == 0){
						    numCompX = 0;
						}else{
						    numCompX = 1;
						}
						numComp = numCompX + numCompY + numCompZ;
						if(numComp > 0){
						    sk1Sq = (int)(sk1[i1+w*j1]*sk1[i1+w*j1] + 0.5f);
						    if(sk1Sq >= rSqTemplate[numComp-1][sk0SqInd])
							notRidgePoint = true;
						}
					    }//if in grid for i1
					    if(notRidgePoint)break;
					}//dx
				    }//if in grid for j1
				    if(notRidgePoint)break;
				}//dy
			    }//if in grid for k1
			    if(notRidgePoint)break;
			}//dz
			if(!notRidgePoint)skNew[ind] = sk[ind];
		    }//if not in background
		}//i
	    }//j
	}//k
	IJ.showStatus("Distance Ridge complete");
	//replace work array s with result of the method, sNew
	s = sNew;
    }
    //For each offset from the origin, (dx,dy,dz), and each radius-squared,
    //rSq, find the smallest radius-squared, r1Squared, such that a ball
    //of radius r1 centered at (dx,dy,dz) includes a ball of radius
    //rSq centered at the origin.  These balls refer to a 3D integer grid.
    //The set of (dx,dy,dz) points considered is a cube center at the origin.
    //The size of the computed array could be considerably reduced by symmetry,
    //but then the time for the calculation using this array would increase
    //(and more code would be needed).
    int[][] createTemplate(int[] distSqValues){
	int[][] t = new int[3][];
	t[0] = scanCube(1,0,0,distSqValues);
	t[1] = scanCube(1,1,0,distSqValues);
	t[2] = scanCube(1,1,1,distSqValues);
	return t;
    }
    //For a list of r^2 values, find the smallest r1^2 values such
    //that a "ball" of radius r1 centered at (dx,dy,dz) includes a "ball"
    //of radius r centered at the origin.  "Ball" refers to a 3D integer grid.
    int[] scanCube(int dx, int dy, int dz,int[] distSqValues){
	int numRadii = distSqValues.length;
	int[] r1Sq = new int[numRadii];
	if((dx==0)&&(dy==0)&&(dz==0)){
	    for(int rSq = 0; rSq < numRadii; rSq++){
		r1Sq[rSq] = Integer.MAX_VALUE;
	    }
	}else{
	    int dxAbs = -(int)Math.abs(dx);
	    int dyAbs = -(int)Math.abs(dy);
	    int dzAbs = -(int)Math.abs(dz);
	    for(int rSqInd = 0; rSqInd < numRadii; rSqInd++){
		int rSq = distSqValues[rSqInd];
		int max = 0;
		int r = 1 + (int)Math.sqrt(rSq);
		int scank,scankj;
		int dk,dkji;
//		int iBall;
		int iPlus;
		for(int k = 0; k <= r; k++){
		    scank = k*k;
		    dk = (k-dzAbs)*(k-dzAbs);
		    for (int j = 0; j <= r; j++){
			scankj = scank + j*j;
			if(scankj <= rSq){
			    iPlus = ((int)Math.sqrt(rSq - scankj)) - dxAbs;
			    dkji = dk + (j-dyAbs)*(j-dyAbs) + iPlus*iPlus;
			    if(dkji > max) max = dkji;
			}
		    }
		}
		r1Sq[rSqInd] = max;
	    }
	}
	return r1Sq;
    }	

    /**
     * <p>DistanceRidgetoLocalThickness</p>
     * <p>Input: Distance Ridge (32-bit stack) (Output from Distance Ridge.java)
     * Output: Local Thickness.  Overwrites the input.</p>
     * <ul>
     * <li>Version 1: September 6, 2006.</li>
     * <li>Version 2: September 25, 2006.  Fixed several bugs that resulted in 
     * non-symmetrical output from symmetrical input.</li>
     * <li>Version 2.1 Oct. 1, 2006.  Fixed a rounding error that caused some points to be missed.</li>
     * <li>Version 3 July 31, 2007.  Parallel processing version.</li>
     * <li>Version 3.1  Multiplies the output by 2 to conform with the definition of local thickness</li>
     * </ul>
     * @param imp
     */
    private void DistanceRidgetoLocalThickness(float[][] s){
	float[] sk;
	//Count the distance ridge points on each slice
	int[] nRidge = new int[d];
	int ind, nr, iR;
	IJ.showStatus("Local Thickness: scanning stack ");
	for (int k = 0; k < d; k++){
	    sk = s[k];
	    nr = 0;
	    for (int j = 0; j < h; j++){
		for (int i = 0; i < w; i++){
		    ind = i + w*j;
		    if(sk[ind] > 0) nr++;
		}
	    }
	    nRidge[k] = nr;
	}
	int[][] iRidge = new int[d][];
	int[][] jRidge = new int[d][];
	float[][] rRidge = new float[d][];
	//Pull out the distance ridge points
	int[] iRidgeK,jRidgeK;
	float[] rRidgeK;
	float sMax = 0;
	for (int k = 0; k < d; k++){
	    nr = nRidge[k];
	    iRidge[k] = new int[nr];
	    jRidge[k] = new int[nr];
	    rRidge[k] = new float[nr];
	    sk = s[k];
	    iRidgeK = iRidge[k];
	    jRidgeK = jRidge[k];
	    rRidgeK = rRidge[k];
	    iR = 0;
	    for (int j = 0; j < h; j++){
		for (int i = 0; i < w; i++){
		    ind = i + w*j;
		    if(sk[ind] > 0){;
		    iRidgeK[iR] = i;
		    jRidgeK[iR] = j;
		    rRidgeK[iR++] = sk[ind];
		    if(sk[ind]>sMax)sMax = sk[ind];
		    sk[ind] = 0;
		    }
		}
	    }
	}
	int nThreads = Runtime.getRuntime().availableProcessors();
	final Object[] resources = new Object[d];//For synchronization
	for(int k = 0; k < d; k++){
	    resources[k] = new Object();
	}
	LTThread[] ltt = new LTThread[nThreads];
	for(int thread = 0; thread < nThreads; thread++){
	    ltt[thread] = new LTThread(thread,nThreads,w,h,d,nRidge,
		    s,iRidge,jRidge,rRidge,resources);
	    ltt[thread].start();
	}
	try{
	    for(int thread = 0; thread< nThreads; thread++){
		ltt[thread].join();
	    }
	}catch(InterruptedException ie){
	    IJ.error("A thread was interrupted .");
	}		

	//Fix the square values and apply factor of 2
	IJ.showStatus("Local Thickness: square root ");
	for (int k = 0; k < d; k++){
	    sk = s[k];
	    for (int j = 0; j < h; j++){
		for (int i = 0; i < w; i++){
		    ind = i + w*j;
		    sk[ind] = (float)(2 * Math.sqrt(sk[ind]));
		}
	    }
	}
	IJ.showStatus("Local Thickness complete");
	return;
    }
    class LTThread extends Thread{
	int thread,nThreads,w,h,d,nR;
	float[][] s;
	int[] nRidge;
	int[][] iRidge, jRidge;
	float[][] rRidge;
	Object[] resources;
	public LTThread(int thread, int nThreads, int w, int h, int d, int[] nRidge,
		float[][] s, int[][] iRidge, int[][] jRidge, float[][] rRidge,
		Object[] resources){
	    this.thread = thread;
	    this.nThreads = nThreads;
	    this.w = w;
	    this.h = h;
	    this.d = d;
	    this.s = s;
	    this.nRidge = nRidge;
	    this.iRidge = iRidge;
	    this.jRidge = jRidge;
	    this.rRidge = rRidge;
	    this.resources = resources;
	}
	public void run(){
	    int i,j;
	    float[] sk1;//sk,sk1;
	    //Loop through ridge points.  For each one, update the local thickness for
	    //the points within its sphere.
	    float r;
	    int rInt,ind1;
	    int iStart,iStop,jStart,jStop,kStart,kStop;
	    float r1SquaredK,r1SquaredJK, r1Squared,s1;
	    int rSquared;
	    int[] iRidgeK,jRidgeK;
	    float[] rRidgeK;
	    for(int k = thread; k < d; k+=nThreads){
		IJ.showStatus("Local Thickness: processing slice "+(k+1)+"/"+d);
		int nR = nRidge[k];
		iRidgeK = iRidge[k];
		jRidgeK = jRidge[k];
		rRidgeK = rRidge[k];
		//sk = s[k];
		for (int iR = 0; iR < nR; iR++){
		    i = iRidgeK[iR];
		    j = jRidgeK[iR];
		    r = rRidgeK[iR];
		    rSquared = (int)(r*r + 0.5f);
		    rInt = (int)r;
		    if(rInt < r)rInt++;
		    iStart = i - rInt;
		    if(iStart < 0)iStart = 0;
		    iStop = i + rInt;
		    if(iStop >= w) iStop = w-1;
		    jStart = j - rInt;
		    if(jStart < 0)jStart = 0;
		    jStop = j + rInt;
		    if(jStop >= h) jStop = h-1;
		    kStart = k - rInt;
		    if(kStart < 0)kStart = 0;
		    kStop = k + rInt;
		    if(kStop >= d) kStop = d-1;
		    for(int k1 = kStart; k1 <= kStop; k1++){
			r1SquaredK = (k1 - k)*(k1 - k);
			sk1 = s[k1];
			for(int j1 = jStart; j1 <= jStop; j1++){
			    r1SquaredJK = r1SquaredK + (j1 - j)*(j1 - j);
			    if(r1SquaredJK <= rSquared){
				for(int i1 = iStart; i1 <= iStop; i1++){
				    r1Squared = r1SquaredJK + (i1 - i)*(i1 - i);
				    if(r1Squared <= rSquared){
					ind1 = i1 + w*j1;
					s1 = sk1[ind1];
					if(rSquared > s1){
					    //Get a lock on sk1 and check again to make sure
					    //that another thread has not increased
					    //sk1[ind1] to something larger than rSquared.
					    //A test shows that this may not be required...
					    synchronized(resources[k1]){
						s1 = sk1[ind1];
						if(rSquared > s1){
						    sk1[ind1] = rSquared;
						}
					    }
					}
				    }//if within shere of DR point
				}//i1
			    }//if k and j components within sphere of DR point
			}//j1
		    }//k1
		}//iR
	    }//k
	}//run
    }//LTThread

    /**
     * <p>LocalThicknesstoCleanedUpLocalThickness</p>
     *
     * <p>Input: 3D Local Thickness map (32-bit stack)</p>
     * <p>Output: Same as input with border voxels corrected for "jaggies." Non-background voxels
     * adjacent to background voxels are have their local thickness values replaced by the average of
     * their non-background neighbors that do not border background points.  Bob Dougherty August 1, 2007</p>
     * 
     * <ul>
     * <li>August 10.  Version 3 This version also multiplies the local thickness by 2 
     * to conform with the official definition of local thickness.</li>
     * </ul>
     * 
     */
    private ImagePlus LocalThicknesstoCleanedUpLocalThickness(float[][] s){
	IJ.showStatus("Cleaning up local thickness...");
	//Create 32 bit floating point stack for output, sNew.
	ImageStack newStack = new ImageStack(w,h);
	sNew = new float[d][];
	for(int k = 0; k < d; k++){
	    ImageProcessor ipk = new FloatProcessor(w,h);
	    newStack.addSlice(null,ipk);
	    sNew[k] = (float[])ipk.getPixels();
	}
	//First set the output array to flags:
	// 0 for a background point
	// -1 for a non-background point that borders a background point
	// s (input data) for an interior non-background point
	for (int k = 0; k < d; k++){
	    for (int j = 0; j < h; j++){
		for (int i = 0; i < w; i++){
		    sNew[k][i + w*j] = setFlag(s, i,j,k);
		}//i
	    }//j
	}//k
	//Process the surface points.  Initially set results to negative values
	//to be able to avoid including them in averages of for subsequent points.
	//During the calculation, positive values in sNew are interior non-background
	//local thicknesses.  Negative values are surface points.  In this case the
	//value might be -1 (not processed yet) or -result, where result is the
	//average of the neighboring interior points.  Negative values are excluded from
	//the averaging.
	for (int k = 0; k < d; k++){
	    for (int j = 0; j < h; j++){
		for (int i = 0; i < w; i++){
		    int ind = i + w*j;
		    if(sNew[k][ind] == -1){
			sNew[k][ind] = -averageInteriorNeighbors(s, i,j,k);
		    }
		}//i
	    }//j
	}//k
	//Fix the negative values and double the results
	for (int k = 0; k < d; k++){
	    for (int j = 0; j < h; j++){
		for (int i = 0; i < w; i++){
		    int ind = i + w*j;
		    sNew[k][ind] = (float)Math.abs(sNew[k][ind]);
		}//i
	    }//j
	}//k
	IJ.showStatus("Clean Up Local Thickness complete");
	String title = stripExtension(baseImp.getTitle());
	ImagePlus impOut = new ImagePlus(title+"_CL",newStack);
	double vW = baseImp.getCalibration().pixelWidth;
	//calibrate the pixel values to pixel width
	//so that thicknesses represent real units (not pixels)
	for (int z = 0; z < d; z++){
	    impOut.setSlice(z+1);
	    impOut.getProcessor().multiply(vW);
	}
	impOut.getProcessor().setMinAndMax(0, 1.5 * impOut.getProcessor().getMax());
	return impOut;
    }
    float setFlag(float[][] s, int i,int j,int k){
	//!!! null pointer here
	if(s[k][i+w*j]==0)return 0;
	//change 1
	if(look(s, i,j,k-1)==0)return -1;
	if(look(s, i,j,k+1)==0)return -1;
	if(look(s, i,j-1,k)==0)return -1;
	if(look(s, i,j+1,k)==0)return -1;
	if(look(s, i-1,j,k)==0)return -1;
	if(look(s, i+1,j,k)==0)return -1;
	//change 1 before plus
	if(look(s, i,j+1,k-1)==0)return -1;
	if(look(s, i,j+1,k+1)==0)return -1;
	if(look(s, i+1,j-1,k)==0)return -1;
	if(look(s, i+1,j+1,k)==0)return -1;
	if(look(s, i-1,j,k+1)==0)return -1;
	if(look(s, i+1,j,k+1)==0)return -1;
	//change 1 before minus
	if(look(s, i,j-1,k-1)==0)return -1;
	if(look(s, i,j-1,k+1)==0)return -1;
	if(look(s, i-1,j-1,k)==0)return -1;
	if(look(s, i-1,j+1,k)==0)return -1;
	if(look(s, i-1,j,k-1)==0)return -1;
	if(look(s, i+1,j,k-1)==0)return -1;
	//change 3, k+1
	if(look(s, i+1,j+1,k+1)==0)return -1;
	if(look(s, i+1,j-1,k+1)==0)return -1;
	if(look(s, i-1,j+1,k+1)==0)return -1;
	if(look(s, i-1,j-1,k+1)==0)return -1;
	//change 3, k-1
	if(look(s, i+1,j+1,k-1)==0)return -1;
	if(look(s, i+1,j-1,k-1)==0)return -1;
	if(look(s, i-1,j+1,k-1)==0)return -1;
	if(look(s, i-1,j-1,k-1)==0)return -1;
	return s[k][i+w*j];
    }
    float averageInteriorNeighbors(float[][] s, int i,int j,int k){
	int n = 0;
	float sum = 0;
	//change 1
	float value = lookNew(i,j,k-1);
	if(value > 0){
	    n++;
	    sum += value;
	}
	value = lookNew(i,j,k+1);
	if(value > 0){
	    n++;
	    sum += value;
	}
	value = lookNew(i,j-1,k);
	if(value > 0){
	    n++;
	    sum += value;
	}
	value = lookNew(i,j+1,k);
	if(value > 0){
	    n++;
	    sum += value;
	}
	value = lookNew(i-1,j,k);
	if(value > 0){
	    n++;
	    sum += value;
	}
	value = lookNew(i+1,j,k);
	if(value > 0){
	    n++;
	    sum += value;
	}
	//change 1 before plus
	value = lookNew(i,j+1,k-1);
	if(value > 0){
	    n++;
	    sum += value;
	}
	value = lookNew(i,j+1,k+1);
	if(value > 0){
	    n++;
	    sum += value;
	}
	value = lookNew(i+1,j-1,k);
	if(value > 0){
	    n++;
	    sum += value;
	}
	value = lookNew(i+1,j+1,k);
	if(value > 0){
	    n++;
	    sum += value;
	}
	value = lookNew(i-1,j,k+1);
	if(value > 0){
	    n++;
	    sum += value;
	}
	value = lookNew(i+1,j,k+1);
	if(value > 0){
	    n++;
	    sum += value;
	}
	//change 1 before minus
	value = lookNew(i,j-1,k-1);
	if(value > 0){
	    n++;
	    sum += value;
	}
	value = lookNew(i,j-1,k+1);
	if(value > 0){
	    n++;
	    sum += value;
	}
	value = lookNew(i-1,j-1,k);
	if(value > 0){
	    n++;
	    sum += value;
	}
	value = lookNew(i-1,j+1,k);
	if(value > 0){
	    n++;
	    sum += value;
	}
	value = lookNew(i-1,j,k-1);
	if(value > 0){
	    n++;
	    sum += value;
	}
	value = lookNew(i+1,j,k-1);
	if(value > 0){
	    n++;
	    sum += value;
	}
	//change 3, k+1
	value = lookNew(i+1,j+1,k+1);
	if(value > 0){
	    n++;
	    sum += value;
	}
	value = lookNew(i+1,j-1,k+1);
	if(value > 0){
	    n++;
	    sum += value;
	}
	value = lookNew(i-1,j+1,k+1);
	if(value > 0){
	    n++;
	    sum += value;
	}
	value = lookNew(i-1,j-1,k+1);
	if(value > 0){
	    n++;
	    sum += value;
	}
	//change 3, k-1
	value = lookNew(i+1,j+1,k-1);
	if(value > 0){
	    n++;
	    sum += value;
	}
	value = lookNew(i+1,j-1,k-1);
	if(value > 0){
	    n++;
	    sum += value;
	}
	value = lookNew(i-1,j+1,k-1);
	if(value > 0){
	    n++;
	    sum += value;
	}
	value = lookNew(i-1,j-1,k-1);
	if(value > 0){
	    n++;
	    sum += value;
	}
	if(n > 0)return sum/n;
	return s[k][i+w*j];
    }
    float look(float[][] s, int i,int j,int k){
	if((i < 0)||(i >= w))return -1;
	if((j < 0)||(j >= h))return -1;
	if((k < 0)||(k >= d))return -1;
	return s[k][i+w*j];
    }
    //A positive result means this is an interior, non-background, point.
    float lookNew(int i,int j,int k){
	if((i < 0)||(i >= w))return -1;
	if((j < 0)||(j >= h))return -1;
	if((k < 0)||(k >= d))return -1;
	return  sNew[k][i+w*j];
    }
    /**
     * Work out some summary stats
     * 
     * @param imp 32-bit thickness image
     */
    private void meanStdDev(ImagePlus imp){
	String units = imp.getCalibration().getUnits();

	ImageStack stack = imp.getStack();
	long pixCount = 0;
	double sumThick = 0;
	long pixCountFiltered = 0;
	double sumThickFiltered = 0;
	double maxThick = 0;
	double minRes = 2 * baseImp.getCalibration().pixelWidth;
	for (int s = 1; s <= stack.getSize(); s++){
	    float[] slicePixels = (float[])stack.getPixels(s);
	    for (int p = 0; p < slicePixels.length; p++){
		double pixVal = slicePixels[p];
		if (pixVal > 0){
		    sumThick += pixVal;
		    maxThick = Math.max(maxThick, pixVal);
		    pixCount++;
		    if (pixVal > minRes){
			sumThickFiltered += pixVal;
			pixCountFiltered++;
		    }
		}
	    }
	}
	double meanThick = sumThick / pixCount;
	double meanThickFiltered = sumThickFiltered / pixCountFiltered;

	double sumSquares = 0;
	double sumSquaresFiltered = 0;
	for (int s = 1; s <= stack.getSize(); s++){
	    float[] slicePixels = (float[])stack.getPixels(s);
	    for (int p = 0; p < slicePixels.length; p++){
		double pixVal = slicePixels[p];
		if (pixVal > 0){
		    double d = meanThick - pixVal;
		    sumSquares += d * d;
		    if (pixVal > minRes){
			double df = meanThickFiltered - pixVal;
			sumSquaresFiltered += df * df;					    
		    }
		}
	    }
	}
	double stDev = Math.sqrt(sumSquares / pixCount);
	double stDevF = Math.sqrt(sumSquaresFiltered / pixCountFiltered);
	
/*	ResultsTable rt = ResultsTable.getResultsTable();
	if (!inverse){ 
	    rt.incrementCounter();
	    //trab thickness
	    rt.addLabel("Label", stripExtension(baseImp.getTitle()));
	    rt.addValue("Tb.Th Mean ("+units+")", meanThick);
	    rt.addValue("Tb.Th Mean F ("+units+")", meanThickFiltered);
	    rt.addValue("Tb.Th Std Dev ("+units+")", stDev);
	    rt.addValue("Tb.Th Std Dev F ("+units+")", stDevF);
	    rt.addValue("Tb.Th Max ("+units+")", maxThick);
	} else {
	    //trab separation
	    rt.addValue("Tb.Sp Mean ("+units+")", meanThick);
	    rt.addValue("Tb.Sp Mean F ("+units+")", meanThickFiltered);
	    rt.addValue("Tb.Sp Std Dev ("+units+")", stDev);
	    rt.addValue("Tb.Sp Std Dev F ("+units+")", stDevF);
	    rt.addValue("Tb.Sp Max ("+units+")", maxThick);
	}
	rt.show("Results");
*/	
	ResultInserter ri = new ResultInserter();
	if (!inverse){ 
	    //trab thickness
	    ri.setResultInRow(baseImp, "Tb.Th Mean ("+units+")", meanThick);
	    ri.setResultInRow(baseImp, "Tb.Th Mean F ("+units+")", meanThickFiltered);
	    ri.setResultInRow(baseImp, "Tb.Th Std Dev ("+units+")", stDev);
	    ri.setResultInRow(baseImp, "Tb.Th Std Dev F ("+units+")", stDevF);
	    ri.setResultInRow(baseImp, "Tb.Th Max ("+units+")", maxThick);
	} else {
	    //trab separation
	    ri.setResultInRow(baseImp, "Tb.Sp Mean ("+units+")", meanThick);
	    ri.setResultInRow(baseImp, "Tb.Sp Mean F ("+units+")", meanThickFiltered);
	    ri.setResultInRow(baseImp, "Tb.Sp Std Dev ("+units+")", stDev);
	    ri.setResultInRow(baseImp, "Tb.Sp Std Dev F ("+units+")", stDevF);
	    ri.setResultInRow(baseImp, "Tb.Sp Max ("+units+")", maxThick);
	}	
	return;
    }
    
    private void checkVoxelSize(ImagePlus imp){
	if(imp.getCalibration().pixelDepth != imp.getCalibration().pixelWidth || 
		imp.getCalibration().pixelDepth != imp.getCalibration().pixelHeight ||
		imp.getCalibration().pixelHeight != imp.getCalibration().pixelWidth){
	    IJ.showMessage("Voxels are anisotropic, please take care with results.");
	}
    }
}