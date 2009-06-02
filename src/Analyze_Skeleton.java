import java.util.ArrayList;

import org.doube.bonej.Skeletonize3D_;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.measure.ResultsTable;
import ij.plugin.filter.PlugInFilter;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;
import ij.process.ShortProcessor;

/**
 * Analyze_Skeleton plugin for ImageJ(C).
 * Copyright (C) 2008,2009 Ignacio Arganda-Carreras 
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation (http://www.gnu.org/licenses/gpl.txt )
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 * 
 */
//TODO use List Iterators rather than for-each when modifying List elements
//TODO branch count too high sometimes
//TODO extra trees found sometimes
/**
 * Main class.
 * This class is a plugin for the ImageJ interface for analyzing
 * 2D/3D skeleton images.
 * <p>
 * For more information, visit the Analyze_Skeleton homepage:
 * http://imagejdocu.tudor.lu/doku.php?id=plugin:analysis:analyzeskeleton:start
 *
 *
 * @version 1.0 05/27/2009
 * @author Ignacio Arganda-Carreras <ignacio.arganda@uam.es>
 *
 */
public class Analyze_Skeleton implements PlugInFilter
{
    /** end point flag */
    private static final byte END_POINT = 30;
    /** junction flag */
    private static final byte JUNCTION = 70;
    /** slab flag */
    private static final byte SLAB = 127;
    /** Square root of 3 */
    private static final double SQRT_3 = Math.sqrt(3);
    /** voxel width in real units */
    private double VOXEL_WIDTH;
    /** voxel height in real units */
    private double VOXEL_HEIGHT;
    /** voxel depth in real units */
    private double VOXEL_DEPTH;

    /** working image plus */
    private ImagePlus imRef;

    /** working image width */
    private int width = 0;
    /** working image height */
    private int height = 0;
    /** working image depth */
    private int depth = 0;
    /** working image stack*/
    private ImageStack inputImage = null;

    /** visit flags */
    private boolean [][][] visited = null;

    // Tree fields
    /** number of branches for every specific tree */
    private int[] numberOfBranches = null;
    /** number of end points voxels of every tree */
    private int[] numberOfEndPoints = null;
    /** number of junctions voxels of every tree*/
    //    private int[] numberOfJunctionVoxels = null;
    /** number of slab voxels of every specific tree */
    private int[] numberOfSlabs = null;	
    /** number of junctions of every specific tree*/
    //    private int[] numberOfJunctions = null;
    /** number of triple points in every tree */
    private int[] numberOfTriplePoints = null;
    /** list of end points in every tree */
    private ArrayList <int[]> endPointsTree [] = null;
    /** list of junction voxels in every tree */
    private ArrayList <int[]> junctionVoxelTree [] = null;
    /** list of special slab coordinates where circular tree starts */
    private ArrayList <int[]> startingSlabTree [] = null;

    /** list of all branch lengths */
    private ArrayList <double[]> listOfBranchLengths = null;

    /** sum of branch lengths for each tree */
    private double[] branchLength = null;

    /** average branch length */
    private double[] averageBranchLength = null;

    /** maximum branch length */
    private double[] maximumBranchLength = null;

    /** list of end point coordinates in the entire image */
    private ArrayList <int[]> listOfEndPoints = new ArrayList<int[]>();
    /** list of junction coordinates in the entire image */
    private ArrayList <int[]> listOfJunctionVoxels = new ArrayList<int[]>();
    /** list of slab coordinates in the entire image */
    private ArrayList <int[]> listOfSlabVoxels = new ArrayList<int[]>();
    /** list of slab coordinates in the entire image */
    private ArrayList <int[]> listOfStartingSlabVoxels = new ArrayList<int[]>();

    /** list of groups of junction voxels that belong to the same tree junction (in every tree) */
    private ArrayList < ArrayList <int[]> > listOfSingleJunctions[] = null;

    /** stack image containing the corresponding skeleton tags (end point, junction or slab) */
    private ImageStack taggedImage = null;

    /** auxiliary temporary point */
    private int[] auxPoint = null;
    /** largest branch coordinates initial point */
    private int[][] initialPoint = null;
    /** largest branch coordinates final point */
    private int[][] finalPoint = null;

    /** number of trees (skeletons) in the image */
    private int numOfTrees = 0;

    /* -----------------------------------------------------------------------*/
    /**
     * This method is called once when the filter is loaded.
     * 
     * @param arg argument specified for this plugin
     * @param imp currently active image
     * @return flag word that specifies the filters capabilities
     */
    public int setup(String arg, ImagePlus imp) 
    {
	this.imRef = imp;
	this.VOXEL_WIDTH = this.imRef.getCalibration().pixelWidth;
	this.VOXEL_HEIGHT = this.imRef.getCalibration().pixelHeight;
	this.VOXEL_DEPTH = this.imRef.getCalibration().pixelDepth;

	if (arg.equals("about")) 
	{
	    showAbout();
	    return DONE;
	}

	return DOES_8G; 
    } /* end setup */

    /* -----------------------------------------------------------------------*/
    /**
     * Process the image: tag skeleton and show results.
     * 
     * @see ij.plugin.filter.PlugInFilter#run(ij.process.ImageProcessor)
     */
    public void run(ImageProcessor ip) 
    {

	this.width = this.imRef.getWidth();
	this.height = this.imRef.getHeight();
	this.depth = this.imRef.getStackSize();
	this.inputImage = this.imRef.getStack();

	// initialize visit flags
	this.visited = new boolean[this.width][this.height][this.depth];

	// Prepare data: classify voxels and tag them.
	GenericDialog gd = new GenericDialog("Prune?");
	gd.addCheckbox("Prune Ends", true);
	gd.showDialog();
	boolean doPrune = gd.getNextBoolean();

	if (doPrune){
	    ImageStack stack1 = tagImage(this.inputImage);
	    this.taggedImage = pruneEndBranches(stack1);
	}
	else
	    this.taggedImage = tagImage(this.inputImage);

	// Mark trees
	ImageStack treeIS = markTrees(this.taggedImage);

	// Ask memory for every tree
	this.numberOfBranches = new int[this.numOfTrees];
	this.numberOfEndPoints = new int[this.numOfTrees];
	//	this.numberOfJunctionVoxels = new int[this.numOfTrees];
	//	this.numberOfJunctions = new int[this.numOfTrees];
	this.numberOfSlabs = new int[this.numOfTrees];
	this.numberOfTriplePoints = new int[this.numOfTrees];
	this.branchLength = new double[this.numOfTrees];
	this.averageBranchLength = new double[this.numOfTrees];
	this.maximumBranchLength = new double[this.numOfTrees];
	this.initialPoint = new int[this.numOfTrees][];
	this.finalPoint = new int[this.numOfTrees][];
	this.endPointsTree = new ArrayList[this.numOfTrees];		
	this.junctionVoxelTree = new ArrayList[this.numOfTrees];
	this.startingSlabTree = new ArrayList[this.numOfTrees];
	this.listOfSingleJunctions = new ArrayList[this.numOfTrees];
	for(int i = 0; i < this.numOfTrees; i++)
	{
	    this.endPointsTree[i] = new ArrayList <int[]>();
	    this.junctionVoxelTree[i] = new ArrayList <int[]>();
	    this.startingSlabTree[i] = new ArrayList <int[]>();
	    this.listOfSingleJunctions[i] = new ArrayList < ArrayList <int[]> > ();
	}

	// Divide groups of end-points and junction voxels
	if(this.numOfTrees > 1)
	    divideVoxelsByTrees(treeIS);
	else
	{
	    this.endPointsTree[0] = this.listOfEndPoints;
	    this.junctionVoxelTree[0] = this.listOfJunctionVoxels;
	}

	// Calculate number of junctions (group neighbor junction voxels)
	groupJunctions(treeIS);

	// Visit skeleton and measure distances.
	for(int i = 0; i < this.numOfTrees; i++)
	    visitSkeleton(taggedImage, treeIS, i+1);

	// Calculate triple points (junctions with exactly 3 branches)
	calculateTriplePoints();

	// Show results table
	showResults();

	// Show tags image.
	ImagePlus tagIP = new ImagePlus("Tagged skeleton", taggedImage);
	tagIP.setCalibration(this.imRef.getCalibration());

	// We apply the Fire LUT and reset the min and max to be between 0-255.
	tagIP.show();
	IJ.run("Fire");
	tagIP.resetDisplayRange();
	tagIP.updateAndDraw();


    } /* end run */

    /* -----------------------------------------------------------------------*/
    /**
     * Divide end point, junction and special (starting) slab voxels in the 
     * corresponding tree lists
     * 
     *  @param treeIS tree image
     */
    private void divideVoxelsByTrees(ImageStack treeIS) 
    {
	// Add end points to the corresponding tree
	for(int i = 0; i < this.listOfEndPoints.size(); i++)
	{
	    final int[] p = this.listOfEndPoints.get(i);
	    this.endPointsTree[getShortPixel(treeIS, p) - 1].add(p);
	}

	// Add junction voxels to the corresponding tree
	for(int i = 0; i < this.listOfJunctionVoxels.size(); i++)
	{
	    final int[] p = this.listOfJunctionVoxels.get(i);			
	    this.junctionVoxelTree[getShortPixel(treeIS, p) - 1].add(p);  //TODO fix AIOOB here
	}

	// Add special slab voxels to the corresponding tree
	for(int i = 0; i < this.listOfStartingSlabVoxels.size(); i++)
	{
	    final int[] p = this.listOfStartingSlabVoxels.get(i);			
	    this.startingSlabTree[getShortPixel(treeIS, p) - 1].add(p);
	}

    } // end divideVoxelsByTrees

    /* -----------------------------------------------------------------------*/
    /**
     * Show results table
     */
    private void showResults() 
    {
	String unit = this.imRef.getCalibration().getUnit();
	ResultsTable rt = new ResultsTable();

	String[] head = {"Skeleton", "# Branches","# Junctions", "# End-point voxels",
		"# Junction voxels","# Slab voxels","Average Branch Length ("+unit+")", 
		"# Triple points", "Maximum Branch Length ("+unit+")", "Sum Branch Length ("+unit+")"};

	for (int i = 0; i < head.length; i++)
	    rt.setHeading(i,head[i]);	

	for(int i = 0 ; i < this.numOfTrees; i++)
	{
	    rt.incrementCounter();

	    rt.addValue(1, this.numberOfBranches[i]);        
	    rt.addValue(2, this.listOfSingleJunctions[i].size());
	    rt.addValue(3, this.numberOfEndPoints[i]);
	    rt.addValue(4, this.junctionVoxelTree[i].size());
	    rt.addValue(5, this.numberOfSlabs[i]);
	    rt.addValue(6, this.averageBranchLength[i]);
	    rt.addValue(7, this.numberOfTriplePoints[i]);
	    rt.addValue(8, this.maximumBranchLength[i]);
	    rt.addValue(9, this.branchLength[i]);

	    rt.show("Results");

	    IJ.log("--- Skeleton #" + (i+1) + " ---");
	    IJ.log("Coordinates of the largest branch:");
	    IJ.log("Initial point: (" + (this.initialPoint[i][0] * this.imRef.getCalibration().pixelWidth) + ", " //TODO fix NPE here 
		    + (this.initialPoint[i][1] * this.imRef.getCalibration().pixelHeight) + ", "
		    + (this.initialPoint[i][2] * this.imRef.getCalibration().pixelDepth) + ")" );
	    IJ.log("Final point: (" + (this.finalPoint[i][0] * this.imRef.getCalibration().pixelWidth) + ", " 
		    + (this.finalPoint[i][1] * this.imRef.getCalibration().pixelHeight) + ", "
		    + (this.finalPoint[i][2] * this.imRef.getCalibration().pixelDepth) + ")" );
	    IJ.log("Euclidean distance: " + this.calculateDistance(this.initialPoint[i], this.finalPoint[i]));
	}
    }

    /* -----------------------------------------------------------------------*/
    /**
     * Visit skeleton from end points and register measures.
     * 
     * @param taggedImage
     * @param treeImage skeleton image with tree classification
     * @param currentTree number of the tree to be visited
     */
    private void visitSkeleton(ImageStack taggedImage, ImageStack treeImage, int currentTree) 
    {

	// IJ.log(" Analyzing tree number " + currentTree);

	// tree index
	final int iTree = currentTree - 1;
	// length of branches
	this.branchLength[iTree] = 0;

	this.maximumBranchLength[iTree] = 0;
	//	this.numberOfEndPoints[iTree] = this.endPointsTree[iTree].size();
	//	this.numberOfJunctionVoxels[iTree] = this.junctionVoxelTree[iTree].size();
	this.numberOfSlabs[iTree] = 0;
	this.listOfBranchLengths = new ArrayList<double[]>();

	// Visit branches starting at end points
	for(int i = 0; i < this.endPointsTree[iTree].size(); i++)
	{			
	    final int[] endPointCoord = this.endPointsTree[iTree].get(i);
	    IJ.log("Checking branch at endPoint "+pointToString(endPointCoord));
	    // Skip when visited
	    if(isVisited(endPointCoord))
	    {
		//if(this.initialPoint[iTree] == null)
		//	IJ.error("WEIRD:" + " (" + endPointCoord[0] + ", " + endPointCoord[1] + ", " + endPointCoord[2] + ")");
		//IJ.log("visited = (" + endPointCoord[0] + ", " + endPointCoord[1] + ", " + endPointCoord[2] + ")");
		IJ.log("...but endPoint has been visited, skipping");
		continue;
	    }

	    // Otherwise, visit branch until next junction or end point.
	    final double length = visitBranch(endPointCoord, iTree);
	    IJ.log("...visited branch of length "+length);
	    // If length is 0, it means the tree is formed by only one voxel.
	    if(length == 0)
	    {
		this.initialPoint[iTree] = this.finalPoint[iTree] = endPointCoord;
		//		IJ.log("Found tree:branch ("+iTree+":"+this.numberOfBranches[iTree]+") " +
		//			"using end point at ("+endPointCoord[0]+", "+endPointCoord[1]+", "+endPointCoord[2]+")" +
		//			" of length "+length);
		continue;
	    }

	    // increase number of branches
	    this.numberOfBranches[iTree]++;
	    IJ.log("Added a branch ("+this.numberOfBranches[iTree]+") to tree "+iTree);
	    this.branchLength[iTree] += length;
	    double[] lA = {length, 0}; //2nd element is the bin. 
	    this.listOfBranchLengths.add(lA);

	    //	    IJ.log("Found tree:branch ("+iTree+":"+this.numberOfBranches[iTree]+") " +
	    //		    "using end point at ("+endPointCoord[0]+", "+endPointCoord[1]+", "+endPointCoord[2]+")" +
	    //		    " of length "+length);

	    // update maximum branch length
	    if(length > this.maximumBranchLength[iTree])
	    {
		this.maximumBranchLength[iTree] = length;
		this.initialPoint[iTree] = endPointCoord;
		this.finalPoint[iTree] = this.auxPoint;
	    }
	}


	// Now visit branches starting at junctions
	//current model works from junction voxels,
	//less problematic to work from junctions
	
	//filter out non-branching junction voxels by marking them visited
	Skeletonize3D_ sk = new Skeletonize3D_();
	int[] LUT = new int[256];
	sk.fillEulerLUT(LUT);
	double[] sumJunc = {0,0,0}; // sum of junction coordinates
	double[] centroidJunc = {0,0,0}; // centroid of junction coordinates (real units)
	for(int i = 0; i < this.junctionVoxelTree[iTree].size(); i++){
	    final int[] junctionCoord = this.junctionVoxelTree[iTree].get(i);
	    int x = junctionCoord[0], y = junctionCoord[1], z = junctionCoord[2];
	    byte[] neighbors = getNeighborhood(this.taggedImage, x,y,z);
	    int nSlabs = 0, nEnds= 0;
	    for (int j = 0; j < 27; j++){
		switch(neighbors[j]){
		case SLAB: nSlabs++;
		case END_POINT: nEnds++;
		}		
	    }
	    //if there are no slabs or ends in neighborhood
	    //go to the next level of filtering
	    if (nSlabs == 0 && nEnds == 0){
		//Junction voxel is surrounded by 0's or junction voxels
		//so see if it is really a branch stub and not in the middle of a junction
		for (int j = 0; j < 27; j++){
		    if (neighbors[j] > 0) neighbors[j] = 1;
		}
		if (!sk.isEulerInvariant(neighbors, LUT)){
		    //junction voxels in the middle of a junction shouldn't contribute to 
		    //branches
		    setVisited(x,y,z, true);
		    IJ.log("Set junction voxel "+pointToString(junctionCoord)+" visited");
		}
	    }
	}
	for(int i = 0; i < this.junctionVoxelTree[iTree].size(); i++)
	{
	    final int[] junctionCoord = this.junctionVoxelTree[iTree].get(i);					
	    IJ.log("Checking branch at junctionCoord "+pointToString(junctionCoord));
	    // Mark junction as visited
	    setVisited(junctionCoord, true);
	    
	    boolean isJunction = false;
	    if (getPixel(this.taggedImage, junctionCoord) == Analyze_Skeleton.JUNCTION)
		isJunction = true;		
	    
	    int[] nextPoint = getNextUnvisitedVoxel(junctionCoord, isJunction);
	    
	    //TODO prevent traversing across a junction between adjacent junction voxels
	   
	    while(nextPoint != null)
	    {
		IJ.log("...heading towards "+pointToString(nextPoint));
		this.branchLength[iTree] += calculateDistance(junctionCoord, nextPoint);								

		double length = visitBranch(nextPoint, iTree);

		this.branchLength[iTree] += length;

		// Increase number of branches
		if(length != 0)
		{
		    this.numberOfBranches[iTree]++;//does this count branches within a junction?
		    IJ.log("Added a branch ("+this.numberOfBranches[iTree]+") to tree "+iTree);
		    double[] lA = {length, 0}; //2nd element is the bin. 
		    this.listOfBranchLengths.add(lA);
		    // update maximum branch length
		    if(length > this.maximumBranchLength[iTree])
		    {
			this.maximumBranchLength[iTree] = length;
			this.initialPoint[iTree] = junctionCoord;
			this.finalPoint[iTree] = this.auxPoint;
		    }
		}
		nextPoint = getNextUnvisitedVoxel(junctionCoord, isJunction);
	    }					
	}

	// Finally visit branches starting at slabs (special case for circular trees)
	if(this.startingSlabTree[iTree].size() == 1)
	{
	    final int[] startCoord = this.startingSlabTree[iTree].get(0);					

	    // visit branch until finding visited voxel.
	    final double length = visitBranch(startCoord, iTree);

	    if(length != 0)
	    {				
		// increase number of branches
		this.numberOfBranches[iTree]++;
		this.branchLength[iTree] += length;

		double[] lA = {length, 0}; //2nd element is the bin. 
		this.listOfBranchLengths.add(lA);

		//		IJ.log("Found tree:branch ("+iTree+":"+this.numberOfBranches[iTree]+") " +
		//			"using slab at ("+startCoord[0]+", "+startCoord[1]+", "+startCoord[2]+")" +
		//			" of length "+length);
		// update maximum branch length
		if(length > this.maximumBranchLength[iTree])
		{
		    this.maximumBranchLength[iTree] = length;
		    this.initialPoint[iTree] = startCoord;
		    this.finalPoint[iTree] = this.auxPoint;
		}
	    }
	}						

	if(this.numberOfBranches[iTree] == 0)
	    return;
	// Average length
	this.averageBranchLength[iTree] = this.branchLength[iTree] / this.numberOfBranches[iTree];

    } /* end visitSkeleton */

    /* -----------------------------------------------------------------------*/
    /**
     * Color the different trees in the skeleton.
     * 
     * @param taggedImage
     * 
     * @return image with every tree tagged with a different number 
     */
    private ImageStack markTrees(ImageStack taggedImage) 
    {
	// Create output image
	ImageStack outputImage = new ImageStack(this.width, this.height, taggedImage.getColorModel());	
	for (int z = 0; z < depth; z++)
	{
	    outputImage.addSlice(taggedImage.getSliceLabel(z+1), new ShortProcessor(this.width, this.height));	
	}

	this.numOfTrees = 0;

	short color = 0;

	// Visit trees starting at end points
	//	IJ.log("Number of endPoints: "+this.listOfEndPoints.size());
	for(int i = 0; i < this.listOfEndPoints.size(); i++)
	{			
	    int[] endPointCoord = this.listOfEndPoints.get(i);

	    if(isVisited(endPointCoord))
		continue;

	    color++;

	    if(color == Short.MAX_VALUE)
	    {
		IJ.error("More than " + (Short.MAX_VALUE-1) +
			" skeletons in the image. AnalyzeSkeleton can only process up to "+ (Short.MAX_VALUE-1));
		return null;
	    }

	    // else, visit the entire tree.
	    int length = visitTree(endPointCoord, outputImage, color);

	    // increase number of trees			
	    this.numOfTrees++;
	}

	// Visit trees starting at junction points 
	// (some circular trees do not have end points)
	IJ.log("Number of junctionVoxels: "+this.listOfJunctionVoxels.size());
	for(int i = 0; i < this.listOfJunctionVoxels.size(); i++)
	{			
	    int[] junctionCoord = this.listOfJunctionVoxels.get(i);
	    if(isVisited(junctionCoord))
		continue;

	    color++;
	    IJ.log("Found a new tree ("+color+") starting at "+pointToString(junctionCoord));

	    if(color == Short.MAX_VALUE)
	    {
		IJ.error("More than " + (Short.MAX_VALUE-1) + " skeletons in the image. AnalyzeSkeleton can only process up to 255");
		return null;
	    }

	    // else, visit branch until next junction or end point.
	    int length = visitTree(junctionCoord, outputImage, color);

	    if(length == 0)
	    {
		color--; // the color was not used
		continue;
	    }

	    // increase number of trees			
	    this.numOfTrees++;
	}

	// Check for unvisited slab voxels
	// (just in case there are circular trees without junctions)
	//	IJ.log("Number of slabVoxels: "+this.listOfSlabVoxels.size());
	for(int i = 0; i < this.listOfSlabVoxels.size(); i++)
	{
	    int[] p = (int[]) this.listOfSlabVoxels.get(i);
	    if(isVisited(p) == false)
	    {
		// Mark that voxel as the start point of the circular skeleton
		this.listOfStartingSlabVoxels.add(p);

		color++;

		if(color == Short.MAX_VALUE)
		{
		    IJ.error("More than " + (Short.MAX_VALUE-1) + " skeletons in the image. AnalyzeSkeleton can only process up to 255");
		    return null;
		}

		// else, visit branch until next junction, end point or visited point
		int length = visitTree(p, outputImage, color);

		if(length == 0)
		{
		    color--; // the color was not used
		    continue;
		}

		// increase number of trees			
		this.numOfTrees++;
	    }
	}



	//	IJ.log("Number of trees = " + this.numOfTrees);

	// Show tree image.

	ImagePlus treesIP = new ImagePlus("Trees skeleton", outputImage);
	treesIP.show();

	// Set same calibration as the input image
	treesIP.setCalibration(this.imRef.getCalibration());

	// We apply the Fire LUT and reset the min and max to be between 0-255.
	IJ.run("Fire");

	//IJ.resetMinAndMax();
	treesIP.resetDisplayRange();
	treesIP.updateAndDraw();


	// Reset visited variable
	this.visited = null;
	this.visited = new boolean[this.width][this.height][this.depth];


	//IJ.log("Number of trees: " + this.numOfTrees + ", # colors = " + color);

	return outputImage;

    } /* end markTrees */


    /* --------------------------------------------------------------*/
    /**
     * Visit tree marking the voxels with a reference tree color
     * 
     * @param startingPoint starting tree point
     * @param outputImage 3D image to visit
     * @param color reference tree color
     * @return number of voxels in the tree
     */
    private int visitTree(int[] startingPoint, ImageStack outputImage,
	    short color) 
    {
	int numOfVoxels = 0;

	IJ.log("visiting " + pointToString(startingPoint) + " color = " + color);

	//	if(isVisited(startingPoint))	
	//	    return 0;
	// Set pixel color
	this.setPixel(outputImage, startingPoint[0], startingPoint[1], startingPoint[2], color);
	setVisited(startingPoint, true);

	ArrayList <int[]> toRevisit = new ArrayList <int []>();

	int[] nextPoint = getNextUnvisitedVoxel(startingPoint);

	while(nextPoint != null || toRevisit.size() != 0)
	{
	    if(nextPoint != null)
	    {
		if(!isVisited(nextPoint))
		{
		    numOfVoxels++;

		    //IJ.log("visiting " + pointToString(nextPoint)+ " color = " + color);

		    // Set color and visit flat
		    this.setPixel(outputImage, nextPoint[0], nextPoint[1], nextPoint[2], color);
		    setVisited(nextPoint, true);

		    // If it is a junction, add it to the revisit list
		    if(isJunction(nextPoint))
			toRevisit.add(nextPoint);

		    // Calculate next point to visit
		    nextPoint = getNextUnvisitedVoxel(nextPoint);
		}				
	    }
	    else // revisit list
	    {				
		nextPoint = toRevisit.get(0);
		//IJ.log("visiting " + pointToString(nextPoint)+ " color = " + color);

		// Calculate next point to visit
		nextPoint = getNextUnvisitedVoxel(nextPoint);
		// Maintain junction in the list until there is no more branches
		if (nextPoint == null)
		    toRevisit.remove(0);									
	    }				
	}

	return numOfVoxels;
    }

    /* -----------------------------------------------------------------------*/
    /**
     * Visit a branch and calculate length
     * 
     * @param startingPoint starting coordinates
     * @return branch length
     * 
     * @deprecated
     */
    private double visitBranch(int[] startingPoint) 
    {
	double length = 0;

	// mark starting point as visited
	setVisited(startingPoint, true);

	// Get next unvisited voxel
	int[] nextPoint = getNextUnvisitedVoxel(startingPoint);

	if (nextPoint == null)
	    return 0;

	int[] previousPoint = startingPoint;

	// We visit the branch until we find an end point or a junction
	while(nextPoint != null && isSlab(nextPoint))
	{
	    // Add length
	    length += calculateDistance(previousPoint, nextPoint);

	    // Mark as visited
	    setVisited(nextPoint, true);

	    // Move in the graph
	    previousPoint = nextPoint;			
	    nextPoint = getNextUnvisitedVoxel(previousPoint);			
	}


	if(nextPoint != null)
	{
	    // Add distance to last point
	    length += calculateDistance(previousPoint, nextPoint);

	    // Mark last point as visited
	    setVisited(nextPoint, true);
	}

	this.auxPoint = previousPoint;

	return length;
    } /* end visitBranch*/

    /* -----------------------------------------------------------------------*/
    /**
     * Visit a branch and calculate length in a specific tree
     * 
     * @param startingPoint starting coordinates
     * @param iTree tree index
     * @return branch length
     */
    private double visitBranch(int[] startingPoint, int iTree) 
    {
	//IJ.log("startingPoint = (" + startingPoint[0] + ", " + startingPoint[1] + ", " + startingPoint[2] + ")");
	double length = 0;

	// mark starting point as visited
	setVisited(startingPoint, true);

	// Get next unvisited voxel
	int[] nextPoint = getNextUnvisitedVoxel(startingPoint);

	if (nextPoint == null)
	    return 0;

	int[] previousPoint = startingPoint;

	// We visit the branch until we find an end point or a junction
	while(nextPoint != null && isSlab(nextPoint))
	{
	    this.numberOfSlabs[iTree]++;

	    // Add length
	    length += calculateDistance(previousPoint, nextPoint);

	    // Mark as visited
	    setVisited(nextPoint, true);

	    // Move in the graph
	    previousPoint = nextPoint;			
	    nextPoint = getNextUnvisitedVoxel(previousPoint);			
	}

	if(nextPoint != null)
	{
	    // Add distance to last point
	    length += calculateDistance(previousPoint, nextPoint);

	    // Mark last point as visited
	    setVisited(nextPoint, true);
	    
	    IJ.log("...and ending at "+pointToString(nextPoint));
	}

	this.auxPoint = previousPoint;

	return length;
    } /* end visitBranch*/	

    /* -----------------------------------------------------------------------*/
    /**
     * Calculate distance between two points in 3D
     * 
     * @param point1 first point coordinates
     * @param point2 second point coordinates
     * @return distance (in the corresponding units)
     */
    private double calculateDistance(int[] point1, int[] point2) 
    {	
	double dx = (point1[0] - point2[0]) * this.VOXEL_WIDTH;
	double dy = (point1[1] - point2[1]) * this.VOXEL_HEIGHT;
	double dz = (point1[2] - point2[2]) * this.VOXEL_DEPTH;
	return Math.sqrt(  dx * dx + dy * dy + dz * dz );
    }

    /* -----------------------------------------------------------------------*/
    /**
     * Calculate number of junction skipping neighbor junction voxels
     * 
     * @param treeIS tree stack
     */
    private void groupJunctions(ImageStack treeIS) 
    {
	for (int iTree = 0; iTree < this.numOfTrees; iTree++)
	{
	    // Visit all junction voxels in tree iTree
	    for(int i = 0; i < this.junctionVoxelTree[iTree].size(); i++)
	    {
		//pi is the coordinates of junction voxel i in tree iTree
		int[] pi = this.junctionVoxelTree[iTree].get(i);
		boolean grouped = false;
		int lowGroup = 0;
		//visit all the junctions j in tree iTree
		for(int j = 0; j < this.listOfSingleJunctions[iTree].size(); j++)
		{
		    //groupOfJunctions is the list of voxels in a junction
		    ArrayList <int[]> groupOfJunctions = this.listOfSingleJunctions[iTree].get(j);

		    //visit all the voxels k in a junction (so far)
		    for(int k = 0; k < groupOfJunctions.size(); k++)
		    {
			int[] pk = groupOfJunctions.get(k);

			// if the test voxel pi is a neighbor of 
			//a voxel pk in junction j, add pi to j   
			if(isNeighbor(pi, pk))
			{
			    groupOfJunctions.add(pi);
			    grouped = true;
			    //can stop checking in this junction now
			    break;
			}
		    }

		    if(grouped){
			//we added pi to junction j
			lowGroup = j;
			//stop checking for first junction group assignment
			break;
		    }
		}
		//check whether pi has neighbors in any other junctions
		//and join j and the other junctions together
		if (grouped){
		    ArrayList <int[]> lowGroupList = this.listOfSingleJunctions[iTree].get(lowGroup);
		    for (int l = lowGroup + 1; l < this.listOfSingleJunctions[iTree].size(); l++){
			ArrayList <int[]> groupOfJunctions = this.listOfSingleJunctions[iTree].get(l);
			//visit all the voxels k in a junction (so far)
			for(int k = 0; k < groupOfJunctions.size(); k++)
			{
			    int[] lpk = groupOfJunctions.get(k);

			    // if the test voxel pi is a neighbor of 
			    //a voxel lpk in junction l, add l to j
			    //and remove l from the list of junctions
			    if(isNeighbor(pi, lpk))
			    {
				grouped = true;
				lowGroupList.addAll(groupOfJunctions);
				this.listOfSingleJunctions[iTree].remove(l);
				//can stop checking junction l now
				break;
			    }
			}
		    }
		}

		if(!grouped)
		{
		    ArrayList <int[]> newGroup = new ArrayList<int[]>();
		    newGroup.add(pi);
		    this.listOfSingleJunctions[iTree].add(newGroup);
		}
	    }
	}
    }	

    /* -----------------------------------------------------------------------*/
    /**
     * Calculate number of triple points in the skeleton. Triple points are
     * junctions with exactly 3 branches.
     */
    private void calculateTriplePoints() 
    {
	for (int iTree = 0; iTree < this.numOfTrees; iTree++)
	{			
	    // Visit the groups of junction voxels
	    for(int i = 0; i < this.listOfSingleJunctions[iTree].size(); i ++)
	    {

		ArrayList <int[]> groupOfJunctions = this.listOfSingleJunctions[iTree].get(i);

		// Count the number of slab neighbors of every voxel in the group
		int nSlab = 0;
		for(int j = 0; j < groupOfJunctions.size(); j++)
		{
		    int[] pj = groupOfJunctions.get(j);

		    // Get neighbors and check the slabs
		    byte[] neighborhood = this.getNeighborhood(this.taggedImage, pj[0], pj[1], pj[2]);
		    for(int k = 0; k < 27; k++)
			if (neighborhood[k] == Analyze_Skeleton.SLAB)
			    nSlab++;
		}
		// If the junction has only 3 slab neighbors, then it is a triple point
		if (nSlab == 3)	
		    this.numberOfTriplePoints[iTree] ++;

	    }		

	}

    }// end calculateTriplePoints


    /* -----------------------------------------------------------------------*/
    /**
     * Calculate if two points are neighbors
     * 
     * @param point1 first point
     * @param point2 second point
     * @return true if the points are neighbors (26-pixel neighborhood)
     */
    private boolean isNeighbor(int[] point1, int[] point2) 
    {		
	int dx = point1[0] - point2[0];
	int dy = point1[1] - point2[1];
	int dz = point1[2] - point2[2];
	return Math.sqrt(dx * dx + dy * dy + dz * dz) <= Analyze_Skeleton.SQRT_3;
    }

    /* -----------------------------------------------------------------------*/
    /**
     * Check if the point is slab
     *  
     * @param point actual point
     * @return true if the point has slab status
     */
    private boolean isSlab(int[] point) 
    {		
	return getPixel(this.taggedImage, point[0], point[1], point[2]) == Analyze_Skeleton.SLAB;
    }

    /* -----------------------------------------------------------------------*/
    /**
     * Check if the point is a junction
     *  
     * @param point actual point
     * @return true if the point has slab status
     */
    private boolean isJunction(int[] point) 
    {		
	return getPixel(this.taggedImage, point[0], point[1], point[2]) == Analyze_Skeleton.JUNCTION;
    }	

    /* -----------------------------------------------------------------------*/
    /**
     * Get next unvisited neighbor voxel 
     * 
     * @param point starting point
     * @return unvisited neighbor or null if all neighbors are visited
     */
    private int[] getNextUnvisitedVoxel(int[] point) 
    {
	int[] unvisitedNeighbor = null;
	
	int xp = point[0], yp = point[1], zp = point[2];
	// Check neighbors status
	for(int x = -1; x < 2; x++)
	    for(int y = -1; y < 2; y++)
		for(int z = -1; z < 2; z++)
		{
		    if(x == 0 && y == 0 && z == 0)
			continue;
		       
		    if(getPixel(this.taggedImage, xp + x, yp + y, zp + z) != 0
			    && isVisited(xp + x, yp + y, zp + z) == false)						
		    {					
			unvisitedNeighbor = new int[]{xp + x, yp + y, zp + z};
			break;
		    }

		}

	return unvisitedNeighbor;
    }/* end getNextUnvisitedVoxel */

    /**
     * Get next unvisited neighbor voxel 
     * Modified to prevent traversing from one junction voxel to another
     *  
     * @param point starting point
     * @param isJunction true if the starting point is a junction voxel
     * @return unvisited neighbor or null if all neighbors are visited
     */
    private int[] getNextUnvisitedVoxel(int[] point, boolean isJunction) 
    {
	int[] unvisitedNeighbor = null;
	int xp = point[0], yp = point[1], zp = point[2];
	// Check neighbors status
	for(int x = -1; x < 2; x++)
	    for(int y = -1; y < 2; y++)
		for(int z = -1; z < 2; z++)
		{
		    if(x == 0 && y == 0 && z == 0)
			continue;
		    byte pixel = getPixel(this.taggedImage, xp + x, yp + y, zp + z);
		    
		    if (pixel == Analyze_Skeleton.JUNCTION && isJunction)
			continue;
		    
		    if(pixel != 0 && isVisited(xp + x, yp + y, zp + z) == false)						
		    {					
			unvisitedNeighbor = new int[]{xp + x, yp + y, zp + z};
			break;
		    }

		}

	return unvisitedNeighbor;
    }/* end getNextUnvisitedVoxel */

    
    /* -----------------------------------------------------------------------*/
    /**
     * Check if a voxel is visited taking into account the borders. 
     * Out of range voxels are considered as visited. 
     * 
     * @param point
     * @return true if the voxel is visited
     */
    private boolean isVisited(int [] point) 
    {
	return isVisited(point[0], point[1], point[2]);
    }

    /* -----------------------------------------------------------------------*/
    /**
     * Check if a voxel is visited taking into account the borders. 
     * Out of range voxels are considered as visited. 
     * 
     * @param x x- voxel coordinate
     * @param y y- voxel coordinate
     * @param z z- voxel coordinate
     * @return true if the voxel is visited
     */
    private boolean isVisited(int x, int y, int z) 
    {
	if(x >= 0 && x < this.width && y >= 0 && y < this.height && z >= 0 && z < this.depth)
	    return this.visited[x][y][z];
	return true;
    }


    /* -----------------------------------------------------------------------*/
    /**
     * Set value in the visited flags matrix
     * 
     * @param x x- voxel coordinate
     * @param y y- voxel coordinate
     * @param z z- voxel coordinate
     * @param b
     */
    private void setVisited(int x, int y, int z, boolean b) 
    {
	if(x >= 0 && x < this.width && y >= 0 && y < this.height && z >= 0 && z < this.depth)
	    this.visited[x][y][z] = b;		
    }

    /* -----------------------------------------------------------------------*/
    /**
     * Set value in the visited flags matrix
     * 
     * @param point voxel coordinates
     * @param b visited flag value
     */
    private void setVisited(int[] point, boolean b) 
    {
	int x = point[0];
	int y = point[1];
	int z = point[2];

	setVisited(x, y, z, b);	
    }

    /* -----------------------------------------------------------------------*/
    /**
     * Tag skeleton dividing the voxels between end points, junctions and slab,
     *  
     * @param inputImage2 skeleton image to be tagged
     * @return tagged skeleton image
     */
    private ImageStack tagImage(ImageStack inputImage2) 
    {
	// Create output image
	ImageStack outputImage = new ImageStack(this.width, this.height, inputImage2.getColorModel());

	// Tag voxels
	for (int z = 0; z < depth; z++)
	{
	    outputImage.addSlice(inputImage2.getSliceLabel(z+1), new ByteProcessor(this.width, this.height));			
	    for (int x = 0; x < width; x++) 
		for (int y = 0; y < height; y++)
		{
		    if(getPixel(inputImage2, x, y, z) != 0)
		    {
			int numOfNeighbors = getNumberOfNeighbors(inputImage2, x, y, z);
			if(numOfNeighbors < 2)
			{
			    setPixel(outputImage, x, y, z, Analyze_Skeleton.END_POINT);
			    //			    this.totalNumberOfEndPoints++;
			    int[] endPoint = new int[]{x, y, z};
			    this.listOfEndPoints.add(endPoint);							
			}
			else if(numOfNeighbors > 2)
			{
			    setPixel(outputImage, x, y, z, Analyze_Skeleton.JUNCTION);
			    int[] junction = new int[]{x, y, z};
			    this.listOfJunctionVoxels.add(junction);	
			    //			    this.totalNumberOfJunctionVoxels++;
			}
			else
			{
			    setPixel(outputImage, x, y, z, Analyze_Skeleton.SLAB);
			    int[] slab = new int[]{x, y, z};
			    this.listOfSlabVoxels.add(slab);
			    //			    this.totalNumberOfSlabs++;
			}
		    }					
		}
	}

	return outputImage;
    }/* end tagImage */

    /*--------------------------------------------------------------------*/
    /**
     * Prune end branches
     * 
     * @param stack ImageStack skeleton image
     * 
     */
    private ImageStack pruneEndBranches(ImageStack stack){
	// Prepare Euler LUT [Lee94]
	Skeletonize3D_ sk = new Skeletonize3D_();
	int eulerLUT[] = new int[256]; 
	sk.fillEulerLUT(eulerLUT);
	int endPoints = this.listOfEndPoints.size();
	while (!this.listOfEndPoints.isEmpty()){
	    IJ.showStatus("Pruning end branches...");
	    IJ.showProgress(endPoints - this.listOfEndPoints.size(), endPoints);
	    prune:
		for (int i = 0; i < this.listOfEndPoints.size(); i++){
		    int[] endPoint = this.listOfEndPoints.get(i);
		    int x = endPoint[0], y = endPoint[1], z = endPoint[2];
		    //if the endpoint is now in a junctionVoxel's position
		    //remove it from the endpoint list  
		    for (int k = 0; k < this.listOfJunctionVoxels.size(); k++){
			int[] junctionVoxel = this.listOfJunctionVoxels.get(k);
			if (junctionVoxel[0] == x && junctionVoxel[1] == y && junctionVoxel[2] == z){
			    this.listOfEndPoints.remove(i);
			    // Check if point is Euler invariant, simple and not an endpoint
			    byte[] neighbors = this.getNeighborhood(stack, x, y, z);
			    byte nNeighbors = 0;
			    for (int l = 0; l < 27; l++){
				if (neighbors[l] > 0){
				    neighbors[l] = 1;
				    nNeighbors++;
				}
			    }
			    if(sk.isEulerInvariant(neighbors, eulerLUT) && 
				    sk.isSimplePoint(neighbors) &&
				    nNeighbors > 2){
				//delete the junction point
				this.listOfJunctionVoxels.remove(k);
				setPixel(stack, x, y, z, (byte) 0);
			    }
			    continue prune;
			}
		    }
		    if (getNumberOfNeighbors(stack, x, y, z) == 1){
			//remove the end voxel
			setPixel(stack, x, y, z, (byte) 0);
			//remove end voxel from list of slabs
			for (int j = 0; j < this.listOfSlabVoxels.size(); j++){
			    int[] slabVoxel = this.listOfSlabVoxels.get(j);
			    if (slabVoxel[0] == x && slabVoxel[1] == y && slabVoxel[2] == z){
				this.listOfSlabVoxels.remove(j);
				break;
			    }
			}
			//get the values of the neighbors 
			byte[] nHood = getNeighborhood(stack, x, y, z);
			//get the coordinates of the single neighbor
			for (int p = 0; p < 27; p++){
			    if (nHood[p] != 0){
				//translate the neighbourhood index 
				//into new endpoint coordinates
				switch(p){
				case  0:   x -= 1; y -= 1; z -= 1; break;
				case  1:   y -= 1; z -= 1; break;
				case  2:   x += 1; y -= 1; z -= 1; break;
				case  3:   x -= 1; z -= 1; break;
				case  4:   z -= 1; break;
				case  5:   x += 1; z -= 1; break;
				case  6:   x -= 1; y += 1; z -= 1; break;
				case  7:   y += 1; z -= 1; break;
				case  8:   x += 1; y += 1; z -= 1; break;
				case  9:   x -= 1; y -= 1; break;
				case 10:   y -= 1; break;
				case 11:   x += 1; y -= 1; break;
				case 12:   x -= 1; break;
				case 13:   break;
				case 14:   x += 1; break;
				case 15:   x -= 1; y += 1; break;
				case 16:   y += 1; break;
				case 17:   x += 1; y += 1; break;
				case 18:   x -= 1; y -= 1; z += 1; break;
				case 19:   y -= 1; z += 1; break;
				case 20:   x += 1; y -= 1; z += 1; break;
				case 21:   x -= 1; z += 1; break;
				case 22:   z += 1; break;
				case 23:   x += 1; z += 1; break;
				case 24:   x -= 1; y += 1; z += 1; break;
				case 25:   y += 1; z += 1; break;
				case 26:   x += 1; y += 1; z += 1; break;			    
				}
				endPoint[0] = x;
				endPoint[1] = y;
				endPoint[2] = z;
				this.listOfEndPoints.set(i, endPoint);
				break;
			    }
			}
		    } else if (getNumberOfNeighbors(stack, x, y, z) > 1){
			this.listOfEndPoints.remove(i);
		    } else {
			//set a remaining endPoint to a slab
			this.listOfEndPoints.remove(i);
			this.listOfSlabVoxels.add(endPoint);
		    }
		}
	}
	//clean up isolated slab points
	for (int i = 0; i < this.listOfSlabVoxels.size(); i++){
	    int[] slab = this.listOfSlabVoxels.get(i);
	    int x = slab[0], y = slab[1], z = slab[2];
	    if (getNumberOfNeighbors(stack, x, y, z) == 0){
		this.listOfSlabVoxels.remove(i);
		this.listOfEndPoints.add(slab);
		setPixel(stack, x, y, z, END_POINT);
	    }
	}
	IJ.log("Number of endPoints after pruning: "+this.listOfEndPoints.size());
	return stack;
    }	


    /* -----------------------------------------------------------------------*/
    /**
     * Get number of neighbors of a voxel in a 3D image (0 border conditions) 
     * 
     * @param image 3D image (ImageStack)
     * @param x x- coordinate
     * @param y y- coordinate
     * @param z z- coordinate (in image stacks the indexes start at 1)
     * @return corresponding 27-pixels neighborhood (0 if out of image)
     */
    private int getNumberOfNeighbors(ImageStack image, int x, int y, int z)
    {
	int n = 0;
	byte[] neighborhood = getNeighborhood(image, x, y, z);

	for(int i = 0; i < 27; i ++)
	    if(neighborhood[i] != 0)
		n++;
	// We return n-1 because neighborhood includes the actual voxel.
	return (n-1);			
    }


    /* -----------------------------------------------------------------------*/
    /**
     * Get neighborhood of a pixel in a 3D image (0 border conditions) 
     * 
     * @param image 3D image (ImageStack)
     * @param x x- coordinate
     * @param y y- coordinate
     * @param z z- coordinate (in image stacks the indexes start at 1)
     * @return corresponding 27-pixels neighborhood (0 if out of image)
     */
    private byte[] getNeighborhood(ImageStack image, int x, int y, int z)
    {
	byte[] neighborhood = new byte[27];

	neighborhood[ 0] = getPixel(image, x-1, y-1, z-1);
	neighborhood[ 1] = getPixel(image, x  , y-1, z-1);
	neighborhood[ 2] = getPixel(image, x+1, y-1, z-1);

	neighborhood[ 3] = getPixel(image, x-1, y,   z-1);
	neighborhood[ 4] = getPixel(image, x,   y,   z-1);
	neighborhood[ 5] = getPixel(image, x+1, y,   z-1);

	neighborhood[ 6] = getPixel(image, x-1, y+1, z-1);
	neighborhood[ 7] = getPixel(image, x,   y+1, z-1);
	neighborhood[ 8] = getPixel(image, x+1, y+1, z-1);

	neighborhood[ 9] = getPixel(image, x-1, y-1, z  );
	neighborhood[10] = getPixel(image, x,   y-1, z  );
	neighborhood[11] = getPixel(image, x+1, y-1, z  );

	neighborhood[12] = getPixel(image, x-1, y,   z  );
	neighborhood[13] = getPixel(image, x,   y,   z  );
	neighborhood[14] = getPixel(image, x+1, y,   z  );

	neighborhood[15] = getPixel(image, x-1, y+1, z  );
	neighborhood[16] = getPixel(image, x,   y+1, z  );
	neighborhood[17] = getPixel(image, x+1, y+1, z  );

	neighborhood[18] = getPixel(image, x-1, y-1, z+1);
	neighborhood[19] = getPixel(image, x,   y-1, z+1);
	neighborhood[20] = getPixel(image, x+1, y-1, z+1);

	neighborhood[21] = getPixel(image, x-1, y,   z+1);
	neighborhood[22] = getPixel(image, x,   y,   z+1);
	neighborhood[23] = getPixel(image, x+1, y,   z+1);

	neighborhood[24] = getPixel(image, x-1, y+1, z+1);
	neighborhood[25] = getPixel(image, x,   y+1, z+1);
	neighborhood[26] = getPixel(image, x+1, y+1, z+1);

	return neighborhood;
    } /* end getNeighborhood */

    /* -----------------------------------------------------------------------*/
    /**
     * Get pixel in 3D image (0 border conditions) 
     * 
     * @param image 3D image
     * @param x x- coordinate
     * @param y y- coordinate
     * @param z z- coordinate (in image stacks the indexes start at 1)
     * @return corresponding pixel (0 if out of image)
     */
    private byte getPixel(ImageStack image, int x, int y, int z)
    {
	if(x >= 0 && x < this.width && y >= 0 && y < this.height && z >= 0 && z < this.depth)
	    return ((byte[]) image.getPixels(z + 1))[x + y * this.width];
	else return 0;
    } /* end getPixel */


    /* -----------------------------------------------------------------------*/
    /**
     * Get pixel in 3D image (0 border conditions) 
     * 
     * @param image 3D image
     * @param x x- coordinate
     * @param y y- coordinate
     * @param z z- coordinate (in image stacks the indexes start at 1)
     * @return corresponding pixel (0 if out of image)
     */
    private short getShortPixel(ImageStack image, int x, int y, int z)
    {
	if(x >= 0 && x < this.width && y >= 0 && y < this.height && z >= 0 && z < this.depth)
	    return ((short[]) image.getPixels(z + 1))[x + y * this.width];
	else return 0;
    } /* end getShortPixel */

    /* -----------------------------------------------------------------------*/
    /**
     * Get pixel in 3D image (0 border conditions) 
     * 
     * @param image 3D image
     * @param point point to be evaluated
     * @return corresponding pixel (0 if out of image)
     */
    private short getShortPixel(ImageStack image, int [] point)
    {
	return getShortPixel(image, point[0], point[1], point[2]);
    } /* end getPixel */	

    /* -----------------------------------------------------------------------*/
    /**
     * Get pixel in 3D image (0 border conditions) 
     * 
     * @param image 3D image
     * @param point point to be evaluated
     * @return corresponding pixel (0 if out of image)
     */
    private byte getPixel(ImageStack image, int [] point)
    {
	return getPixel(image, point[0], point[1], point[2]);
    } /* end getPixel */


    /* -----------------------------------------------------------------------*/
    /**
     * Set pixel in 3D image 
     * 
     * @param image 3D image
     * @param x x- coordinate
     * @param y y- coordinate
     * @param z z- coordinate (in image stacks the indexes start at 1)
     * @param value pixel value
     */
    private void setPixel(ImageStack image, int x, int y, int z, byte value)
    {
	if(x >= 0 && x < this.width && y >= 0 && y < this.height && z >= 0 && z < this.depth)
	    ((byte[]) image.getPixels(z + 1))[x + y * this.width] = value;
    } /* end getPixel */

    /* -----------------------------------------------------------------------*/
    /**
     * Set pixel in 3D (short) image 
     * 
     * @param image 3D image
     * @param x x- coordinate
     * @param y y- coordinate
     * @param z z- coordinate (in image stacks the indexes start at 1)
     * @param value pixel value
     */
    private void setPixel(ImageStack image, int x, int y, int z, short value)
    {
	if(x >= 0 && x < this.width && y >= 0 && y < this.height && z >= 0 && z < this.depth)
	    ((short[]) image.getPixels(z + 1))[x + y * this.width] = value;
    } /* end getPixel */	


    /**
     * 
     */
    String pointToString(int[] p)
    {
	return new String("(" + p[0] + ", " + p[1] + ", " + p[2] + ")");
    }

    /* -----------------------------------------------------------------------*/
    /**
     * Show plug-in information.
     * 
     */
    void showAbout() 
    {
	IJ.showMessage(
		"About AnalyzeSkeleton...",
	"This plug-in filter analyzes a 2D/3D image skeleton.\n");
    } /* end showAbout */

}