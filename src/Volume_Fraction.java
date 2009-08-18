import java.awt.Rectangle;

import org.doube.bonej.ResultInserter;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;
import ij.plugin.filter.PlugInFilter;
import ij.gui.*;

public class Volume_Fraction implements PlugInFilter {
    ImagePlus imp;
    
    short minT, maxT;

    protected ImageStack stack;

    boolean showMask = true;

    public int setup(String arg, ImagePlus imp) {
	if (imp == null) {
	    IJ.noImage();
	    return DONE;
	}
	stack = imp.getStack();
	this.imp = imp;
	setThreshold(imp);
	return DOES_8G + DOES_16 + DOES_32;
    }

    public void run(ImageProcessor ip) {
	int startSlice = 1;
	int endSlice = stack.getSize();
	GenericDialog gd = new GenericDialog("Limit Slices");
	gd.addNumericField("Start Slice:", startSlice, 0);
	gd.addNumericField("End Slice:", endSlice, 0);
	gd.showDialog();
	if (gd.wasCanceled()) {
	    IJ.error("PlugIn canceled!");
	    return;
	}
	startSlice = (int) gd.getNextNumber();
	endSlice = (int) gd.getNextNumber();
	Rectangle r = ip.getRoi();
	ImageProcessor mask = ip.getMask();
	boolean hasMask = (mask != null);
	if (hasMask && showMask) {
	    (new ImagePlus("The Mask", mask)).show();
	}

	int rLeft = r.x;
	int rTop = r.y;
	int rRight = rLeft + r.width;
	int rBottom = rTop + r.height;

	long volTotal = 0;
	long volBone = 0;
	for (int s = startSlice; s <= endSlice; s++) {
	    ImageProcessor ipSlice = stack.getProcessor(s);
	    for (int v = rTop; v < rBottom; v++) {
		for (int u = rLeft; u < rRight; u++) {
		    if (!hasMask || mask.getPixel(u - rLeft, v - rTop) > 0) {
			volTotal++;
			if (ipSlice.getPixel(u, v) >= minT
				&& ipSlice.getPixel(u, v) <= maxT) {
			    volBone++;
			}
		    }
		}
	    }
	}
	double p = (double)volBone / (double)volTotal;	
	ResultInserter ri = new ResultInserter();
	ri.setResultInRow(imp, "BV/TV", p);
	ri.updateTable();
    }
    private void setThreshold(ImagePlus imp){
	if (imp != null	&& (imp.getType() == ImagePlus.GRAY8 || imp.getType() == ImagePlus.COLOR_256)) {
	    ImageStatistics stats = imp.getStatistics();
	    if (stats.histogram[0] + stats.histogram[255] == stats.pixelCount) {
		minT = 128;
		maxT = 255;
	    } else {
		IJ.run("Threshold...");
		new WaitForUserDialog("Set the threshold, then click OK.").show();
		minT = (short) imp.getProcessor().getMinThreshold();
		maxT = (short) imp.getProcessor().getMaxThreshold();		
	    }
	} else {
	    IJ.run("Threshold...");
		new WaitForUserDialog("Set the threshold, then click OK.").show();
		minT = (short) imp.getProcessor().getMinThreshold();
		maxT = (short) imp.getProcessor().getMaxThreshold();
	}
    }
}
