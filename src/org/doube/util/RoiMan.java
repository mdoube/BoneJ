package org.doube.util;

import java.awt.*;
import java.util.ArrayList;

import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Roi;
import ij.measure.Calibration;
import ij.plugin.frame.RoiManager;
import ij.process.ImageProcessor;

/**
 * Do useful things with ImageJ's ROI Manager
 *
 * @author Michael Doube
 */
public class RoiMan {
    public static final int NO_SLICE_NUMBER = -1;

    /**
     * Get the calibrated 3D coordinates of point ROIs from the ROI manager
     *
     * @param imp
     * @param roiMan
     * @return double[n][3] containing n (x, y, z) coordinates or null if there
     * are no points
     */
    public static double[][] getRoiManPoints(ImagePlus imp, RoiManager roiMan) {
        Calibration cal = imp.getCalibration();
        double vW = cal.pixelWidth;
        double vH = cal.pixelHeight;
        double vD = cal.pixelDepth;
        int nPoints = 0;
        Roi[] roiList = roiMan.getRoisAsArray();
        for (int i = 0; i < roiMan.getCount(); i++) {
            Roi roi = roiList[i];
            if (roi.getType() == 10) {
                nPoints++;
            }
        }
        if (nPoints == 0)
            return null;
        double[][] dataPoints = new double[nPoints][3];
        int j = 0;
        for (int i = 0; i < roiMan.getCount(); i++) {
            Roi roi = roiList[i];
            if (roi.getType() == 10) {
                Rectangle xy = roi.getBounds();
                dataPoints[j][0] = xy.getX() * vW;
                dataPoints[j][1] = xy.getY() * vH;
                dataPoints[j][2] = roi.getPosition() * vD;
                j++;
            }
        }
        return dataPoints;
    }

    /**
     * Return a list of ROIs that are active in the given slice, sliceNumber. ROIs without
     * a slice number are assumed to be active in all slices.
     *
     * @param roiMan
     * @param stack
     * @param sliceNumber
     * @return A list of active ROIs on the slice. Returns an empty list if
     * roiMan == null or stack == null Returns an empty list if slice
     * number is out of bounds
     */
    public static ArrayList<Roi> getSliceRoi(RoiManager roiMan, ImageStack stack, int sliceNumber) {
        ArrayList<Roi> roiList = new ArrayList<Roi>();

        if (roiMan == null || stack == null) {
            return roiList;
        }

        if (sliceNumber < 1 || sliceNumber > stack.getSize()) {
            return roiList;
        }

        Roi[] rois = roiMan.getRoisAsArray();
        for (Roi roi : rois) {
            String roiName = roi.getName();
            if (roiName == null) {
                continue;
            }
            int roiSliceNumber = roiMan.getSliceNumber(roiName);
            int roiPosition = roi.getPosition();
            if (roiSliceNumber == sliceNumber || roiSliceNumber == NO_SLICE_NUMBER || roiPosition == sliceNumber) {
                roiList.add(roi);
            }
        }
        return roiList;
    }

    /**
     * Find the x, y and z limits of the ROIs in the ROI Manager
     *
     * @param roiMan A collection of ROIs
     * @param stack  The stack inside which the ROIs must fit (max limits).
     * @return int[] containing x min, x max, y min, y max, z min and z max.
     * Returns null if roiMan == null or if roiMan.getCount() == 0
     * Returns null if stack == null If any of the ROIs contains no
     * slice information, z min is set to 1 and z max is set to
     * stack.getSize()
     */
    public static int[] getLimits(RoiManager roiMan, ImageStack stack) {
        if (roiMan == null || roiMan.getCount() == 0) {
            return null;
        }

        if (stack == null) {
            return null;
        }

        final int LAST_SLIDE = stack.getSize();
        final int DEFAULT_Z_MIN = 1;
        final int DEFAULT_Z_MAX = LAST_SLIDE;

        int xMin = Integer.MAX_VALUE;
        int xMax = 0;
        int yMin = Integer.MAX_VALUE;
        int yMax = 0;
        int zMin = DEFAULT_Z_MAX;
        int zMax = DEFAULT_Z_MIN;
        boolean noZRoi = false;
        boolean noValidRois = true;
        Roi[] rois = roiMan.getRoisAsArray();

        for (Roi roi : rois) {
            Rectangle r = roi.getBounds();
            boolean valid = getSafeRoiBounds(r, stack.getWidth(), stack.getHeight());

            if (!valid) {
                continue;
            }

            xMin = Math.min(r.x, xMin);
            xMax = Math.max(r.x + r.width, xMax);
            yMin = Math.min(r.y, yMin);
            yMax = Math.max(r.y + r.height, yMax);

            int slice = roiMan.getSliceNumber(roi.getName());
            if (slice >= 1 && slice <= LAST_SLIDE) {
                zMin = Math.min(slice, zMin);
                zMax = Math.max(slice, zMax);
                noValidRois = false;
            } else if (isActiveOnAllSlices(slice)) {
                noZRoi = true; // found a ROI with no Z info
                noValidRois = false;
            }
        }

        if (noValidRois) {
            return null;
        }

        if (noZRoi) {
            int[] limits = {xMin, xMax, yMin, yMax, DEFAULT_Z_MIN, DEFAULT_Z_MAX};
            return limits;
        } else {
            int[] limits = {xMin, xMax, yMin, yMax, zMin, zMax};
            return limits;
        }
    }

    public static boolean isActiveOnAllSlices(RoiManager roiManager, Roi roi) {
        if (roi.getName() == null) {
            return false;
        }

        int sliceNumber = roiManager.getSliceNumber(roi.getName());
        return isActiveOnAllSlices(sliceNumber);
    }

    /**
     * Crops the given rectangle to the area [0, 0, width, height]
     *
     * @param bounds The rectangle to be fitted
     * @param width  Maximum width of the rectangle
     * @param height Maximum height of the rectangle
     * @return false if the height or width of the fitted rectangle is 0
     * (Couldn't be cropped inside the area).
     */
    public static boolean getSafeRoiBounds(Rectangle bounds, int width, int height) {
        int xMin = clamp(bounds.x, 0, width);
        int xMax = clamp(bounds.x + bounds.width, 0, width);
        int yMin = clamp(bounds.y, 0, height);
        int yMax = clamp(bounds.y + bounds.height, 0, height);
        int newWidth = xMax - xMin;
        int newHeight = yMax - yMin;

        bounds.setBounds(xMin, yMin, newWidth, newHeight);

        return newWidth > 0 && newHeight > 0;
    }

    /**
     * Crop a stack to the limits of the ROIs in the ROI Manager and optionally
     * fill the background with a single pixel value.
     *
     * @param roiMan         ROI Manager containing ROIs
     * @param sourceStack    input stack
     * @param fillBackground if true, background will be set to value.
     *                       NB Will not set the background of the areas copied from sourceStack
     * @param fillColor      value to set background to
     * @param padding        empty pixels to pad faces of cropped stack with
     * @return cropped copy of input stack
     */
    public static ImageStack cropStack(RoiManager roiMan, ImageStack sourceStack,
                                       boolean fillBackground, int fillColor, int padding) {
        int[] limits = getLimits(roiMan, sourceStack);

        if (limits == null) {
            return null;
        }

        final int xMin = limits[0];
        final int xMax = limits[1];
        final int yMin = limits[2];
        final int yMax = limits[3];
        final int zMin = limits[4];
        final int zMax = limits[5];

        final int croppedWidth = xMax - xMin + 2 * padding;
        final int croppedHeight = yMax - yMin + 2 * padding;

        ImageStack targetStack = new ImageStack(croppedWidth, croppedHeight);

        // copy
        ImageProcessor sourceProcessor;
        ImageProcessor targetProcessor;
        ArrayList<Roi> sliceRois;

        for (int sourceZ = zMin; sourceZ <= zMax; sourceZ++) {
            sliceRois = getSliceRoi(roiMan, sourceStack, sourceZ);
            if (sliceRois.size() == 0) {
                continue;
            }

            sourceProcessor = sourceStack.getProcessor(sourceZ);
            targetProcessor = sourceProcessor.createProcessor(croppedWidth, croppedHeight);

            if (fillBackground) {
                targetProcessor.setColor(fillColor);
                targetProcessor.fill();
            }

            copySlice(sourceProcessor, targetProcessor, sliceRois, padding);
            targetStack.addSlice("", targetProcessor);
        }

        // z padding
        targetProcessor = targetStack.getProcessor(1).createProcessor(croppedWidth, croppedHeight);
        if (fillBackground) {
            targetProcessor.setColor(fillColor);
            targetProcessor.fill();
        }
        for (int i = 0; i < padding; i++) {
            targetStack.addSlice("", targetProcessor.duplicate(), 0);
            targetStack.addSlice(targetProcessor.duplicate());
        }

        return targetStack;
    }

    /**
     * Remove all ROIs from the ROI manager
     */
    public static void deleteAll(RoiManager roiMan) {
        Roi[] rois = roiMan.getRoisAsArray();
        for (int i = 0; i < rois.length; i++) {
            if (roiMan.getCount() == 0)
                break;
            roiMan.select(i);
            roiMan.runCommand("delete");
        }
    }

    /**
     * Copies the pixels in the given ROI from the source image to the target
     * image. Copies only those pixels where the color of the given mask > 0.
     *
     * @param sourceProcessor Copy source
     * @param targetProcessor Copy target
     * @param minX            Horizontal start of the copy area 0 <= minX < width
     * @param minY            Vertical start of the copy area 0 <= minY < height
     * @param maxX            Horizontal end of the copy area 0 <= maxX <= width
     * @param maxY            Vertical end of the copy area 0 <= maxY <= height
     * @param padding         Number pixels added to each side of the copy target
     */
    public static void copyRoiWithMask(ImageProcessor sourceProcessor, ImageProcessor targetProcessor, final int minX,
                                       final int minY, final int maxX, final int maxY, final int padding) {
        ImageProcessor mask = sourceProcessor.getMask();
        if (mask == null) {
            copyRoi(sourceProcessor, targetProcessor, minX, minY, maxX, maxY, padding);
            return;
        }

        int targetY = padding;
        for (int sourceY = minY; sourceY < maxY; sourceY++) {
            int targetX = padding;
            for (int sourceX = minX; sourceX < maxX; sourceX++) {
                int maskColor = mask.get(sourceX, sourceY);
                if (maskColor > 0) {
                    int sourceColor = sourceProcessor.get(sourceX, sourceY);
                    targetProcessor.set(targetX, targetY, sourceColor);
                }
                targetX++;
            }
            targetY++;
        }
    }

    private static int clamp(int value, int min, int max) {
        if (Integer.compare(value, min) < 0) {
            return min;
        }
        if (Integer.compare(value, max) > 0) {
            return max;
        }
        return value;
    }

    private static boolean isActiveOnAllSlices(int sliceNumber) {
        return sliceNumber == NO_SLICE_NUMBER;
    }

    /**
     * Copies pixels under all the ROIs on a slide
     *
     * @param sourceProcessor The source image slide
     * @param targetProcessor The target slide
     * @param sliceRois       List of all the ROIs on the source slide
     * @param padding         Number of pixels added on each side of the target slide
     */
    private static void copySlice(ImageProcessor sourceProcessor, ImageProcessor targetProcessor,
                                  ArrayList<Roi> sliceRois, int padding) {
        for (Roi sliceRoi : sliceRois) {
            Rectangle rectangle = sliceRoi.getBounds();
            boolean valid = getSafeRoiBounds(rectangle, sourceProcessor.getWidth(), sourceProcessor.getHeight());

            if (!valid) {
                continue;
            }

            int minY = rectangle.y;
            int minX = rectangle.x;
            int maxY = rectangle.y + rectangle.height;
            int maxX = rectangle.x + rectangle.width;

            ImageProcessor mask = sourceProcessor.getMask();
            if (mask == null) {
                copyRoi(sourceProcessor, targetProcessor, minX, minY, maxX, maxY, padding);
            } else {
                copyRoiWithMask(sourceProcessor, targetProcessor, minX, minY, maxX, maxY, padding);
            }
        }
    }

    /**
     * Copies the pixels in the given ROI from the source image to the target
     * image.
     *
     * @param sourceProcessor Copy source
     * @param targetProcessor Copy target
     * @param minX            Horizontal start of the copy area 0 <= minX < width
     * @param minY            Vertical start of the copy area 0 <= minY < height
     * @param maxX            Horizontal end of the copy area 0 <= maxX <= width
     * @param maxY            Vertical end of the copy area 0 <= maxY <= height
     * @param padding         Number pixels added to each side of the copy target
     */
    private static void copyRoi(ImageProcessor sourceProcessor, ImageProcessor targetProcessor, final int minX,
                                final int minY, final int maxX, final int maxY, final int padding) {
        int targetY = padding;
        for (int sourceY = minY; sourceY < maxY; sourceY++) {
            int targetX = padding;
            for (int sourceX = minX; sourceX < maxX; sourceX++) {
                int sourceColor = sourceProcessor.get(sourceX, sourceY);
                targetProcessor.set(targetX, targetY, sourceColor);
                targetX++;
            }
            targetY++;
        }
    }
}
