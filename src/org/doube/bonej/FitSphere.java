package org.doube.bonej;

/**
 * FitSphere class for ImageJ
 * Copyright 2009 Michael Doube
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
import ij.gui.Roi;
import ij.measure.Calibration;
import ij.plugin.frame.RoiManager;
import ij.IJ;
import ij.ImagePlus;

import java.awt.List;
import java.awt.Rectangle;

import org.doube.jama.Matrix;

/**
 *<p>
 * Takes point selections from ROI manager and returns the centroid and radius
 * of a best fit sphere Ported from Angelo Tardugno's C++
 * </p>
 * 
 * 
 *@author Michael Doube and Angelo Tardugno
 *@version 0.1
 */
public class FitSphere {

    /**
     * Get the calibrated 3D coordinates of point ROIs from the ROI manager
     * 
     * @param imp
     * @param roiMan
     * @return double[n][3] containing n (x, y, z) coordinates
     */
    public double[][] getRoiManPoints(ImagePlus imp, RoiManager roiMan) {
	Calibration cal = imp.getCalibration();
	double vW = cal.pixelWidth;
	double vH = cal.pixelHeight;
	double vD = cal.pixelDepth;
	int nPoints = 0;
	List listRoi = roiMan.getList();
	Roi[] roiList = roiMan.getRoisAsArray();
	for (int i = 0; i < roiMan.getCount(); i++) {
	    Roi roi = roiList[i];
	    if (roi.getType() == 10) {
		nPoints++;
	    }
	}
	double[][] dataPoints = new double[nPoints][3];
	int j = 0;
	for (int i = 0; i < roiMan.getCount(); i++) {
	    Roi roi = roiList[i];
	    if (roi.getType() == 10) {
		String label = listRoi.getItem(i);
		Rectangle xy = roi.getBoundingRect();
		dataPoints[j][0] = xy.getX() * vW;
		dataPoints[j][1] = xy.getY() * vH;
		dataPoints[j][2] = roiMan.getSliceNumber(label) * vD;
		j++;
	    }
	}
	return dataPoints;
    }

    /**
     * Fit a sphere to 3D coordinates
     * 
     * @param points
     *            double[n][3] containing n (x, y, z) coordinates
     * @return double[4] containing (x, y, z) centre and radius
     */
    public double[] fitSphere(double[][] points) {
	int nPoints = points.length;
	if (nPoints < 5) {
	    IJ.error("Too few points to calculate a sphere");
	    double[] error = { -1, -1, -1, -1 };
	    return error;
	}
	double xSum = 0, ySum = 0, zSum = 0;
	for (int i = 0; i < points.length; i++) {
	    xSum += points[i][0];
	    ySum += points[i][1];
	    zSum += points[i][2];
	}

	double x = xSum / nPoints;
	double y = ySum / nPoints;
	double z = zSum / nPoints;

	double[] radii = new double[nPoints];
	double g_new = 100.0;
	double g_old = 1.0;
	double r = 0;

	for (int i = 0; i < nPoints; i++) {
	    double pX = points[i][0] - x;
	    double pY = points[i][1] - y;
	    double pZ = points[i][2] - z;
	    r += Math.sqrt(pX * pX + pY * pY + pZ * pZ);
	}
	r /= nPoints;

	while (Math.abs(g_new - g_old) > 1e-10) {
	    Matrix J = new Matrix(nPoints, 4);
	    double[][] Jp = J.getArray();
	    Matrix D = new Matrix(nPoints, 1);
	    double[][] dp = D.getArray(); // dp is a pointer to d's values
	    g_old = g_new;
	    for (int i = 0; i < nPoints; i++) {
		double pX = points[i][0] - x;
		double pY = points[i][1] - y;
		double pZ = points[i][2] - z;
		radii[i] = Math.sqrt(pX * pX + pY * pY + pZ * pZ);
		dp[i][0] = radii[i] - r;
		Jp[i][0] = -pX / radii[i];
		Jp[i][1] = -pY / radii[i];
		Jp[i][2] = -pZ / radii[i];
		Jp[i][3] = -1;
	    }
	    D = D.times(-1);
	    Matrix J1 = J;
	    J = J.transpose();
	    Matrix J2 = J.times(J1);
	    Matrix Jd = J.times(D);
	    Matrix X = J2.inverse().times(Jd);
	    double[][] xp = X.getArray();
	    x += xp[0][0];
	    y += xp[1][0];
	    z += xp[2][0];
	    r += xp[3][0];
	    D = D.times(-1);
	    Matrix G = J.times(D);
	    double[][] Gp = G.getArray();
	    g_new = 0.0;
	    for (int i = 0; i < 4; i++)
		g_new += Gp[i][0];
	}
	double[] centreRadius = {x, y, z, r};
	return centreRadius;
    }
}
