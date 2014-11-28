package org.doube.geometry;

import ij.IJ;

import java.util.Arrays;

import org.doube.jama.EigenvalueDecomposition;
import org.doube.jama.Matrix;
import org.doube.util.ArrayHelper;

/**
 * <p>
 * Represents an ellipsoid defined as <i>ax</i><sup>2</sup> +
 * <i>by</i><sup>2</sup> + <i>cz</i><sup>2</sup> + 2<i>dxy</i> + 2<i>exz</i> +
 * 2<i>fyz</i> + 2<i>gx</i> + 2<i>hy</i> + 2<i>iz</i> = 1 <br />
 * </p>
 * 
 * @author Michael Doube
 */
public class Ellipsoid {

	// centre
	private double cx;
	private double cy;
	private double cz;

	// radii
	private double ra;
	private double rb;
	private double rc;

	// ellipsoid equation coefficients
	private double a;
	private double b;
	private double c;
	private double d;
	private double e;
	private double f;
	private double g;
	private double h;
	private double i;

	// calculated volume
	private double volume;

	/** Eigenvector matrix */
	// private Matrix V;
	private double[][] ev;

	/** Eigenvalue matrix */
	// private Matrix D;
	private double[][] ed;

	/** 3x3 matrix describing shape of ellipsoid */
	// private Matrix H;
	private double[][] eh;

	/**
	 * Instantiate an ellipsoid from the result of FitEllipsoid
	 * 
	 * @param ellipsoid
	 */
	public Ellipsoid(Object[] ellipsoid) {
		// Object[] ellipsoid = { centre, radii, eigenVectors, equation, E };

		double[] centre = (double[]) ellipsoid[0];
		this.cx = centre[0];
		this.cy = centre[1];
		this.cz = centre[2];

		double[] radii = (double[]) ellipsoid[1];
		this.ra = radii[0];
		this.rb = radii[1];
		this.rc = radii[2];

		if (Double.isNaN(ra) || Double.isNaN(rb) || Double.isNaN(rc))
			throw new IllegalArgumentException("Radius is NaN");

		// this.V = new Matrix(3, 3);
		this.ev = new double[3][3];
		// this.D = new Matrix(3, 3);
		this.ed = new double[3][3];
		// this.H = new Matrix(3, 3);
		this.eh = new double[3][3];
		setRotation((double[][]) ellipsoid[2]);
		setEigenvalues();
		setVolume();

		double[] equation = (double[]) ellipsoid[3];
		this.a = equation[0];
		this.b = equation[1];
		this.c = equation[2];
		this.d = equation[3];
		this.e = equation[4];
		this.f = equation[5];
		this.g = equation[6];
		this.h = equation[7];
		this.i = equation[8];
	}

	/**
	 * Construct an Ellipsoid from the radii (a,b,c), centroid (cx, cy, cz) and
	 * Eigenvectors.
	 * 
	 * @param a
	 * @param b
	 * @param c
	 * @param cx
	 * @param cy
	 * @param cz
	 * @param eigenVectors
	 */
	public Ellipsoid(double a, double b, double c, double cx, double cy,
			double cz, double[][] eigenVectors) {

		this.ra = a;
		this.rb = b;
		this.rc = c;
		this.cx = cx;
		this.cy = cy;
		this.cz = cz;
		// this.V = new Matrix(3, 3);
		this.ev = new double[3][3];
		// this.D = new Matrix(3, 3);
		this.ed = new double[3][3];
		// this.H = new Matrix(3, 3);
		this.eh = new double[3][3];
		setRotation(eigenVectors);
		setEigenvalues();
		// TODO update equation variables
		setVolume();
	}

	public double getVolume() {
		return (new Double(volume)).doubleValue();
	}

	private void setVolume() {
		volume = Math.PI * ra * rb * rc * 4 / 3;
	}

	/**
	 * @return the semiaxis lengths a, b and c. Note these are not ordered by
	 *         size, but the order does relate to the 0th, 1st and 2nd columns
	 *         of the rotation matrix respectively.
	 */
	public double[] getRadii() {
		double[] radii = { ra, rb, rc };
		return radii.clone();
	}

	public double[] getEquation() {
		double[] equation = { a, b, c, d, e, f, g, h, i };
		return equation.clone();
	}

	/**
	 * Method based on the inequality
	 * 
	 * (X-X0)H(X-X0)^T <= 1
	 * 
	 * Where X is the test point, X0 is the centroid, H is the ellipsoid's 3x3
	 * matrix
	 * 
	 * @param x
	 * @param y
	 * @param z
	 * @return true if the point (x,y,z) lies inside or on the ellipsoid, false
	 *         otherwise
	 */
	public boolean contains(double x, double y, double z) {
		// calculate vector between point and centroid
		final double vx = x - cx;
		final double vy = y - cy;
		final double vz = z - cz;

		// calculate distance from centroid
		final double length = Math.sqrt(vx * vx + vy * vy + vz * vz);

		double[] radii = { ra, rb, rc };
		Arrays.sort(radii);

		// if further from centroid than major semiaxis length
		// must be outside
		if (length > radii[2])
			return false;

		// if length closer than minor semiaxis length
		// must be inside
		if (length <= radii[0])
			return true;

		// final double[][] h = H.getArray();
		final double[][] h = eh;

		final double dot0 = vx * h[0][0] + vy * h[1][0] + vz * h[2][0];
		final double dot1 = vx * h[0][1] + vy * h[1][1] + vz * h[2][1];
		final double dot2 = vx * h[0][2] + vy * h[1][2] + vz * h[2][2];

		final double dot = dot0 * vx + dot1 * vy + dot2 * vz;

		if (dot <= 1)
			return true;

		return false;
	}

	public double solve(double x, double y, double z) {
		return a * x * x + b * y * y + c * z * z + 2
				* (d * x * y + e * x * z + f * y * z + g * x + h * y + i * z);
	}

	public double[] getCentre() {
		double[] centre = { cx, cy, cz };
		return centre.clone();
	}

	public double[][] getSurfacePoints(final int nPoints) {

		// get regularly-spaced points on the unit sphere
		double[][] vectors = Vectors.regularVectors(nPoints);
		return getSurfacePoints(vectors);

	}

	public double[][] getSurfacePoints(final double[][] vectors) {
		final int nPoints = vectors.length;
		for (int p = 0; p < nPoints; p++) {
			// stretch the unit sphere into an ellipsoid
			final double x = ra * vectors[p][0];
			final double y = rb * vectors[p][1];
			final double z = rc * vectors[p][2];
			// rotate and translate the ellipsoid into position
			vectors[p][0] = x * ev[0][0] + y * ev[0][1] + z * ev[0][2] + cx;
			vectors[p][1] = x * ev[1][0] + y * ev[1][1] + z * ev[1][2] + cy;
			vectors[p][2] = x * ev[2][0] + y * ev[2][1] + z * ev[2][2] + cz;
		}
		return vectors;
	}

	/**
	 * Dilate all three axes by a fractional increment
	 * 
	 * @param increment
	 */
	public void dilate(double increment) {
		dilate(this.ra * increment, this.rb * increment, this.rc * increment);
	}

	/**
	 * Constrict all three axes by a fractional increment
	 * 
	 * @param increment
	 */
	public void contract(double increment) {
		dilate(-increment);
	}

	/**
	 * Contract the semiaxes by independent absolute amounts
	 * 
	 * @param ca
	 * @param cb
	 * @param cd
	 */
	public void contract(double ca, double cb, double cd) {
		dilate(-ca, -cb, -cd);
	}

	/**
	 * Dilate the ellipsoid semiaxes by independent absolute amounts
	 * 
	 * @param da
	 * @param db
	 * @param dc
	 */
	public void dilate(double da, double db, double dc) {
		setRadii(this.ra + da, this.rb + db, this.rc + dc);
	}

	/**
	 * Translate the ellipsoid
	 * 
	 * @param dx
	 *            shift in x
	 * @param dy
	 *            shift in y
	 * @param dz
	 *            shift in z
	 */
	public void translate(double dx, double dy, double dz) {
		this.cx += dx;
		this.cy += dy;
		this.cz += dz;
	}

	/**
	 * Translate the ellipsoid to a given new centroid
	 * 
	 * @param x
	 *            new centroid x-coordinate
	 * @param y
	 *            new centroid y-coordinate
	 * @param z
	 *            new centroid z-coordinate
	 */
	public void setCentroid(double x, double y, double z) {
		this.cx = x;
		this.cy = y;
		this.cz = z;
	}

	/**
	 * Rotate the ellipsoid by the given 3x3 Matrix
	 * 
	 * @param R
	 *            a 3x3 rotation matrix
	 */
	// public void rotate(Matrix R) {
	// setRotation(this.V.times(R));
	// }

	/**
	 * Rotate the ellipsoid by the given 3x3 Matrix
	 * 
	 * @param R
	 *            a 3x3 rotation matrix
	 */
	public void rotate(double[][] rotation) {
		setRotation(ArrayHelper.times(this.ev, rotation));
	}

	/**
	 * Set the rotation to the supplied eigenvector matrix
	 * 
	 * @param R
	 *            3x3 eigenvector matrix
	 */
	// public void setRotation(Matrix R) {
	// final double detDiff = Math.abs(1 - R.det());
	// if (!is3x3Matrix(R) || detDiff > 1E-10)
	// throw new IllegalArgumentException("Not a 3x3 rotation matrix");
	// this.V = R.copy();
	// update3x3Matrix();
	// }

	/**
	 * Faster version which avoids Matrix instantiation and which does no error
	 * checking
	 */
	public void setRotation(double[][] rotation) {
		this.ev = rotation.clone();
		update3x3Matrix();
	}

	/**
	 * Calculate the intercepts of the x, y and z axes
	 * 
	 * @return array containing 6 intercepts, ordered +x, -x, +y, -y, +z, -z
	 */
	public double[] intercepts() {
		double plusX = (-g + Math.sqrt(g * g + 4 * a)) / (2 * a);
		double minusX = (-g - Math.sqrt(g * g + 4 * a)) / (2 * a);
		double plusY = (-h + Math.sqrt(h * h + 4 * b)) / (2 * b);
		double minusY = (-h - Math.sqrt(h * h + 4 * b)) / (2 * b);
		double plusZ = (-i + Math.sqrt(i * i + 4 * c)) / (2 * c);
		double minusZ = (-i - Math.sqrt(i * i + 4 * c)) / (2 * c);
		return new double[] { plusX, minusX, plusY, minusY, plusZ, minusZ };
	}

	/**
	 * Set the radii (semiaxes). No ordering is assumed, except with regard to
	 * the columns of the eigenvector rotation matrix (i.e. a relates to the 0th
	 * eigenvector column, b to the 1st and c to the 2nd)
	 * 
	 * @param a
	 * @param b
	 * @param c
	 */
	public void setRadii(double a, double b, double c) {
		if (a <= 0 || b <= 0 || c <= 0) {
			throw new IllegalArgumentException(
					"Ellipsoid cannot have semiaxis <= 0");
		}
		this.ra = a;
		this.rb = b;
		this.rc = c;
		setEigenvalues();
		setVolume();
	}

	/**
	 * Calculates eigenvalues from current radii
	 */
	private void setEigenvalues() {
		// this.D.set(0, 0, 1 / (this.ra * this.ra));
		// this.D.set(1, 1, 1 / (this.rb * this.rb));
		// this.D.set(2, 2, 1 / (this.rc * this.rc));
		this.ed[0][0] = 1 / (this.ra * this.ra);
		this.ed[1][1] = 1 / (this.rb * this.rb);
		this.ed[2][2] = 1 / (this.rc * this.rc);
		update3x3Matrix();
	}

	/**
	 * Needs to be run any time the eigenvalues or eigenvectors change
	 */
	private void update3x3Matrix() {
		// this.H = (this.V.times(this.D)).times(this.V.transpose());
		this.eh = ArrayHelper.times(ArrayHelper.times(ev, ed),
				ArrayHelper.transpose(ev));
	}

	// private boolean is3x3Matrix(Matrix rotation) {
	// return (rotation.getRowDimension() == 3 && rotation
	// .getColumnDimension() == 3);
	// }

	/**
	 * Calculate the matrix representation of the ellipsoid (centre,
	 * eigenvalues, eigenvectors) from the equation variables
	 * <i>ax</i><sup>2</sup> + <i>by</i><sup>2</sup> + <i>cz</i><sup>2</sup> +
	 * 2<i>dxy</i> + 2<i>exz</i> + 2<i>fyz</i> + 2<i>gx</i> + 2<i>hy</i> +
	 * 2<i>iz</i> = 1 <br />
	 * 
	 * @param a
	 * @param b
	 * @param c
	 * @param d
	 * @param e
	 * @param f
	 * @param g
	 * @param h
	 * @param i
	 * @return Object[] array containing centre (double[3]), eigenvalues
	 *         (double[3][3]), eigenvectors (double[3][3]), and the
	 *         EigenvalueDecomposition
	 */
	public static Object[] matrixFromEquation(double a, double b, double c,
			double d, double e, double f, double g, double h, double i) {

		// the fitted equation
		double[][] v = { { a }, { b }, { c }, { d }, { e }, { f }, { g },
				{ h }, { i } };
		Matrix V = new Matrix(v);

		// 4x4 based on equation variables
		double[][] aa = { { a, d, e, g }, { d, b, f, h }, { e, f, c, i },
				{ g, h, i, -1 }, };
		Matrix A = new Matrix(aa);

		// find the centre
		Matrix C = (A.getMatrix(0, 2, 0, 2).times(-1).inverse()).times(V
				.getMatrix(6, 8, 0, 0));

		// using the centre and 4x4 calculate the
		// eigendecomposition
		Matrix T = Matrix.eye(4);
		T.setMatrix(3, 3, 0, 2, C.transpose());
		Matrix R = T.times(A.times(T.transpose()));
		double r33 = R.get(3, 3);
		Matrix R02 = R.getMatrix(0, 2, 0, 2);
		EigenvalueDecomposition E = new EigenvalueDecomposition(R02.times(-1
				/ r33));

		double[] centre = C.getColumnPackedCopy();
		double[][] eigenVectors = E.getV().getArrayCopy();
		double[][] eigenValues = E.getD().getArrayCopy();
		Object[] result = { centre, eigenValues, eigenVectors, E };
		return result;
	}

	/**
	 * Perform a deep copy of this Ellipsoid
	 * 
	 * @return
	 */
	public Ellipsoid copy() {
		Ellipsoid copy = new Ellipsoid(this.ra, this.rb, this.rc, this.cx,
				this.cy, this.cz, this.ev);
		return copy;
	}

	/**
	 * 
	 * @param centre
	 * @param eigenValues
	 * @param eigenVectors
	 * @return
	 */
	public static double[] equationFromMatrix(double[] centre,
			double[][] eigenValues, double[][] eigenVectors) {

		// orientation of ellipsoid
		Matrix P = new Matrix(eigenVectors);

		// size of ellipsoid
		// related to radii by r = sqrt(1/eVal)
		Matrix D = new Matrix(eigenValues);

		// B = P^-1DP
		// where B is a square matrix, D is eigenvalues, P is eigenvectors
		Matrix B = (P.inverse().times(D)).times(P);

		// now B = R02.times(-1/r33) in the above equation

		// make a 4x4 matrix with B in the top left and 1 in the bottom right, 0
		// elsewhere indicating no translation
		Matrix E = Matrix.eye(4);
		E.setMatrix(0, 2, 0, 2, B);

		// now set up a 4x4 translation matrix
		Matrix C = new Matrix(3, 1);
		C.set(0, 0, centre[0]);
		C.set(1, 0, centre[1]);
		C.set(2, 0, centre[2]);
		Matrix T = Matrix.eye(4);
		T.setMatrix(0, 2, 3, 3, C);

		// above leaves bottom row zeros

		// work out the translated ellipsoid
		E = E.times(T);

		// get the negative inverse of the bottom right corner
		final double e33 = -1 / E.get(3, 3);

		// work out the scaled ellipsoid
		E = E.times(e33);

		// pack the matrix into the equation form
		double[] e = E.getColumnPackedCopy();
		double[] equation = { e[0], e[5], e[10], e[1], e[2], e[6], e[3], e[7],
				e[11] };

		return equation;
	}

	/**
	 * Generate a string of useful information about this Ellipsoid
	 * 
	 * @return
	 */
	public String debugOutput() {
		String string = "Ellipsoid variables:\n";
		string = string + "a = " + a + "\n";
		string = string + "b = " + b + "\n";
		string = string + "c = " + c + "\n";
		string = string + "d = " + d + "\n";
		string = string + "e = " + e + "\n";
		string = string + "f = " + f + "\n";
		string = string + "g = " + g + "\n";
		string = string + "h = " + h + "\n";
		string = string + "i = " + i + "\n";

		string = string + "\nCentre: \n";
		string = string + "cx: " + this.cx + "\n";
		string = string + "cy: " + this.cy + "\n";
		string = string + "cz: " + this.cz + "\n";

		string = string + "\nRadii: \n";
		string = string + "ra: " + this.ra + "\n";
		string = string + "rb: " + this.rb + "\n";
		string = string + "rc: " + this.rc + "\n";

		string = string + "\nEigenvalues: \n";
		string = string + "eVal0 = " + ed[0][0] + "\n";
		string = string + "eVal1 = " + ed[1][1] + "\n";
		string = string + "eVal2 = " + ed[2][2] + "\n";

		string = string + "\nEigenvectors: \n";
		string = string + "eV00 = " + ev[0][0] + "\n";
		string = string + "eV01 = " + ev[0][1] + "\n";
		string = string + "eV02 = " + ev[0][2] + "\n";
		string = string + "eV10 = " + ev[1][0] + "\n";
		string = string + "eV11 = " + ev[1][1] + "\n";
		string = string + "eV12 = " + ev[1][2] + "\n";
		string = string + "eV20 = " + ev[2][0] + "\n";
		string = string + "eV21 = " + ev[2][1] + "\n";
		string = string + "eV22 = " + ev[2][2] + "\n";
		return string;
	}
}
