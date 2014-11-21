package org.doube.geometry;

import java.util.Arrays;

import org.doube.jama.EigenvalueDecomposition;
import org.doube.jama.Matrix;

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

	// eigenvectors of axes
	private double[][] eigenVectors;
	// unpacked for convenience
	private double eV00;
	private double eV01;
	private double eV02;
	private double eV10;
	private double eV11;
	private double eV12;
	private double eV20;
	private double eV21;
	private double eV22;

	/** Eigenvector matrix */
	private Matrix V;

	/** Eigenvalue matrix */
	private Matrix D;

	/** 3x3 matrix describing shape of ellipsoid */
	private Matrix H;

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

		double[][] eigenVectors = (double[][]) ellipsoid[2];
		this.V = new Matrix(3, 3);
		this.D = new Matrix(3, 3);
		this.H = new Matrix(3, 3);
		setEigenVectors(eigenVectors);
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
	 * Eigenvectors. Assumes that a > b > c and enforces it with a sort().
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

		double[] radii = { a, b, c };
		Arrays.sort(radii);

		this.ra = radii[2];
		this.rb = radii[1];
		this.rc = radii[0];
		this.cx = cx;
		this.cy = cy;
		this.cz = cz;
		this.V = new Matrix(3, 3);
		this.D = new Matrix(3, 3);
		this.H = new Matrix(3, 3);
		setEigenVectors(eigenVectors);
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

	public double[] getRadii() {
		double[] radii = { ra, rb, rc };
		return radii.clone();
	}

	public double[] getEquation() {
		double[] equation = { a, b, c, d, e, f, g, h, i };
		return equation.clone();
	}

	public boolean contains(double x, double y, double z) {

		// calculate vector between point and centroid
		double vx = x - cx;
		double vy = y - cy;
		double vz = z - cz;

		// calculate distance from centroid
		final double length = Math.sqrt(vx * vx + vy * vy + vz * vz);

		// if further from centroid than major semiaxis length
		// must be outside
		if (length > ra)
			return false;

		// if length closer than minor semiaxis length
		// must be inside
		if (length <= rc)
			return true;

		// calculate unit vector (normalise)
		vx /= length;
		vy /= length;
		vz /= length;
		System.out.println("Unit vector is [" + vx + ", " + vy + ", " + vz
				+ "]");

		// get eigenvector matrix
		Matrix eV = new Matrix(eigenVectors);
		eV.print(8, 5);
		// invert it
		Matrix eVinv = eV.inverse();
		eVinv.print(8, 5);
		double[][] dv = eVinv.getArrayCopy();
		// calculate the derotated unit vector
		double dx = vx * dv[0][0] + vy * dv[0][1] + vz * dv[0][2];
		double dy = vx * dv[1][0] + vy * dv[1][1] + vz * dv[1][2];
		double dz = vx * dv[2][0] + vy * dv[2][1] + vz * dv[2][2];
		System.out.println("Derotated vector is [" + dx + ", " + dy + ", " + dz
				+ "]");

		// find the size of the ellipsoid in this direction using semiaxis
		// lengths
		dx = dx * ra;
		dy = dy * rb;
		dz = dz * rc;

		System.out.println("Ellipsoid point is (" + dx + ", " + dy + ", " + dz
				+ ")");

		// returns true if the ellipsoid is bigger in this direction
		// than the test point
		double ellipsoidLength = Math.sqrt(dx * dx + dy * dy + dz * dz);
		System.out.println("Point length = " + length + ", ellipsoid length = "
				+ ellipsoidLength + ", ratio = " + ellipsoidLength / length);
		return (ellipsoidLength > length);

	}

	/**
	 * Alternative contains() method based on the inequality
	 * 
	 * (X-X0)'H(X-X0) <= 1
	 * 
	 * Where X is the test point, X0 is the centroid, H is the ellipsoid's 3x3
	 * matrix
	 * 
	 * @param x
	 * @param y
	 * @param z
	 * @param cx
	 * @param cy
	 * @param cz
	 * @param H
	 * @return
	 */
	public static boolean contains(double x, double y, double z, double cx,
			double cy, double cz, Matrix H) {

		Matrix X = vector3(x, y, z);
		Matrix X0 = vector3(cx, cy, cz);

		Matrix XX0 = X.minus(X0);
		if (XX0.transpose().times(H).times(XX0).get(0, 0) <= 1)
			return true;

		return false;
	}

	public double solve(double x, double y, double z) {
		return a * x * x + b * y * y + c * z * z + 2
				* (d * x * y + e * x * z + f * y * z + g * x + h * y + i * z);
	}

	public double getMajorRadius() {
		double val = ra;
		return val;
	}

	public double[] getCentre() {
		double[] centre = { cx, cy, cz };
		return centre.clone();
	}

	public double[][] getSurfacePoints(final int nPoints) {

		// get regularly-spaced points on the unit sphere
		double[][] points = Vectors.regularVectors(nPoints);

		for (int p = 0; p < nPoints; p++) {
			// stretch the unit sphere into an ellipsoid
			final double x = ra * points[p][0];
			final double y = rb * points[p][1];
			final double z = rc * points[p][2];
			// rotate and translate the ellipsoid into position
			points[p][0] = x * eV00 + y * eV01 + z * eV02 + cx;
			points[p][1] = x * eV10 + y * eV11 + z * eV12 + cy;
			points[p][2] = x * eV20 + y * eV21 + z * eV22 + cz;
		}
		return points;
	}

	private void setEigenVectors(double[][] eigenVectors) {
		this.eigenVectors = eigenVectors;
		this.eV00 = this.eigenVectors[0][0];
		this.eV01 = this.eigenVectors[0][1];
		this.eV02 = this.eigenVectors[0][2];
		this.eV10 = this.eigenVectors[1][0];
		this.eV11 = this.eigenVectors[1][1];
		this.eV12 = this.eigenVectors[1][2];
		this.eV20 = this.eigenVectors[2][0];
		this.eV21 = this.eigenVectors[2][1];
		this.eV22 = this.eigenVectors[2][2];
		this.V.set(0, 0, this.eV00);
		this.V.set(0, 1, this.eV01);
		this.V.set(0, 2, this.eV02);
		this.V.set(1, 0, this.eV10);
		this.V.set(1, 1, this.eV11);
		this.V.set(1, 2, this.eV12);
		this.V.set(2, 0, this.eV20);
		this.V.set(2, 1, this.eV21);
		this.V.set(2, 2, this.eV22);
		update3x3Matrix();
	}

	/**
	 * Dilate all three axes by an increment
	 * 
	 * @param increment
	 */
	public void dilate(double increment) {
		this.ra += increment;
		this.rb += increment;
		this.rc += increment;

		setVolume();

		setEigenvalues();
	}

	/**
	 * Constrict all three axes by an increment
	 * 
	 * @param increment
	 */
	public void contract(double increment) {
		dilate(-increment);
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
	 * @param rotation
	 *            a 3x3 rotation matrix
	 */
	public void rotate(Matrix rotation) {
		if (!is3x3Matrix(rotation))
			throw new IllegalArgumentException("Not a 3x3 rotation matrix");
		setRotation(this.V.times(rotation));
	}

	/**
	 * Set the rotation to the supplied eigenvector matrix
	 * 
	 * @param rotation
	 *            3x3 eigenvector matrix
	 */
	public void setRotation(Matrix rotation) {
		if (!is3x3Matrix(rotation))
			throw new IllegalArgumentException("Not a 3x3 rotation matrix");
		setEigenVectors(rotation.getArray());
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
	 * Calculates eigenvalues from current radii
	 */
	private void setEigenvalues() {
		this.D.set(0, 0, 1 / (this.ra * this.ra));
		this.D.set(1, 1, 1 / (this.rb * this.rb));
		this.D.set(2, 2, 1 / (this.rc * this.rc));
		update3x3Matrix();
	}

	/**
	 * Needs to be run any time the eigenvalues or eigenvectors change
	 */
	private void update3x3Matrix() {
		this.H = (this.V.inverse().times(this.D)).times(this.V);
	}


	private boolean is3x3Matrix(Matrix rotation) {
		return (rotation.getRowDimension() == 3 && rotation
				.getColumnDimension() == 3);
	}

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
		Matrix C = vector3(centre[0], centre[1], centre[2]);
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

	private static Matrix vector3(double j, double k, double l) {
		Matrix C = new Matrix(3, 1);
		C.set(0, 0, j);
		C.set(1, 0, k);
		C.set(2, 0, l);
		return C;
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
		string = string + "eVal0 = " + D.get(0, 0) + "\n";
		string = string + "eVal1 = " + D.get(1, 1) + "\n";
		string = string + "eVal2 = " + D.get(2, 2) + "\n";

		string = string + "\nEigenvectors: \n";
		string = string + "eV00 = " + eV00 + "\n";
		string = string + "eV01 = " + eV01 + "\n";
		string = string + "eV02 = " + eV02 + "\n";
		string = string + "eV10 = " + eV10 + "\n";
		string = string + "eV11 = " + eV11 + "\n";
		string = string + "eV12 = " + eV12 + "\n";
		string = string + "eV20 = " + eV20 + "\n";
		string = string + "eV21 = " + eV21 + "\n";
		string = string + "eV22 = " + eV22 + "\n";
		return string;
	}

}
