package org.doube.geometry;

import ij.IJ;

import org.doube.jama.EigenvalueDecomposition;
import org.doube.jama.Matrix;

public class FitEllipsoid {
	/**
	 * Calculate the best-fit ellipsoid by least squares
	 * 
	 * @param coOrdinates
	 *            n x 3 array containing 3D coordinates
	 * @return 10 x 1 array containing constants for the ellipsoid equation
	 *         <i>ax</i><sup>2</sup> + <i>by</i><sup>2</sup> +
	 *         <i>cz</i><sup>2</sup> + 2<i>fyz</i> + 2<i>gxz</i> + 2<i>hxy</i> +
	 *         2<i>px</i> + 2<i>qy</i> + 2<i>rz</i> + <i>d</i> = 0
	 * 
	 * @see <p>
	 *      Qingde Li, Griffiths J (2004) Least squares ellipsoid specific
	 *      fitting. Geometric Modeling and Processing, 2004. Proceedings. pp.
	 *      335-340. <a
	 *      href="http://dx.doi.org/10.1109/GMAP.2004.1290055">doi:10.1109
	 *      /GMAP.2004.1290055</a>
	 *      </p>
	 */
	public double[] liGriffiths(double[][] coOrdinates, int maxPoints) {
		final int nPoints = coOrdinates.length;
		IJ.showStatus("Fitting ellipsoid");
		double[][] coOrd = new double[maxPoints][3];
		if (nPoints > maxPoints) {
			// randomly subsample
			for (int n = 0; n < maxPoints; n++) {
				final int point = (int) Math.floor(Math.random() * nPoints);
				coOrd[n][0] = coOrdinates[point][0];
			}
		}
		double[][] Darray = new double[maxPoints][10];
		for (int n = 0; n < maxPoints; n++) {
			// populate each column of D with ten-element Xi
			final double xi = coOrd[n][0];
			final double yi = coOrd[n][1];
			final double zi = coOrd[n][2];
			Darray[n][0] = xi * xi;
			Darray[n][1] = yi * yi;
			Darray[n][2] = zi * zi;
			Darray[n][3] = yi * zi * 2;
			Darray[n][4] = xi * zi * 2;
			Darray[n][5] = xi * yi * 2;
			Darray[n][6] = xi * 2;
			Darray[n][7] = yi * 2;
			Darray[n][8] = zi * 2;
			Darray[n][9] = 1;
		}
		Matrix D = new Matrix(Darray);
		Matrix S = D.times(D.transpose());
		// Define 6x6 Matrix C1 (Eq. 7)
		double[][] C1array = { { -1, 1, 1, 0, 0, 0 }, { 1, -1, 1, 0, 0, 0 },
				{ 1, 1, -1, 0, 0, 0 }, { 0, 0, 0, -4, 0, 0 },
				{ 0, 0, 0, 0, -4, 0 }, { 0, 0, 0, 0, 0, -4 } };
		Matrix C1 = new Matrix(C1array);
		// Define S11, S12, S22 from S (Eq. 11)
		Matrix S11 = S.getMatrix(0, 5, 0, 5);
		Matrix S12 = S.getMatrix(0, 5, 6, 9);
		Matrix S22 = S.getMatrix(6, 9, 6, 9);
		// Eq. 15
		Matrix C1inv = C1.inverse();
		Matrix S12T = S12.transpose();
		Matrix S22inv = S22.inverse();
		Matrix Eq15 = C1inv.times(S11.minus(S12.times(S22inv.times(S12T))));
		// However, the matrix (S11 - S12*S22^-1*S12^T ) can be singular in some
		// cases. In this situations, the corresponding u1 can be replaced
		// with the eigenvector associated with the largest eigenvalue.

		EigenvalueDecomposition E = new EigenvalueDecomposition(Eq15);
		Matrix eigenValues = E.getD();
		eigenValues.printToIJLog("eigenValues");
		Matrix eigenVectors = E.getV();
		eigenVectors.printToIJLog("eigenVectors");
		Matrix v1 = new Matrix(6, 1);
		double[][] EigenValueMatrix = eigenValues.getArray();
		double posVal = 1 - Double.MAX_VALUE;
		final int evl = EigenValueMatrix.length;
		for (int p = 0; p < evl; p++) {
			if (EigenValueMatrix[p][p] > posVal) {
				posVal = EigenValueMatrix[p][p];
				v1 = eigenVectors.getMatrix(0, 5, p, p);
			}
		}
		v1.printToIJLog("v1");
		Matrix v2 = S22.inverse().times(S12.transpose()).times(v1).times(-1);
		v2.printToIJLog("v2");
		Matrix v = new Matrix(10, 1);
		int[] c = { 0 };
		v.setMatrix(0, 5, c, v1);
		v.setMatrix(6, 9, c, v2);
		double[] ellipsoid = v.getColumnPackedCopy();
		IJ.log("Ellipsoid equation: " + ellipsoid[0] + " x^2 + " + ellipsoid[1]
				+ " y^2 + " + ellipsoid[2] + " z^2 + " + ellipsoid[3]
				+ " 2yz + " + ellipsoid[4] + " 2xz + " + ellipsoid[5]
				+ " 2xy + " + ellipsoid[6] + " 2x + " + ellipsoid[7] + " 2y + "
				+ ellipsoid[8] + " 2z + " + ellipsoid[9] + " = 0");
		return ellipsoid;
	}

	/**
	 * Implement the ellipsoid fitting method by Yury Petrov
	 * 
	 * @see <p>
	 *      <a href=
	 *      "http://www.mathworks.com/matlabcentral/fileexchange/24693-ellipsoid-fit"
	 *      >MATLAB script</a>
	 *      </p>
	 * 
	 * @param coOrdinates
	 * @return
	 */
	public static Object[] yuryPetrov(double[][] coOrdinates) {
		// function [ center, radii, evecs, v ] = ellipsoid_fit( X, flag, equals
		// )
		// %
		// % Fit an ellispoid/sphere to a set of xyz data points:
		// %
		// % [center, radii, evecs, pars ] = ellipsoid_fit( X )
		// %
		// % Parameters:
		// % * X, [x y z] - Cartesian data, n x 3 matrix
		// %
		// % Output:
		// % * center - ellispoid center coordinates [xc; yc; zc]
		// % * ax - ellipsoid radii [a; b; c]
		// % * evecs - ellipsoid radii directions as columns of the 3x3 matrix
		// % * v - the 9 parameters describing the ellipsoid algebraically:
		// % Ax^2 + By^2 + Cz^2 + 2Dxy + 2Exz + 2Fyz + 2Gx + 2Hy + 2Iz = 1
		// %
		// % Author:
		// % Yury Petrov, Northeastern University, Boston, MA
		// %
		
		// % need nine or more data points
		// if length( x ) < 9 && flag == 0
		// error( 'Must have at least 9 points to fit a unique ellipsoid' );
		// end
		final int nPoints = coOrdinates.length;
		if (nPoints < 9) {
			System.out
					.print("Too few points; need at least 9 to a unique ellipsoid");
			return null;
		}

		// if flag == 0
		// % fit ellipsoid in the form Ax^2 + By^2 + Cz^2 + 2Dxy + 2Exz + 2Fyz +
		// 2Gx + 2Hy + 2Iz = 1
		// D = [ x .* x, ...
		// y .* y, ...
		// z .* z, ...
		// 2 * x .* y, ...
		// 2 * x .* z, ...
		// 2 * y .* z, ...
		// 2 * x, ...
		// 2 * y, ...
		// 2 * z ];
		double[][] d = new double[nPoints][9];
		for (int i = 0; i < nPoints; i++) {
			final double x = coOrdinates[i][0];
			final double y = coOrdinates[i][1];
			final double z = coOrdinates[i][2];
			d[i][0] = x * x;
			d[i][1] = y * y;
			d[i][2] = z * z;
			d[i][3] = 2 * x * y;
			d[i][4] = 2 * x * z;
			d[i][5] = 2 * y * z;
			d[i][6] = 2 * x;
			d[i][7] = 2 * y;
			d[i][8] = 2 * z;
		}

		
		// % solve the normal system of equations
		// v = ( D' * D ) \ ( D' * ones( size( x, 1 ), 1 ) );
		Matrix D = new Matrix(d);
		Matrix ones = Matrix.ones(nPoints, 1);
		Matrix V = ((D.transpose().times(D)).inverse()).times(D.transpose()
				.times(ones));
		double[] v = V.getColumnPackedCopy();

		//
		// % find the ellipsoid parameters
		// if flag == 0
		// % form the algebraic form of the ellipsoid
		// A = [ v(1) v(4) v(5) v(7); ...
		// v(4) v(2) v(6) v(8); ...
		// v(5) v(6) v(3) v(9); ...
		// v(7) v(8) v(9) -1 ];
		double[][] a = { { v[0], v[3], v[4], v[6] },
				{ v[3], v[1], v[5], v[7] }, { v[4], v[5], v[2], v[8] },
				{ v[6], v[7], v[8], -1 }, };
		Matrix A = new Matrix(a);
		// % find the center of the ellipsoid
		// center = -A( 1:3, 1:3 ) \ [ v(7); v(8); v(9) ];
		Matrix C = (A.getMatrix(0, 2, 0, 2).times(-1).inverse()).times(V
				.getMatrix(6, 8, 0, 0));
		C.printToIJLog("C");
		// % form the corresponding translation matrix
		// T = eye( 4 );
		Matrix T = Matrix.eye(4);
		// T( 4, 1:3 ) = center';
		T.setMatrix(3, 3, 0, 2, C.transpose());
		// % translate to the center
		// R = T * A * T';
		Matrix R = T.times(A.times(T.transpose()));
		// % solve the eigenproblem
		double r33 = R.get(3, 3);
		Matrix R02 = R.getMatrix(0, 2, 0, 2);
		// [ evecs evals ] = eig( R( 1:3, 1:3 ) / -R( 4, 4 ) );
		EigenvalueDecomposition E = new EigenvalueDecomposition(R02.times(-1
				/ r33));
		
		// radii = sqrt( 1 ./ diag( evals ) );
		Matrix eVal = E.getD();
		eVal.printToIJLog("eVal");
		Matrix diagonal = eVal.diag();
		final int nEvals = diagonal.getRowDimension();
		double[] radii = new double[nEvals];
		IJ.log("Radii");
		for (int i = 0; i < nEvals; i++) {
			radii[i] = Math.sqrt(1 / diagonal.get(i, 0));
			IJ.log(""+radii[i]);
		}
		// else
		// v = [ v(1) v(1) v(1) 0 0 0 v(2) v(3) v(4) ];
		// end
		// center = ( -v( 7:9 ) ./ v( 1:3 ) )';
		// gam = 1 + ( v(7)^2 / v(1) + v(8)^2 / v(2) + v(9)^2 / v(3) );
		// radii = ( sqrt( gam ./ v( 1:3 ) ) )';
		// evecs = eye( 3 );
		// end
		double[] centre = C.getColumnPackedCopy();
		double[][] eigenVectors = E.getV().getArrayCopy();
		double[] equation = v;
		Object[] ellipsoid = {centre, radii, eigenVectors, equation};
		return ellipsoid;
	}
}
