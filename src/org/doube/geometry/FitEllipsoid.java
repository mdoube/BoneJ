package org.doube.geometry;

import ij.IJ;

import org.doube.jama.EigenvalueDecomposition;
import org.doube.jama.Matrix;

public class FitEllipsoid {
	public int maxPoints = 1000;

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
	public double[] fitLiGriffiths(double[][] coOrdinates) {
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
	public Object[] yuryPetrov(double[][] coOrdinates) {
		// function [ center, radii, evecs, v ] = ellipsoid_fit( X, flag, equals
		// )
		// %
		// % Fit an ellispoid/sphere to a set of xyz data points:
		// %
		// % [center, radii, evecs, pars ] = ellipsoid_fit( X )
		// % [center, radii, evecs, pars ] = ellipsoid_fit( [x y z] );
		// % [center, radii, evecs, pars ] = ellipsoid_fit( X, 1 );
		// % [center, radii, evecs, pars ] = ellipsoid_fit( X, 2, 'xz' );
		// % [center, radii, evecs, pars ] = ellipsoid_fit( X, 3 );
		// %
		// % Parameters:
		// % * X, [x y z] - Cartesian data, n x 3 matrix or three n x 1 vectors
		// % * flag - 0 fits an arbitrary ellipsoid (default),
		// % - 1 fits an ellipsoid with its axes along [x y z] axes
		// % - 2 followed by, say, 'xy' fits as 1 but also x_rad = y_rad
		// % - 3 fits a sphere
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
		//
		// else
		// x = X( :, 1 );
		// y = X( :, 2 );
		// z = X( :, 3 );
		// end
		//
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
		// if length( x ) < 6 && flag == 1
		// error( 'Must have at least 6 points to fit a unique oriented
		// ellipsoid' );
		// end
		// if length( x ) < 5 && flag == 2
		// error( 'Must have at least 5 points to fit a unique oriented
		// ellipsoid with two axes equal' );
		// end
		// if length( x ) < 3 && flag == 3
		// error( 'Must have at least 4 points to fit a unique sphere' );
		// end
		//
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
		// 2 * z ]; % ndatapoints x 9 ellipsoid parameters
		// elseif flag == 1
		// % fit ellipsoid in the form Ax^2 + By^2 + Cz^2 + 2Gx + 2Hy + 2Iz = 1
		// D = [ x .* x, ...
		// y .* y, ...
		// z .* z, ...
		// 2 * x, ...
		// 2 * y, ...
		// 2 * z ]; % ndatapoints x 6 ellipsoid parameters
		// elseif flag == 2
		// % fit ellipsoid in the form Ax^2 + By^2 + Cz^2 + 2Gx + 2Hy + 2Iz = 1,
		// % where A = B or B = C or A = C
		// if strcmp( equals, 'yz' ) || strcmp( equals, 'zy' )
		// D = [ y .* y + z .* z, ...
		// x .* x, ...
		// 2 * x, ...
		// 2 * y, ...
		// 2 * z ];
		// elseif strcmp( equals, 'xz' ) || strcmp( equals, 'zx' )
		// D = [ x .* x + z .* z, ...
		// y .* y, ...
		// 2 * x, ...
		// 2 * y, ...
		// 2 * z ];
		// else
		// D = [ x .* x + y .* y, ...
		// z .* z, ...
		// 2 * x, ...
		// 2 * y, ...
		// 2 * z ];
		// end
		// else
		// % fit sphere in the form A(x^2 + y^2 + z^2) + 2Gx + 2Hy + 2Iz = 1
		// D = [ x .* x + y .* y + z .* z, ...
		// 2 * x, ...
		// 2 * y, ...
		// 2 * z ]; % ndatapoints x 4 sphere parameters
		// end
		//
		// % solve the normal system of equations
		// v = ( D' * D ) \ ( D' * ones( size( x, 1 ), 1 ) );
		//
		// % find the ellipsoid parameters
		// if flag == 0
		// % form the algebraic form of the ellipsoid
		// A = [ v(1) v(4) v(5) v(7); ...
		// v(4) v(2) v(6) v(8); ...
		// v(5) v(6) v(3) v(9); ...
		// v(7) v(8) v(9) -1 ];
		// % find the center of the ellipsoid
		// center = -A( 1:3, 1:3 ) \ [ v(7); v(8); v(9) ];
		// % form the corresponding translation matrix
		// T = eye( 4 );
		// T( 4, 1:3 ) = center';
		// % translate to the center
		// R = T * A * T';
		// % solve the eigenproblem
		// [ evecs evals ] = eig( R( 1:3, 1:3 ) / -R( 4, 4 ) );
		// radii = sqrt( 1 ./ diag( evals ) );
		// else
		// if flag == 1
		// v = [ v(1) v(2) v(3) 0 0 0 v(4) v(5) v(6) ];
		// elseif flag == 2
		// if strcmp( equals, 'xz' ) || strcmp( equals, 'zx' )
		// v = [ v(1) v(2) v(1) 0 0 0 v(3) v(4) v(5) ];
		// elseif strcmp( equals, 'yz' ) || strcmp( equals, 'zy' )
		// v = [ v(2) v(1) v(1) 0 0 0 v(3) v(4) v(5) ];
		// else % xy
		// v = [ v(1) v(1) v(2) 0 0 0 v(3) v(4) v(5) ];
		// end
		// else
		// v = [ v(1) v(1) v(1) 0 0 0 v(2) v(3) v(4) ];
		// end
		// center = ( -v( 7:9 ) ./ v( 1:3 ) )';
		// gam = 1 + ( v(7)^2 / v(1) + v(8)^2 / v(2) + v(9)^2 / v(3) );
		// radii = ( sqrt( gam ./ v( 1:3 ) ) )';
		// evecs = eye( 3 );
		// end

		Object[] ellipsoid = {};
		return ellipsoid;
	}
}
