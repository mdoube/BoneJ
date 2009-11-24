package org.doube.bonej;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.Plot;
import ij.plugin.PlugIn;

import org.doube.bonej.FastFourierTransform;

/**
 * Calculate and analyse the 3D autocorrelation function
 * 
 * @author Michael Doube
 * 
 */
public class AutoCorrelation implements PlugIn {

	public void run(String arg) {
		if (!ImageCheck.checkIJVersion())
			return;
		ImagePlus imp = IJ.getImage();

		ImagePlus acf = getACF(imp);

		Object result = getAnisotropy(acf);

		return;
	}

	/**
	 * Calculate anisotropy from an autocorrelation function. Also return the
	 * width of the peaks, for bone this relates to Tb.Th and Tb.Sp
	 * 
	 * @param acf
	 *            ImagePlus containing the autocorrelation function
	 * @return results in an object array
	 */
	private Object getAnisotropy(ImagePlus acf) {
		/*
		 * >> > > Now, I can only _guess_ that the result (basically, the >> > >
		 * autocorrelation with respect to all possible offsets) is what they >>
		 * > > refer to as the "anisotropy tensor", and I would expect that the
		 * >> > > analysis boils down to determining the PCA of the offset
		 * vectors >> > > (weighted by the corresponding autocorrelation value).
		 * > > > > I get lost at that part; I think I will have to write to the
		 * authors.
		 * 
		 * The result of the FFTs followed by the point-wise product and the
		 * reverse FFT would be a stack consisting of the cross-correlation
		 * values for every offset vector. I.e. for all dx,dy,dz, you would have
		 * 
		 * CC(dx,dy,dz) = \sum A(x,y,z)A(x+dx,y+dy,z+dz)
		 * 
		 * Now, I think that this is what they refer to as the
		 * "anisotropy tensor".
		 * 
		 * By taking the principal components of this field, i.e the most
		 * prominent directions of cross-correlation, you should be able to come
		 * up with a measure for anisotropy.
		 * 
		 * But this is a little hand-waving, I am sorry...
		 * 
		 * Ciao, Dscho
		 */
		// TODO Auto-generated method stub
		return null;
	}

	/**
	 * Calculate the autocorrelation function from an input image
	 * 
	 * @param imp
	 * @return
	 */
	private ImagePlus getACF(ImagePlus imp) {
		FastFourierTransform fft = new FastFourierTransform();
		fft.run("");
		/*
		 * the Fourier transform of the complex conjugate of a function f(x) is
		 * F*(-s), the reflection of the conjugate of the transform
		 * 
		 * The autocorrelation of an image with itself can be described as a >>
		 * > > convolution of the image with the flipped version of itself
		 * (where >> > > "flipped" means "flipped in all available axis"). > > >
		 * convolution can be >> > > described as a point-wise multiplication in
		 * Fourier space, 
		 * 
		 * Actually, you need to do point-wise multiplications with the FFT of
		 * the reflected version (this is not the same as the reflected FFT).
		 * 
		 * AFAIR there is a very simple transformation of the FFT that turns it
		 * into the FFT of the reflected image. Basically, it boils down to the
		 * equations
		 * 
		 * cos(-x) = cos(x) sin(-x) = -sin(x)
		 * 
		 * So, basically you have to negate the imaginary coefficients of all
		 * odd, and the real coefficients of all even factors.
		 * 
		 * It should be easily verified by taking the FFT of an image and then
		 * an FFT of the flipped (AKA reflected) image, and then having a look
		 * at the actual values.
		 * 
		 * >> > > You have to extend the dimensions, though, as Fourier wraps
		 * around at >> > > the borders. For example, you have to embed a cube
		 * into a cube that >> > > is 2x2x2 times as large, filling the rest
		 * with zeroes. > > > > Extend the dimensions of the cube of original
		 * data or its Fourier > > transform?
		 * 
		 * The original data.
		 * 
		 * Remember, your correlation would look something like this (in 2D, for
		 * simplicity):
		 * 
		 * \sum A(x,y)A(x+dx,y+dy)
		 * 
		 * in the range where both x,y and x+dx,y+dy are inside the image.
		 * Extending the original data by doubling all Cartesian axes naturally
		 * gives the same result, as the additional summands will be 0 due to
		 * one of A(x,y) and A(x+dx,y+dy) being outside of the original image.
		 * 
		 * If you do not do that, the wrap-around nature of FFT will lead to
		 * summands that would look something like this:
		 * 
		 * A(x,y)A(x+dx-w,y+dy-h)
		 * 
		 * which is obviously not what you want.
		 */
		// TODO Auto-generated method stub
		return null;
	}

	/**
	 * 1D autocorrelation
	 * 
	 * @param x
	 * @return
	 */
	private float[] autoCorrelation(float[] x) {
		int size = x.length;
		float[] R = new float[size];
		float sum;

		for (int i = 0; i < size; i++) {
			sum = 0;
			for (int j = 0; j < size - i; j++) {
				sum += x[j] * x[j + i];
			}
			R[i] = sum;
		}
		return R;
	}

}
