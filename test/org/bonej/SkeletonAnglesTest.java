package org.bonej;

import static org.junit.Assert.*;
import ij.ImagePlus;

import org.doube.geometry.TestDataMaker;
import org.junit.Test;

public class SkeletonAnglesTest {

	private final double pi4 = Math.PI / 4;
	private final double pi2 = Math.PI / 2;
	
	private final double[][][] circleCrossResult = {{ 
			{pi4, pi4, pi2},
			null,
			{pi4, pi2, pi4},
			{pi4, pi2, pi4},
			{pi4, pi4, pi2}
		}};
	
	private final double[][][] boxFrameResult = {{ 
		{pi2, pi2, pi2},
		{pi2, pi2, pi2},
		{pi2, pi2, pi2},
		{pi2, pi2, pi2},
		{pi2, pi2, pi2},
		{pi2, pi2, pi2},
		{pi2, pi2, pi2},
		{pi2, pi2, pi2}
	}};

	@Test
	public void testCalculateTriplePointAnglesCrossedCircle() {
		ImagePlus imp = TestDataMaker.crossedCircle(256);
		double[][][] result = (new SkeletonAngles())
				.calculateTriplePointAngles(imp,
						SkeletonAngles.VERTEX_TO_VERTEX);
		for (int g = 0; g < circleCrossResult.length; g++)
			for (int v = 0; v < circleCrossResult[g].length; v++)
				assertArrayEquals(circleCrossResult[g][v], result[g][v], 1e-12);
	}

	@Test
	public void testCalculateTriplePointAnglesBoxFrame() {
		ImagePlus imp = TestDataMaker.boxFrame(128, 128, 128);
		double[][][] result = (new SkeletonAngles())
				.calculateTriplePointAngles(imp,
						SkeletonAngles.VERTEX_TO_VERTEX);
		for (int g = 0; g < boxFrameResult.length; g++)
			for (int v = 0; v < boxFrameResult[g].length; v++)
				assertArrayEquals(boxFrameResult[g][v], result[g][v], 1e-12);
	}
	
	
}
