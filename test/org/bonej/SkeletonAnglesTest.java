package org.bonej;

import static org.junit.Assert.*;
import ij.ImagePlus;

import org.doube.geometry.TestDataMaker;
import org.junit.Test;

public class SkeletonAnglesTest {

	private final double pi4 = Math.PI / 4;
	private final double pi2 = Math.PI / 2;
	/** Not quite pi / 2 because vertex isn't exactly at corner */
	private final double pi2ish = 1.58555518815571;
	
	private final double[][][] circleCrossResult = {{ 
			{pi4, pi4, pi2},
			null,
			{pi4, pi2, pi4},
			{pi4, pi2, pi4},
			{pi4, pi4, pi2}
		}};
	
	private final double[][][] circleCrossResultNth = {{ 
		{pi2, pi2, Math.PI},
		null,
		{pi2, Math.PI, pi2},
		{pi2, Math.PI, pi2},
		{pi2, pi2, Math.PI}
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
	
	private final double[][][] boxFrameResultNth = {{ 
		{pi2ish, pi2ish, pi2ish},
		{pi2ish, pi2ish, pi2ish},
		{pi2ish, pi2ish, pi2ish},
		{pi2ish, pi2ish, pi2ish},
		{pi2ish, pi2ish, pi2ish},
		{pi2ish, pi2ish, pi2ish},
		{pi2ish, pi2ish, pi2ish},
		{pi2ish, pi2ish, pi2ish}
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
	public void testCalculateTriplePointAnglesCrossedCircleNth() {
		ImagePlus imp = TestDataMaker.crossedCircle(256);
		double[][][] result = (new SkeletonAngles())
				.calculateTriplePointAngles(imp, 8);
		for (int g = 0; g < circleCrossResultNth.length; g++)
			for (int v = 0; v < circleCrossResultNth[g].length; v++)
				assertArrayEquals(circleCrossResultNth[g][v], result[g][v], 1e-12);
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
	
	@Test
	public void testCalculateTriplePointAnglesBoxFrameNth() {
		ImagePlus imp = TestDataMaker.boxFrame(128, 128, 128);
		double[][][] result = (new SkeletonAngles())
				.calculateTriplePointAngles(imp, 32);
		for (int g = 0; g < boxFrameResultNth.length; g++)
			for (int v = 0; v < boxFrameResultNth[g].length; v++)
				assertArrayEquals(boxFrameResultNth[g][v], result[g][v], 1e-12);
	}
}
