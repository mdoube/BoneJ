package org.doube.util;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.when;

import org.junit.Test;

import ij.ImagePlus;
import ij.measure.Calibration;

/**
 * Unit tests for the org.doube.util.ImageCheck class
 *
 * Richard Domander
 */
public class ImageCheckTest {
    @Test
    public void testIsVoxelIsotropicReturnsFalseIfImageIsNull() throws Exception {
        boolean result = ImageCheck.isVoxelIsotropic(null);
        assertFalse("Null image should not be isotropic", result);
    }

    @Test
    public void testIsVoxelIsotropic() throws Exception {
        ImagePlus testImage = mock(ImagePlus.class);
        Calibration anisotropicCalibration = new Calibration();

        // 2D anisotropic image with 0 tolerance
        anisotropicCalibration.pixelWidth = 2;
        anisotropicCalibration.pixelHeight = 1;

        when(testImage.getCalibration()).thenReturn(anisotropicCalibration);
        when(testImage.getStackSize()).thenReturn(1);

        boolean result = ImageCheck.isVoxelIsotropic(testImage, 0.0);
        assertFalse("Image where width > height should not be isotropic", result);

        // 2D image where anisotropy is within tolerance
        result = ImageCheck.isVoxelIsotropic(testImage, 1.0);
        assertTrue("Image should be isotropic if anisotropy is within tolerance", result);

        // 3D image where depth anisotropy is beyond tolerance
        anisotropicCalibration.pixelDepth = 1000;
        when(testImage.getStackSize()).thenReturn(100);

        result = ImageCheck.isVoxelIsotropic(testImage, 1.0);
        assertFalse("Pixel depth too great to be anisotropic within tolerance", result);
    }
}