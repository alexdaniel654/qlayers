import numpy as np
import nibabel as nib
import pytest
import numpy.testing as npt

from qlayers import QLayers


class TestQLayers:
    basic_data = np.zeros((32, 32, 32))
    basic_data[8:24, 8:24, 8:24] = 1
    basic_img = nib.Nifti1Image(basic_data, np.eye(4))

    basic_data_with_cyst = np.zeros((32, 32, 32))
    basic_data_with_cyst[8:24, 8:24, 8:24] = 1
    basic_data_with_cyst[10:14, 10:14, 10:14] = 0
    basic_img_with_cyst = nib.Nifti1Image(basic_data_with_cyst, np.eye(4))

    kidneys_with_pelvis = np.zeros((256, 256, 17))
    kidneys_with_pelvis[96:160, 64:96, 5:12] = 1
    kidneys_with_pelvis[123:133, 86:96, 6:10] = 0
    kidneys_with_pelvis[96:160, 160:192, 5:12] = 1
    kidneys_with_pelvis[123:133, 160:170, 6:10] = 0
    aff = np.array(
        [
            [1.5, 0, 0, 0],
            [0, 1.5, 0, 0],
            [0, 0, 5.5, 0],
            [0, 0, 0, 1],
        ]
    )
    kidneys_with_pelvis_img = nib.Nifti1Image(kidneys_with_pelvis, aff)

    def test_basic_depth(self):
        qlayers = QLayers(self.basic_img)
        depth = qlayers.get_depth()
        npt.assert_almost_equal(depth.max(), 8, decimal=2)
        npt.assert_almost_equal(depth.min(), 0, decimal=2)
        npt.assert_almost_equal(depth.sum(), 9161.554, decimal=2)
        assert np.sum(depth > 0) == 4096

    def test_basic_layers(self):
        expected_layers = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0])
        qlayers = QLayers(self.basic_img, thickness=1)
        layers = qlayers.get_layers()
        assert layers.max() == 9
        assert layers.min() == 0
        assert layers.sum() == 10526
        assert np.sum(layers > 0) == 4096
        npt.assert_array_equal(np.unique(layers), expected_layers)

    def test_thickness(self):
        expected_layers = np.array([0.0, 2.0, 4.0, 6.0, 8.0, 10.0])
        qlayers = QLayers(self.basic_img, thickness=2)
        layers = qlayers.get_layers()
        assert layers.max() == 10
        assert layers.min() == 0
        assert layers.sum() == 12960
        assert np.sum(layers > 0) == 4096
        npt.assert_array_equal(np.unique(layers), expected_layers)

    def test_fill_cysts(self):
        qlayers = QLayers(self.basic_img_with_cyst)
        depth = qlayers.get_depth()
        npt.assert_almost_equal(depth.max(), 8, decimal=2)
        npt.assert_almost_equal(depth.min(), 0, decimal=2)
        npt.assert_almost_equal(depth.sum(), 8951.619, decimal=2)
        assert np.sum(depth > 0) == 4032

    def test_fill_ml(self):
        # Fill the cysts
        qlayers = QLayers(self.basic_img_with_cyst, fill_ml=0.065)
        depth = qlayers.get_depth()
        npt.assert_almost_equal(depth.max(), 8, decimal=2)
        npt.assert_almost_equal(depth.min(), 0, decimal=2)
        npt.assert_almost_equal(depth.sum(), 8951.619, decimal=2)
        assert np.sum(depth > 0) == 4032

        # Don't fill the cysts (treat them as their own surface)
        qlayers = QLayers(self.basic_img_with_cyst, fill_ml=0.063)
        depth = qlayers.get_depth()
        npt.assert_almost_equal(depth.max(), 7.222, decimal=2)
        npt.assert_almost_equal(depth.min(), 0, decimal=2)
        npt.assert_almost_equal(depth.sum(), 8574.536, decimal=2)
        assert np.sum(depth > 0) == 4032

    # TODO pelvis segmentation tests
    def test_segment_pelvis(self):
        qlayers = QLayers(self.kidneys_with_pelvis_img)
        qlayers._segment_pelvis()
        pelvis = qlayers.pelvis
        assert pelvis.sum() == 640
        # assert hash(pelvis.tostring()) == 8678233953073131695

    # TODO get_df tests
    # TODO save tests
    # TODO add_map tests
    # TODO remove_all_maps test
