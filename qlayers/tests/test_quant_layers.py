import numpy as np
import nibabel as nib
import numpy.testing as npt
import os
import pytest
import shutil

from hashlib import sha1
from qlayers import QLayers, load_pickle


class TestQLayers:
    # Generate a 16 x 16 x 16 cube in a 32 x 32 x 32 image
    basic_data = np.zeros((32, 32, 32))
    basic_data[8:24, 8:24, 8:24] = 1
    basic_img = nib.Nifti1Image(basic_data, np.eye(4))

    # Generate a 32 x 32 x 32 cube with values increasing in the y direction
    # to simulate a quantitative map
    basic_map_data, _, _ = np.meshgrid(np.arange(32), np.ones(32), np.ones(32))
    basic_map_img = nib.Nifti1Image(basic_map_data.astype(np.int32), np.eye(4))

    # Generate a 16 x 16 x 16 cube in a 32 x 32 x 32 image where half the
    # cube is 1 and the other half is 2. This is to simulate tissue labels
    basic_tissue_data = np.zeros((32, 32, 32))
    basic_tissue_data[8:16, 8:24, 8:24] = 1
    basic_tissue_data[16:24, 8:24, 8:24] = 2
    basic_tissue_img = nib.Nifti1Image(basic_tissue_data, np.eye(4))

    # A single 2D slice of basic_map_img
    basic_map_data_2d = basic_map_data[:, :, 16]
    basic_map_img_2d = nib.Nifti1Image(
        basic_map_data_2d.astype(np.int32), np.eye(4)
    )

    # Generate a low resolution version of the map above (with 2 mm voxels)
    basic_map_low_res_data, _, _ = np.meshgrid(
        np.arange(16), np.ones(16), np.ones(16)
    )
    basic_map_low_res_img = nib.Nifti1Image(
        basic_map_low_res_data.astype(np.int32), np.eye(4) * 2
    )

    # Generate a 16 x 16 x 16 cube in a 32 x 32 x 32 image with a 4 x 4 x 4
    # cyst in the middle
    basic_data_with_cyst = np.zeros((32, 32, 32))
    basic_data_with_cyst[8:24, 8:24, 8:24] = 1
    basic_data_with_cyst[10:14, 10:14, 10:14] = 0
    basic_img_with_cyst = nib.Nifti1Image(basic_data_with_cyst, np.eye(4))

    # Generate a 256 x 256 x 17 image with two kidneys and a pelvis
    # Kidneys are oblongs rather than the usual kidney shape. Voxel size is
    # 1.5 x 1.5 x 5.5 mm typical of anatomical scans
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

    # Generate a 256 x 256 x 17 image with two kidneys and a small pelvis
    # Kidneys are oblongs rather than the usual kidney shape. Voxel size is
    # 1.5 x 1.5 x 5.5 mm typical of anatomical scans
    kidneys_with_small_pelvis = np.zeros((256, 256, 17))
    kidneys_with_small_pelvis[96:160, 64:96, 5:12] = 1
    kidneys_with_small_pelvis[130:133, 86:96, 6:10] = 0
    kidneys_with_small_pelvis[96:160, 160:192, 5:12] = 1
    kidneys_with_small_pelvis[130:133, 160:170, 6:10] = 0
    aff = np.array(
        [
            [1.5, 0, 0, 0],
            [0, 1.5, 0, 0],
            [0, 0, 5.5, 0],
            [0, 0, 0, 1],
        ]
    )
    kidneys_with_small_pelvis_img = nib.Nifti1Image(kidneys_with_small_pelvis, aff)

    def test_basic_depth(self):
        qlayers = QLayers(self.basic_img)
        depth = qlayers.get_depth()
        npt.assert_almost_equal(depth.max(), 8, decimal=2)
        npt.assert_almost_equal(depth.min(), 0, decimal=2)
        npt.assert_almost_equal(depth.sum(), 9161.554, decimal=2)
        assert np.sum(depth > 0) == 4096

    def test_basic_layers(self):
        expected_layers = np.array(
            [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]
        )
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
        qlayers = QLayers(self.basic_img_with_cyst, fill_ml=1)
        depth = qlayers.get_depth()
        npt.assert_almost_equal(depth.max(), 8, decimal=2)
        npt.assert_almost_equal(depth.min(), 0, decimal=2)
        npt.assert_almost_equal(depth.sum(), 8951.619, decimal=2)
        assert np.sum(depth > 0) == 4032

        # Don't fill the cysts (treat them as their own surface)
        qlayers = QLayers(self.basic_img_with_cyst, fill_ml=0.001)
        depth = qlayers.get_depth()
        npt.assert_almost_equal(depth.max(), 7.222, decimal=2)
        npt.assert_almost_equal(depth.min(), 0, decimal=2)
        npt.assert_almost_equal(depth.sum(), 8574.536, decimal=2)
        assert np.sum(depth > 0) == 4032

    def test_segment_pelvis(self):
        qlayers = QLayers(self.kidneys_with_pelvis_img)
        qlayers._segment_pelvis()
        pelvis = qlayers.pelvis
        assert pelvis.sum() == 640
        assert (
                sha1(pelvis).hexdigest()
                == "524bddb7e7ff85c88a284aaf11a268e26c3c8a73"
        )

    def test_pelvis_dist(self):
        # 20 mm pelvis distance
        qlayers = QLayers(self.kidneys_with_pelvis_img, pelvis_dist=20)
        layers = qlayers.get_layers()
        pelvis = qlayers.pelvis
        assert pelvis.sum() == 640
        assert layers.max() == 21
        assert layers.min() == 0
        assert layers.sum() == 121868
        assert np.sum(layers > 0) == 18984
        assert (
                sha1(layers).hexdigest()
                == "0b23f5f58a19d901902f135c83d388df2290268d"
        )

    def test_iterative_pelvis_dist(self):
        # Pelvis is smaller than the default noise threshold of 2.5 ml
        # Segmetation is repeated with a lower threshold until one is found
        qlayers = QLayers(self.kidneys_with_small_pelvis_img, pelvis_dist=20)
        layers = qlayers.get_layers()
        pelvis = qlayers.pelvis
        assert pelvis.sum() == 136
        assert (
                sha1(pelvis).hexdigest()
                == "57c4b6c815a4817612f1dfde766d0454284581e3"
        )

    def test_save(self):
        qlayers = QLayers(self.kidneys_with_pelvis_img)
        if os.path.exists("test_output"):
            shutil.rmtree("test_output")
        os.makedirs("test_output", exist_ok=True)

        qlayers.save_depth("test_output/depth.nii.gz")
        qlayers.save_layers("test_output/layers.nii.gz")
        qlayers.save_pelvis("test_output/pelvis.nii.gz")
        qlayers.save_surface("test_output/surface.stl")
        qlayers.to_pickle("test_output/qlayers.pkl")
        output_files = os.listdir("test_output")
        assert len(output_files) == 5
        assert "depth.nii.gz" in output_files
        assert "layers.nii.gz" in output_files
        assert "pelvis.nii.gz" in output_files
        assert "surface.stl" in output_files
        assert "qlayers.pkl" in output_files

        for f in os.listdir("test_output"):
            os.remove(os.path.join("test_output", f))
        shutil.rmtree("test_output")

    def test_load_pickle(self):
        qlayers = QLayers(self.kidneys_with_pelvis_img)
        if os.path.exists("test_output"):
            shutil.rmtree("test_output")
        os.makedirs("test_output", exist_ok=True)

        qlayers.to_pickle("test_output/qlayers.pkl")
        qlayers_loaded = load_pickle("test_output/qlayers.pkl")
        assert qlayers_loaded.space == qlayers.space
        assert qlayers_loaded.layers.shape == qlayers.layers.shape
        assert qlayers_loaded.depth.shape == qlayers.depth.shape

        for f in os.listdir("test_output"):
            os.remove(os.path.join("test_output", f))
        shutil.rmtree("test_output")

    def test_not_a_space(self):
        with pytest.raises(ValueError):
            QLayers(self.basic_img, space="invalid")

    def test_add_map_layer_space(self):
        qlayers = QLayers(self.basic_img, space="layers")
        qlayers.add_map(self.basic_map_img, "t1")
        assert len(qlayers.maps) == 1
        assert qlayers.maps[0] == "t1"
        qlayers.add_map(self.basic_map_img, "t2")
        assert len(qlayers.maps) == 2
        assert qlayers.maps[1] == "t2"
        qlayers.add_map(self.basic_map_img, "swi", norm=True)
        assert len(qlayers.maps) == 3
        assert qlayers.maps[2] == "swi"

    def test_add_map_map_space(self):
        qlayers = QLayers(self.basic_img, space="map")
        qlayers.add_map(self.basic_map_img, "t1")
        assert len(qlayers.maps) == 1
        assert qlayers.maps[0] == "t1"
        qlayers.add_map(self.basic_map_img, "t2")
        assert len(qlayers.maps) == 2
        assert qlayers.maps[1] == "t2"
        qlayers.add_map(self.basic_map_img, "swi", norm=True)
        assert len(qlayers.maps) == 3
        assert qlayers.maps[2] == "swi"

    def test_add_map_2d(self):
        # Layers space
        qlayers = QLayers(self.basic_img, space="layers")
        qlayers.add_map(self.basic_map_img_2d, "t1")
        assert len(qlayers.maps) == 1
        assert qlayers.maps[0] == "t1"
        qlayers.add_map(self.basic_map_img_2d, "t2")
        assert len(qlayers.maps) == 2
        assert qlayers.maps[1] == "t2"

        # Map space
        qlayers = QLayers(self.basic_img, space="map")
        qlayers.add_map(self.basic_map_img_2d, "t1")
        assert len(qlayers.maps) == 1
        assert qlayers.maps[0] == "t1"
        qlayers.add_map(self.basic_map_img_2d, "t2")
        assert len(qlayers.maps) == 2
        assert qlayers.maps[1] == "t2"

    def test_remove_maps(self):
        # Without tissues
        qlayers = QLayers(self.basic_img, space="layers")
        qlayers.add_map(self.basic_map_img, "t1")
        qlayers.add_map(self.basic_map_img, "t2")
        qlayers.remove_all_maps()
        assert len(qlayers.maps) == 0

        # With tissues
        qlayers = QLayers(self.basic_img, space="layers")
        qlayers.add_tissue(self.basic_tissue_img)
        qlayers.add_map(self.basic_map_img, "t1")
        qlayers.add_map(self.basic_map_img, "t2")
        qlayers.remove_all_maps()
        assert len(qlayers.maps) == 0

    def test_get_df(self):
        qlayers = QLayers(self.basic_img, space="layers")
        qlayers.add_map(self.basic_map_img, "t1")
        qlayers.add_map(self.basic_map_img, "t2")
        df = qlayers.get_df(format="long")
        assert df.shape == (8192, 4)
        assert df.columns.tolist() == ["depth", "layer", "measurement",
                                       "value"]
        assert df["measurement"].unique().tolist() == ["t1", "t2"]
        npt.assert_array_almost_equal(
            df.mean(numeric_only=True).values,
            np.array([2.2367, 2.5698, 15.500]),
            decimal=2,
        )
        df = qlayers.get_df(format="wide")
        assert df.shape == (4096, 4)
        assert df.columns.tolist() == ["depth", "layer", "t1", "t2"]
        npt.assert_array_almost_equal(
            df.mean(numeric_only=True).values,
            np.array([2.2367, 2.5698, 15.500, 15.500]),
            decimal=2,
        )
        with pytest.raises(ValueError):
            qlayers.get_df(format="invalid")

    def test_map_resampling_layer_space(self):
        qlayers = QLayers(self.basic_img, space="layers")
        qlayers.add_map(self.basic_map_low_res_img, "t1")
        assert len(qlayers.maps) == 1
        assert qlayers.maps[0] == "t1"
        df_long = qlayers.get_df(format="long")
        assert df_long.shape == (4096, 4)
        assert df_long.columns.tolist() == [
            "depth",
            "layer",
            "measurement",
            "value",
        ]
        assert df_long["measurement"].unique().tolist() == ["t1"]
        df_wide = qlayers.get_df(format="wide")
        assert df_wide.shape == (4096, 3)
        assert df_wide.columns.tolist() == ["depth", "layer", "t1"]

        qlayers.add_map(self.basic_map_img, "t2")
        assert len(qlayers.maps) == 2
        assert qlayers.maps == ["t1", "t2"]
        df_long = qlayers.get_df(format="long")
        assert df_long.shape == (8192, 4)
        assert df_long.columns.tolist() == [
            "depth",
            "layer",
            "measurement",
            "value",
        ]
        assert df_long["measurement"].unique().tolist() == ["t1", "t2"]
        df_wide = qlayers.get_df(format="wide")
        assert df_wide.shape == (4096, 4)
        assert df_wide.columns.tolist() == ["depth", "layer", "t1", "t2"]

    def test_map_resampling_map_space(self):
        qlayers = QLayers(self.basic_img, space="map")
        qlayers.add_map(self.basic_map_low_res_img, "t1")
        assert len(qlayers.maps) == 1
        assert qlayers.maps[0] == "t1"
        df_long = qlayers.get_df(format="long")
        assert df_long.shape == (512, 4)
        assert df_long.columns.tolist() == [
            "depth",
            "layer",
            "measurement",
            "value",
        ]
        assert df_long["measurement"].unique().tolist() == ["t1"]

        qlayers.add_map(self.basic_map_img, "t2")
        assert len(qlayers.maps) == 2
        assert qlayers.maps == ["t1", "t2"]
        df_long = qlayers.get_df(format="long")
        assert df_long.shape == (512 + 4096, 4)
        assert df_long.columns.tolist() == [
            "depth",
            "layer",
            "measurement",
            "value",
        ]
        assert df_long["measurement"].unique().tolist() == ["t1", "t2"]
        with pytest.raises(NotImplementedError):
            df_wide = qlayers.get_df(format="wide")

    def test_add_tissue_layer_space(self):
        # Adding the tissue, then the map
        qlayers = QLayers(self.basic_img, space="layers")
        qlayers.add_tissue(self.basic_tissue_img)
        qlayers.add_map(self.basic_map_img, "t1")
        qlayers.add_map(self.basic_map_low_res_img, "t2")
        assert len(qlayers.maps) == 2
        df = qlayers.get_df(format="long")
        assert df.shape == (8192, 5)
        assert df.columns.tolist() == [
            "depth",
            "layer",
            "tissue",
            "measurement",
            "value",
        ]
        assert df["measurement"].unique().tolist() == ["t1", "t2"]
        npt.assert_array_almost_equal(
            df.mean(numeric_only=True).values,
            np.array([2.2367, 2.5698, 1.5, 11.62]),
            decimal=2,
        )
        df = qlayers.get_df(format="wide")
        assert df.shape == (4096, 5)
        assert df.columns.tolist() == ["depth", "layer", "tissue", "t1", "t2"]
        npt.assert_array_almost_equal(
            df.mean(numeric_only=True).values,
            np.array([2.2367, 2.5698, 1.5, 15.500, 7.75]),
            decimal=2,
        )
        with pytest.raises(ValueError):
            qlayers = QLayers(self.basic_img, space="layers")
            qlayers.add_map(self.basic_map_img, "t1")
            qlayers.add_tissue(self.basic_tissue_img)

    def test_add_tissue_map_space(self):
        # Adding the tissue, then the map
        qlayers = QLayers(self.basic_img, space="map")
        qlayers.add_tissue(self.basic_tissue_img)
        qlayers.add_map(self.basic_map_img, "t1")
        qlayers.add_map(self.basic_map_low_res_img, "t2")
        assert len(qlayers.maps) == 2
        df = qlayers.get_df(format="long")
        assert df.shape == (4608, 5)
        assert df.columns.tolist() == [
            "depth",
            "layer",
            "tissue",
            "measurement",
            "value",
        ]
        assert df["measurement"].unique().tolist() == ["t1", "t2"]
        npt.assert_array_almost_equal(
            df.mean(numeric_only=True).values,
            np.array([2.2367, 2.5698, 1.5, 14.61]),
            decimal=2,
        )
        with pytest.raises(ValueError):
            qlayers = QLayers(self.basic_img, space="layers")
            qlayers.add_map(self.basic_map_img, "t1")
            qlayers.add_tissue(self.basic_tissue_img)

    def test_tissue_text_labels(self):
        # Layers space
        qlayers = QLayers(self.basic_img, space="layers")
        qlayers.add_tissue(
            self.basic_tissue_img, tissue_labels=["cortex", "medulla"]
        )
        qlayers.add_map(self.basic_map_img, "t1")
        df_long = qlayers.get_df(format="long")
        assert df_long["tissue"].unique().tolist() == ["cortex", "medulla"]
        df_wide = qlayers.get_df(format="wide")
        assert df_wide["tissue"].unique().tolist() == ["cortex", "medulla"]

        # Map space
        qlayers = QLayers(self.basic_img, space="map")
        qlayers.add_tissue(
            self.basic_tissue_img, tissue_labels=["cortex", "medulla"]
        )
        qlayers.add_map(self.basic_map_img, "t1")
        df_long = qlayers.get_df(format="long")
        assert df_long["tissue"].unique().tolist() == ["cortex", "medulla"]

        # Not enough labels
        with pytest.raises(ValueError):
            qlayers = QLayers(self.basic_img, space="layers")
            qlayers.add_tissue(self.basic_tissue_img, tissue_labels=["cortex"])

        # Too many labels
        with pytest.raises(ValueError):
            qlayers = QLayers(self.basic_img, space="layers")
            qlayers.add_tissue(
                self.basic_tissue_img,
                tissue_labels=["inner cortex", "outer cortex", "medulla"],
            )

        # Labels not a list
        with pytest.raises(ValueError):
            qlayers = QLayers(self.basic_img, space="layers")
            qlayers.add_tissue(self.basic_tissue_img, tissue_labels="cortex")

    def test_norm_map(self):
        qlayers = QLayers(self.basic_img, space="map")
        map_norm = qlayers._normalise_data(self.basic_map_data)
        npt.assert_almost_equal(map_norm.mean(), 0, decimal=5)
        npt.assert_almost_equal(map_norm.std(), 1, decimal=5)
