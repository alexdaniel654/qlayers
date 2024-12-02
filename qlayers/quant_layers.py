import nibabel as nib
import numpy as np
import pandas as pd
import pickle
import trimesh
import warnings
from nibabel.processing import resample_from_to
from skimage.measure import marching_cubes
from skimage.morphology import (opening, remove_small_holes,
                                remove_small_objects)
from tqdm import tqdm
from trimesh import smoothing

from .utils import convex_hull_objects, pad_dimensions


class QLayers:
    """
    This class takes in a mask of a kidney and calculates the depth of each
    voxel from the surface of the kidney. The kidney is then divided into
    layers of a specified thickness that can be used as regions of interest.
    It can then be used to add maps of other quantitative parameters to the
    kidney and generate dataframes of these parameters with depth/layer.

    Parameters
    ----------
    mask_img : nibabel.nifti1.Nifti1Image
        Binary mask of kidney
    thickness : float, optional
        Default 1
        Thickness of layers to use when quantising depth, in millimeters
    fill_ml : float, optional
        Default 10
        Volume of holes in the mask to fill (usually cycst), in cubic
        centimeters (ml)
    pelvis_dist : float, optional
        Default 0
        Voxels within `pelvis_dist` of the renal pelvis are excluded from
        the resulting depth/layer calculations.
        If 0, no pelvis segmentation is performed
    space : {"map", "layers"}, optional
        Default "map"
        If "map", the depth map is resampled to the space of the
        depth/layers are resampled to the space of the quantitative map. If
        "layers", the quantitative map is resampled to the space of the
        layers/depth. "map" gives more accurate quantitative results as the
        map is not resampled however "layers" allows for the output of a
        wide dataframe where each row contains the depth, layer and all
        quantitative measurements for a voxel.
    """

    # TODO add verbose option

    def __init__(
            self, mask_img, thickness=1, fill_ml=10, pelvis_dist=0, space="map"
    ):
        self.mask_img = mask_img
        self.mask = mask_img.get_fdata() > 0.5
        self.zoom = mask_img.header.get_zooms()
        self.affine = mask_img.affine
        self.thickness = thickness
        self.fill_ml = fill_ml
        self.pelvis_dist = pelvis_dist
        if space not in ["map", "layers"]:
            raise ValueError("space must be 'map' or 'layers'")
        self.space = space
        self.depth = self._calculate_depth()
        self.layers = np.ceil(self.depth * (1 / self.thickness)) / (
                1 / self.thickness
        )
        self.layers_list = np.unique(self.layers)
        self.df_long = pd.DataFrame()
        self.maps = []
        self._tissue_data_ls = None
        self._tissue_labels = None
        self._tissue_img = None
        if self.space == "layers":
            self.df_wide = pd.DataFrame(columns=["depth", "layer"])
            self.df_wide["depth"] = self.depth[self.mask]
            self.df_wide["layer"] = self.layers[self.mask]

    def add_map(self, map_img, name, norm=False):
        """
        Add a quantitative map to the object. The either map or
        layers will be resampled to a common space depending on the
        `space` parameter of the `QLayers` object.

        Parameters
        ----------
        map_img : nibabel.nifti1.Nifti1Image
            Quantitative map to be added to the object
        name : str
            Name of the quantitative map to be added. This name will be used
            as a column name in the output dataframe.
        norm: bool, optional
            Default False
            If True, the map will be normalized to have a mean of 0 and a
            standard deviation of 1 before being added to the object.
        """
        self.maps.append(name)
        if map_img.ndim == 2:
            map_img = pad_dimensions(map_img)

        if self.space == "layers":
            # Resample map into space of layers. Doing cval as a big and
            # unusual number as cval=np.nan doesn't work
            map_img = resample_from_to(map_img, self.mask_img,
                                       cval=2 ** 16 - 2)
            map_data = map_img.get_fdata()
            map_data[map_data == 2 ** 16 - 2] = np.nan

            if norm:
                map_data = self._normalise_data(map_data, self.mask)

            self.df_wide[name] = map_data[self.mask]

            if self._tissue_img is None:
                sub_df = pd.DataFrame(
                    columns=["depth", "layer", "measurement", "value"]
                )
            else:
                sub_df = pd.DataFrame(
                    columns=[
                        "depth",
                        "layer",
                        "tissue",
                        "measurement",
                        "value",
                    ]
                )
                sub_df["tissue"] = self._tissue_data_ls[self.mask]
            sub_df["depth"] = self.depth[self.mask]
            sub_df["layer"] = self.layers[self.mask]
            sub_df["measurement"] = name
            sub_df["value"] = map_data[self.mask]
            if self.df_long.empty:
                self.df_long = sub_df
            else:
                self.df_long = pd.concat([self.df_long, sub_df])

        if self.space == "map":
            # Resample layers into space of map
            layers_img = nib.Nifti1Image(self.layers, self.affine)
            layers_img_rs = resample_from_to(layers_img, map_img, order=0)
            layers_rs = layers_img_rs.get_fdata()

            depth_img = nib.Nifti1Image(self.depth, self.affine)
            depth_img_rs = resample_from_to(depth_img, map_img)
            depth_rs = depth_img_rs.get_fdata()

            mask_img_rs = resample_from_to(self.mask_img, map_img, order=0)
            mask_rs = mask_img_rs.get_fdata() > 0.5
            map_data = map_img.get_fdata()

            if norm:
                map_data = self._normalise_data(map_data, mask_rs)

            if self._tissue_img is None:
                sub_df = pd.DataFrame(
                    columns=["depth", "layer", "measurement", "value"]
                )
            else:
                sub_df = pd.DataFrame(
                    columns=["depth", "layer", "tissue", "measurement",
                             "value"]
                )
                tissue_rs = resample_from_to(
                    self._tissue_img, map_img, order=0
                ).get_fdata()
                sub_df["tissue"] = tissue_rs[mask_rs]
            sub_df["depth"] = depth_rs[mask_rs]
            sub_df["layer"] = layers_rs[mask_rs]
            sub_df["measurement"] = name
            sub_df["value"] = map_data[mask_rs]
            self.df_long = pd.concat([self.df_long, sub_df])

    def add_tissue(self, tissue_img, tissue_labels=None):
        """
        Add a tissue segmentation image to the object. This segmentation
        should contain integer labels for each tissue type in the kidney
        where 0 is background. The labels will be used in the output
        dataframe. The tissue segmentation will be resampled to the space of
        the layers if the `space` parameter of the `QLayers` object is set
        to "layers" or the space of the quantitative maps added if the
        `space` parameter is set to "map".

        Parameters
        ----------
        tissue_img : nibabel.nifti1.Nifti1Image
            Tissue segmentation image to be added to the object
        tissue_labels : list of str, optional
            Default None
            List of text labels for each tissue type in the segmentation
            e.g. ["cortex", "medulla"] if 1 in `tissue_img` represents
            renal cortex and 2 represents medulla. If None, the integer labels
            from `tissue_img` will be used.
        """
        if len(self.maps) > 0:
            raise ValueError(
                "Tissue labels must be added before any quantitative maps."
            )

        self._tissue_img = tissue_img
        self._tissue_labels = tissue_labels
        self._tissue_data_ls = resample_from_to(
            tissue_img, self.mask_img, order=0, cval=2 ** 16 - 2
        ).get_fdata()
        self._tissue_data_ls[self._tissue_data_ls == 2 ** 16 - 2] = np.nan
        self._tissue_data_ls[self._tissue_data_ls == 0] = np.nan

        if self._tissue_labels is not None:
            if len(np.unique(np.nan_to_num(self._tissue_data_ls))) - 1 != len(
                    self._tissue_labels
            ):
                raise ValueError(
                    "Number of tissue labels must equal number of unique "
                    "non-zero values in tissue image."
                )
        if self.space == "layers":
            self.df_wide["tissue"] = self._tissue_data_ls[self.mask]

    def get_df(self, format="long"):
        """
        Returns a dataframe of all the quantitative maps added to the object
        with the depth/layer of each voxel.

        Parameters
        ----------
        format : {"wide", "long"}, optional
            Default "long"
            If "long", the dataframe is returned in a long
            format where each row contains the depth, layer and a single
            quantitative measurement for a voxel. If "wide", the dataframe is
            returned in a wide format where each row contains the depth,
            layer and all quantitative measurements for a voxel. This option
            is only available if the `space` parameter of the `QLayers`
            object is set to "layers".

        Returns
        -------
        pandas.DataFrame
            Dataframe of all quantitative maps added to the object with the
            depth/layer of each voxel.
        """
        if self._tissue_labels is not None:
            if len(self._tissue_labels) > 1:
                _tissue_labels_type = type(self._tissue_labels[0])
            else:
                _tissue_labels_type = type(self._tissue_labels)

            if len(self.df_long) > 0:
                self.df_long = self.df_long.dropna(subset=["tissue"])
                self.df_long["tissue"] = self.df_long["tissue"].astype(
                    _tissue_labels_type
                )
                for ind, label in zip(
                        np.sort(self.df_long["tissue"].unique()),
                        self._tissue_labels
                ):
                    self.df_long.loc[
                        self.df_long["tissue"] == ind, "tissue"
                    ] = label
            if self.space == "layers":
                self.df_wide = self.df_wide.dropna(subset=["tissue"])
                self.df_wide["tissue"] = self.df_wide["tissue"].astype(
                    _tissue_labels_type
                )
                for ind, label in zip(
                        np.sort(self.df_wide["tissue"].unique()),
                        self._tissue_labels,
                ):
                    self.df_wide.loc[
                        self.df_wide["tissue"] == ind, "tissue"
                    ] = label

        if format == "wide":
            if self.space == "map":
                raise NotImplementedError(
                    "Data cannot be retrieved in wide "
                    "format when space is 'map'. This "
                    "is because each quantitative map "
                    "is in a different space so each "
                    "row of a wide dataframe would not "
                    "come from the same space in the "
                    "kidney."
                )
            else:
                return self.df_wide.loc[self.df_wide["depth"] > 0]
        elif format == "long":
            if "tissue" in self.df_long.columns:
                self.df_long = self.df_long[
                    ["depth", "layer", "tissue", "measurement", "value"]
                ]
            return self.df_long.loc[self.df_long["depth"] > 0]
        else:
            raise ValueError("format must be 'wide' or 'long'")

    def get_depth(self):
        """
        Returns distance from each voxel to the surface of the kidney as a
        numpy array.

        Returns
        -------
        numpy.ndarray
            Distance from each voxel to the surface of the kidney
        """
        return self.depth

    def get_layers(self):
        """
        Returns layer of each voxel in the kidney as a numpy array.

        Returns
        -------
        numpy.ndarray
            Layer of each voxel
        """
        return self.layers

    def remove_all_maps(self):
        """
        Removes all quantitative maps from the object and resets the
        dataframe to only contain depth/layer information.
        """
        self.maps = []
        if self._tissue_img is None:
            self.df_long = pd.DataFrame(
                columns=["depth", "layer", "measurement", "value"]
            )
        else:
            self.df_long = pd.DataFrame(
                columns=["depth", "layer", "tissue", "measurement", "value"]
            )

        if self.space == "layers":
            self.df_wide = pd.DataFrame(columns=["depth", "layer"])
            self.df_wide["depth"] = self.depth[self.mask]
            self.df_wide["layer"] = self.layers[self.mask]
            if self._tissue_img is not None:
                self.df_wide["tissue"] = self._tissue_data_ls[self.mask]

    def save_depth(self, fname):
        """
        Saves the depth map as a nifti file.

        Parameters
        ----------
        fname : str
            Filename to save the depth map to.
        """
        depth_img = nib.Nifti1Image(self.depth, self.affine)
        nib.save(depth_img, fname)

    def save_layers(self, fname):
        """
        Saves the layers as a nifti file.

        Parameters
        ----------
        fname : str
            Filename to save the layers to.
        """
        layer_img = nib.Nifti1Image(self.layers, self.affine)
        nib.save(layer_img, fname)

    def save_pelvis(self, fname):
        """
        Saves the pelvis segmentation as a nifti file.

        Parameters
        ----------
        fname : str
            Filename to save the pelvis segmentation to.
        """
        if not hasattr(self, "pelvis"):
            self._segment_pelvis()
        pelvis_img = nib.Nifti1Image(self.pelvis.astype(np.int32), self.affine)
        nib.save(pelvis_img, fname)

    def save_surface(self, fname):
        """
        Saves the kidney surface as a stl file.

        Parameters
        ----------
        fname : str
            Filename to save the kidney surface to.
        """
        if not hasattr(self, "smooth_mesh"):
            self._calculate_depth()
        self.smooth_mesh.export(fname)

    def to_pickle(self, fname):
        """
        Saves the QLayers object as a pickle file.

        Warning: The ability to import pickle files between Python versions
        is not guaranteed, as such this is not recommended for long term
        storage of QLayers objects!

        Parameters
        ----------
        fname : str
            Filename to save the QLayers to.
        """
        with open(fname, "wb") as f:
            pickle.dump(self, f)

    def _calculate_depth(self):
        """
        Calculates the distance from each voxel to the surface of the kidney.

        Returns
        -------
        numpy.ndarray
            Distance from each voxel to the surface of the kidney
        """
        # Fill any holes in the mask with volume less than fill_ml
        # (measured in millileters)
        fill_vox = int(self.fill_ml / (np.prod(self.zoom) / 1000))
        mask_filled = remove_small_holes(self.mask, fill_vox)

        # Convert the voxel mask into a mesh using the marching cubes
        # algorithm and trimesh
        print("Making Mesh")
        verts, faces, normals, _ = marching_cubes(
            mask_filled.astype(np.uint8),
            spacing=self.zoom,
            level=0.5,
            step_size=1.0,
        )
        mesh = trimesh.Trimesh(
            vertices=verts, faces=faces, vertex_normals=normals
        )

        # Smooth the resulting mesh
        print("Smoothing Mesh")
        mesh = smoothing.filter_mut_dif_laplacian(mesh, lamb=1, iterations=50)
        self.smooth_mesh = mesh

        # Generate a pointcloud of query points
        x, y, z = np.meshgrid(
            (np.arange(self.mask.shape[0]) * self.zoom[0]),
            (np.arange(self.mask.shape[1]) * self.zoom[1]),
            (np.arange(self.mask.shape[2]) * self.zoom[2]),
            indexing="ij",
        )
        x, y, z = x.reshape(-1), y.reshape(-1), z.reshape(-1)
        points = np.array([x, y, z]).T

        # Find the nearest surface to each point inside the kidney
        points = points[self.mask.reshape(-1) > 0.5]
        batch_size = 10000
        batches = np.ceil(len(points) / batch_size).astype(int)
        distances = np.zeros(len(points))
        i = 0
        # print("Calculating Distances")
        with tqdm(total=batches, desc="Distance Calculation") as pbar:
            while i < len(points):
                (closest_points, dist, triangle_id) = mesh.nearest.on_surface(
                    points[i:i + batch_size]
                )
                distances[i:i + batch_size] = dist
                i += batch_size
                pbar.update()

        # Write these distances to voxels in the shape of the original image
        depth = np.zeros(self.mask.shape)
        depth[self.mask > 0.5] = distances
        if self.pelvis_dist != 0:
            noise_ml = 2.5  # Start with pretty high noise exclusion
            self._segment_pelvis(noise_ml=2.5)

            # If no pelvis is found, reduce the noise exclusion until a pelvis
            # is found or noise_ml is 0
            while (np.sum(self.pelvis) == 0) and (noise_ml > 0):
                noise_ml -= 0.5
                warnings.warn(f"No pelvis found, reducing noise exclusion to {noise_ml} ml")
                self._segment_pelvis(noise_ml=noise_ml)

            if np.sum(self.pelvis) > 0:
                verts_p, faces_p, normals_p, _ = marching_cubes(
                    self.pelvis.astype(np.uint8),
                    spacing=self.zoom,
                    level=0.5,
                    step_size=1.0,
                )

                mesh_p = trimesh.Trimesh(
                    vertices=verts_p, faces=faces_p, vertex_normals=normals_p
                )

                i = 0
                distances_p = np.zeros(len(points))
                with tqdm(total=batches,
                          desc="Pelvis Distance Calculation") as pbar:
                    while i < len(points):
                        (
                            closest_points_p,
                            dist_p,
                            triangle_id_p,
                        ) = mesh_p.nearest.on_surface(points[i:i + batch_size])

                        distances_p[i:i + batch_size] = dist_p
                        i += batch_size
                        pbar.update()
                depth_p = np.zeros(self.mask.shape)
                depth_p[self.mask > 0.5] = distances_p
                depth[depth_p < self.pelvis_dist] = 0
            else:
                warnings.warn("No pelvis found, returning depth without pelvis exclusion.")

        return depth

    def _segment_pelvis(self, noise_ml=2.5):
        """
        Segments the renal pelvis from the kidney mask and saves it as an
        attribute of the object.

        Parameters
        ----------
        noise_ml : float, optional
            Default 2.5
            Volume of noise in the mask to remove, in cubic centimeters (ml)
        """
        fill_vox = int(self.fill_ml / (np.prod(self.zoom) / 1000))
        mask_filled = remove_small_holes(self.mask, fill_vox)
        mask_ch = convex_hull_objects(mask_filled)
        hulls = (mask_ch ^ mask_filled) & mask_ch
        hulls = opening(hulls)

        noise_vox = int(noise_ml / (np.prod(self.zoom) / 1000))
        self.pelvis = remove_small_objects(hulls, noise_vox)

    @staticmethod
    def _normalise_data(data, mask=None):
        """
        Normalise data to have a mean of 0 and a standard deviation of 1.

        Parameters
        ----------
        data : numpy.ndarray
            Data to be normalised

        Returns
        -------
        numpy.ndarray
            Normalised data
        """
        if mask is None:
            mask = np.ones(data.shape, dtype=bool)
        return (data - np.nanmean(data[mask])) / np.nanstd(data[mask])


def load_pickle(fname):
    """
    Load a pickle file containing a QLayers object.

    Parameters
    ----------
    fname : str
        Filename of the pickle file to load.

    Returns
    -------
    QLayers
        The QLayers object stored in the pickle file.
    """
    with open(fname, "rb") as f:
        qlayers = pickle.load(f)
    return qlayers
