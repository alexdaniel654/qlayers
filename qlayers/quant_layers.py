import nibabel as nib
import numpy as np
import pandas as pd
import trimesh
from nibabel.processing import resample_from_to
from skimage.measure import marching_cubes
from skimage.morphology import opening, remove_small_holes, remove_small_objects
from trimesh import smoothing

from .utils import convex_hull_objects, pad_dimensions


class QLayers:
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
        self.space = space
        self.depth = self._calculate_depth()
        self.layers = np.ceil(self.depth * (1 / self.thickness)) / (
            1 / self.thickness
        )
        self.layers_list = np.unique(self.layers)
        self.df_long = pd.DataFrame(
            columns=["depth", "layer", "measurement", "value"]
        )
        self.maps = []
        if self.space == "layers":
            self.df_wide = pd.DataFrame(columns=["depth", "layer"])
            self.df_wide["depth"] = self.depth[self.mask]
            self.df_wide["layer"] = self.layers[self.mask]

    def add_map(self, map_img, name):
        self.maps.append(name)
        if map_img.ndim == 2:
            map_img = pad_dimensions(map_img)

        if self.space == "layers":
            # Resample map into space of layers
            # Doing cval as a big and unusual number as cval=np.nan doesn't
            # work
            map_img = resample_from_to(
                map_img, self.mask_img, cval=np.pi * np.e * 1e10
            )
            map_data = map_img.get_fdata()
            map_data[map_data == np.pi * np.e * 1e10] = np.nan
            self.df_wide[name] = map_data[self.mask]
            sub_df = pd.DataFrame(
                columns=["depth", "layer", "measurement", "value"]
            )
            sub_df["depth"] = self.depth[self.mask]
            sub_df["layer"] = self.layers[self.mask]
            sub_df["measurement"] = name
            sub_df["value"] = map_data[self.mask]
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

            sub_df = pd.DataFrame(
                columns=["depth", "layer", "measurement", "value"]
            )
            sub_df["depth"] = depth_rs[mask_rs]
            sub_df["layer"] = layers_rs[mask_rs]
            sub_df["measurement"] = name
            sub_df["value"] = map_data[mask_rs]
            self.df_long = pd.concat([self.df_long, sub_df])

    def add_tissue(self, tissue_img):
        # TODO add_tissue method
        raise NotImplementedError("Not yet implemented")

    def get_df(self, format="long"):
        if format == "wide":
            if self.space == "map":
                raise NotImplementedError
            else:
                return self.df_wide
        elif format == "long":
            return self.df_long.dropna()
        else:
            raise NotImplementedError

    def get_depth(self):
        return self.depth

    def get_layers(self):
        return self.layers

    def remove_all_maps(self):
        self.df_long = pd.DataFrame(
            columns=["depth", "layer", "measurement", "value"]
        )
        self.maps = []
        if self.space == "layers":
            self.df_wide = pd.DataFrame(columns=["depth", "layer"])
            self.df_wide["depth"] = self.depth[self.mask]
            self.df_wide["layer"] = self.layers[self.mask]

    def save_depth(self, fname):
        depth_img = nib.Nifti1Image(self.depth, self.affine)
        nib.save(depth_img, fname)

    def save_layers(self, fname):
        layer_img = nib.Nifti1Image(self.layers, self.affine)
        nib.save(layer_img, fname)

    def save_pelvis(self, fname):
        if not hasattr(self, "pelvis"):
            self._segment_pelvis()
        pelvis_img = nib.Nifti1Image(self.pelvis.astype(np.int32), self.affine)
        nib.save(pelvis_img, fname)

    def _calculate_depth(self):
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
        print("Calculating Distances")
        (closest_points, distances, triangle_id) = mesh.nearest.on_surface(
            points
        )

        # Write these distances to voxels in the shape of the original image
        depth = np.zeros(self.mask.shape)
        depth[self.mask > 0.5] = distances
        if self.pelvis_dist != 0:
            print("Masking Pelvise")
            self._segment_pelvis()
            verts_p, faces_p, normals_p, _ = marching_cubes(
                self.pelvis.astype(np.uint8),
                spacing=self.zoom,
                level=0.5,
                step_size=1.0,
            )
            mesh_p = trimesh.Trimesh(
                vertices=verts_p, faces=faces_p, vertex_normals=normals_p
            )
            (
                closest_points_p,
                distances_p,
                triangle_id_p,
            ) = mesh_p.nearest.on_surface(points)
            depth_p = np.zeros(self.mask.shape)
            depth_p[self.mask > 0.5] = distances_p
            depth[depth_p < self.pelvis_dist] = 0

        return depth

    def _segment_pelvis(self):
        fill_vox = int(self.fill_ml / (np.prod(self.zoom) / 1000))
        mask_filled = remove_small_holes(self.mask, fill_vox)
        mask_ch = convex_hull_objects(mask_filled)
        hulls = (mask_ch ^ mask_filled) & mask_ch
        hulls = opening(hulls)
        noise_ml = 2.5
        noise_vox = int(noise_ml / (np.prod(self.zoom) / 1000))
        self.pelvis = remove_small_objects(hulls, noise_vox)
