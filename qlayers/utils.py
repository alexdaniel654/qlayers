import nibabel as nib
import numpy as np
import skimage
from skimage.morphology import convex_hull_image


def convex_hull_objects(mask):
    """
    Compute convex hull of objects in a binary image. Each object is a
    different integer label in the input image.

    Note: If the convex hulls of multiple objects overlap, only the highest
    numbered object will be returned for the overlapping voxels.

    Parameters
    ----------
    mask : ndarray
        Binary input image.

    Returns
    -------
    mask_ch : ndarray
        Binary image with convex hull of objects in `mask`.
    """
    mask_labelled = skimage.measure.label(mask)
    labels = np.unique(mask_labelled)[np.unique(mask_labelled) > 0]
    mask_ch = []
    for label in labels:
        sub_mask = mask_labelled == label
        sub_ch = convex_hull_image(
            sub_mask, offset_coordinates=False, include_borders=False
        )
        mask_ch.append(sub_ch)
    mask_ch = np.sum(np.array(mask_ch), axis=0) > 0
    return mask_ch


def pad_dimensions(map_img):
    """
    Pad the dimensions of a 2D Nifti image to 3D.

    Parameters
    ----------
    map_img : Nifti1Image
        2D input image.

    Returns
    -------
    padded : Nifti1Image
        Padded image.
    """
    data = map_img.get_fdata()
    data = np.expand_dims(data, 2)
    return nib.Nifti1Image(data, map_img.affine, map_img.header)
