import nibabel as nib
import numpy as np
import skimage
from skimage.morphology import convex_hull_image


def convex_hull_objects(mask):
    """
    Compute convex hull of objects in a binary image.

    This is essentially a combination of
    skimage.morphology.convex_hull_image (which only supports a single
    object but works in 3D) and skimage.morphology.convex_hull_object (which
    only works in 2D but supports multiple objects).

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
    if map_img.ndim == 2:
        data = map_img.get_fdata()
        data = np.expand_dims(data, 2)
        return nib.Nifti1Image(data, map_img.affine, map_img.header)
    elif map_img.ndim == 3:
        return map_img
    else:
        raise ValueError("Input image must be 2D or 3D")
