import os

import numpy as np
import numpy.typing as npt
from scipy.ndimage import uniform_filter, uniform_filter1d


def moving_std1d(arr, windows, axis=-1, **kwargs):
    """Rolling standard deviation on 1 dimension of the input array.

    Parameters
    ----------
    arr : array
        Input array
    windows : int
        windows over wich the standard deviation is calculated.
    axis : int
        axis along which the standard deviation is calculated. (the default is -1).
    **kwargs :
        Optional parameters passed to `scipy.ndimage.uniform_filter1d`

    Returns
    -------
    array
        Array of standard deviation

    """
    # based on https://stackoverflow.com/a/18422519/9530017
    c1 = uniform_filter1d(arr, windows, mode="constant", axis=axis, **kwargs)
    c2 = uniform_filter1d(arr * arr, windows, mode="constant", axis=axis, **kwargs)
    return (c2 - c1 * c1) ** 0.5


def moving_std(arr, size, **kwargs):
    """Rolling standard deviation on n-dimensional arrays.

    Parameters
    ----------
    arr : array
        Input array
    size : int
        The sizes of the uniform filter are given for each axis as a sequence, or as a single number, in which case the size is equal for all axes.
    **kwargs :
        Optional parameters passed to `scipy.ndimage.uniform_filter`

    Returns
    -------
    array
        Array of standard deviation

    """
    # based on https://stackoverflow.com/a/18422519/9530017
    c1 = uniform_filter(arr, size, mode="constant", **kwargs)
    c2 = uniform_filter(arr * arr, size, mode="constant", **kwargs)
    return (c2 - c1 * c1) ** 0.5


def compute_common_part(
    vec1: np.ndarray, vec2: np.ndarray
) -> tuple[np.ndarray, float, float]:
    """
    Compute the common part between two vector, and return it. If it is not exactly colocated, it returns a interpolation with the same number of points betwen the bounds common to the two arrays.

    Parameters
    ----------
    vec1 : np.ndarray
        first 1D array.
    vec2 : np.ndarray
        second 1D array.

    Returns
    -------
    z_interp: np.ndarray
        interpolated vector
    z_min: float
        lower bound
    z_max: float
        higher bound
    """
    z_min = np.nanmax([np.nanmin(vec1), np.nanmin(vec2)])
    z_max = np.min([np.nanmax(vec1), np.nanmax(vec2)])
    n_pts = (
        ((vec1 <= z_max) & (vec1 >= z_min)).sum()
        + ((vec2 <= z_max) & (vec2 >= z_min)).sum()
    ) / 2
    z_interp = np.linspace(z_min, z_max, int(n_pts))
    return z_interp, z_min, z_max


def create_arboresence(path):
    # Create all arborescence present in path if directories does not exist. Otherwise leave them unaltered.
    if not os.path.exists(path):
        if path[-1] != os.sep:
            path += os.sep
        os.makedirs(path[: path.rindex(os.path.sep)], exist_ok=True)
