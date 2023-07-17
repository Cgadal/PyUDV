from scipy.ndimage import uniform_filter1d, uniform_filter
import os


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
    c1 = uniform_filter1d(arr, windows, mode='constant', axis=axis, **kwargs)
    c2 = uniform_filter1d(arr*arr, windows, mode='constant', axis=axis, **kwargs)
    return ((c2 - c1*c1)**.5)


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
    c1 = uniform_filter(arr, size, mode='constant', **kwargs)
    c2 = uniform_filter(arr*arr, size, mode='constant', **kwargs)
    return ((c2 - c1*c1)**.5)


def create_arboresence(path):
    # Create all arborescence present in path if directories does not exist. Otherwise leave them unaltered.
    if not os.path.exists(path):
        if path[-1] != os.sep:
            path += os.sep
        os.makedirs(path[:path.rindex(os.path.sep)], exist_ok=True)
