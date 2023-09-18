import numpy as np
import scipy.integrate as scipyint
import scipy.interpolate as scipyinter

from pyudv.attenuation.direct_models import correction_factor

r"""
Explicit inversion scheme --
"""


def _explicit_solution_form(f, B, Xi, r):
    integral = scipyint.cumulative_trapezoid(f, r, initial=0)
    return f / (B - 4 * Xi * integral), f


def explicit_inversion(
    MSV: np.ndarray,
    r: np.ndarray,
    Xi: float,
    alpha_w: float,
    psi: np.ndarray,
    C0: float | np.ndarray,
    r0: float,
    delta_r: float = None,
    n_interp: int = 1000,
    kind: str = "cubic",
) -> tuple[np.ndarray, np.ndarray]:
    r"""
    Use the following explicit solution to infer particle concentration from Mean Square Voltage (MSV, :math:`\langle V^{2} \rangle`) data:

    .. math::

        C(r) = \frac{f}{B - 4\xi\displaystyle\int^{r}f},

    where

    .. math::

        f(r) = \langle V^{2} \rangle \left(\frac{r\psi_{(r)}}{K_{\rm s} K_{\rm t}}\right) e^{4r\alpha_{\rm w}} = Ce^{-4\alpha_{\rm s}}.

    Parameters
    ----------
    MSV : np.ndarray
        Input array containing :math:`\langle V^{2} \rangle` data. If 1D, must be the same size as `r`. If 2D, its last dimension must match the size of `r`.
    r : np.ndarray
        Radial distance to the probe array. Must be 1D, and match the last dimension of `MSV`.
    Xi : float
        sediment asborption coefficient.
    alpha_w : float
        water absorption coefficient.
    psi : np.ndarray
        Near-field function. Must be 1D, and match the last dimension of `MSV`.
    C0 : float | np.ndarray
        Reference concentration. Can be a float, or a 1D array matching the first dimension of `MSV`.
    r0 : float
        Reference distance corresponding to `C0`.
    delta_r : float, optional
        If not None, the reference concentration is then interpreted as an average between :math:`r-\delta r` and :math:`r+\delta r` (by default None).
    n_interp : int, optional
        Number of point used in the interpolation of the input signal, by default 1000. It is advise to vary this parameter to find its lowest value above which the result becomes independant of it.
    kind : str, optional
        Interpolation method used, by default "cubic". Refer to the documentation of :func:`interp1d <scipy.interpolate.interp1d>`.

    Returns
    -------
    C : np.ndarray
        Inferred concentration.
    f : np.ndarray
        Integrated array during the procedure. Returned to check the order of magnitudes over wich the integration is performed, and check possible numerical precision errors.
    """

    r_expanded = np.expand_dims(
        np.squeeze(r), axis=tuple(np.arange(len(MSV.shape[:-1])))
    )
    psi_expanded = np.expand_dims(
        np.squeeze(psi), axis=tuple(np.arange(len(MSV.shape[:-1])))
    )
    C0_expanded = C0 if np.isscalar(C0) else np.squeeze(C0)[..., None]
    #
    # ## calculating MSV reduced of factor
    factor = correction_factor(r_expanded, alpha_w, 1, 1, psi_expanded)
    f = MSV / factor
    interp_f = scipyinter.interp1d(r, f, axis=-1, kind=kind)
    # ## calculating adjustable constant with imposed point
    r_new = np.linspace(r.min(), r0, n_interp)
    if delta_r is None:  # if single fixed point
        constant_int = scipyint.trapezoid(interp_f(r_new), r_new)[..., None]
        constant_MSV = interp_f(r0)[..., None]
        B = constant_MSV / C0_expanded + 4 * Xi * constant_int
    else:  # if averaged within delta_r
        r_new = np.linspace(r.min(), r0 + delta_r, n_interp)
        constant_int_pos = scipyint.trapezoid(interp_f(r_new), r_new)[..., None]
        #
        r_new = np.linspace(r.min(), r0 - delta_r, n_interp)
        constant_int_neg = scipyint.trapezoid(interp_f(r_new), r_new)[..., None]
        #
        B = (
            4
            * Xi
            * (
                constant_int_pos
                - np.exp(-8 * delta_r * Xi * C0_expanded) * constant_int_neg
            )
            / (1 - np.exp(-8 * delta_r * Xi * C0_expanded))
        )
    return _explicit_solution_form(f, B, Xi, r_expanded)


#
# def explicit_inversion(MSV, r, Xi, alpha_w, psi, C0, r0):
#     r_expanded = np.expand_dims(r, axis=tuple(np.arange(len(MSV.shape[:-1]))))
#     psi_expanded = np.expand_dims(psi, axis=tuple(np.arange(len(MSV.shape[:-1]))))
#     C0_expanded = np.expand_dims(C0, axis=-1)
#     #
#     factor = correction_factor(r_expanded, alpha_w, 1, 1, psi_expanded)
#     MSV_reduced = MSV/factor
#     #
#     integral = cumulative_trapezoid(MSV_reduced, r_expanded, initial=0)
#     if np.isscalar(r0):
#         ind_r0 = np.arange(r.size)[r <= r0][-1]
#         constant_int = trapezoid(MSV_reduced[..., :ind_r0], r_expanded[..., :ind_r0])
#         return C0_expanded*MSV_reduced/(MSV_reduced[..., ind_r0:ind_r0 + 1] - C0_expanded*4*Xi*(integral - constant_int[..., None]))
#     else:
#         r_temp = np.broadcast_to(r_expanded, MSV_reduced.shape)
#         MSV_zeroed = np.where(r_temp <= r0[..., None], MSV_reduced, 0)
#         del r_temp
#         constant_int = trapezoid(MSV_zeroed, r_expanded)
#         del MSV_zeroed
#         #
#         ind_r0 = np.argmin(np.abs(r0[..., None] - r_expanded), axis=-1)
#         constant_MSV = np.take_along_axis(MSV_reduced, ind_r0[..., None], axis=-1)
#         del ind_r0
#         return C0_expanded*MSV_reduced/(constant_MSV - C0_expanded*4*Xi*(integral - constant_int[..., None]))
