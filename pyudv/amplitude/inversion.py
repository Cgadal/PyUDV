import numpy as np
import scipy.integrate as scipyint
from .direct_models import correction_factor

r"""
Explicit inversion scheme --
"""


def _explicit_solution_form(f, B, Xi, r):
    integral = scipyint.cumulative_trapezoid(f, r, initial=0)
    return f / (B - 4 * Xi * integral)


def explicit_inversion(MSV, r, Xi, alpha_w, psi, C0, r0, delta_r=None):
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
    # ## calculating adjustable constant with imposed point
    r0_brod = (
        np.broadcast_to(r0, f.shape)
        if np.isscalar(r0)
        else np.broadcast_to(r0.squeeze()[..., None], f.shape)
    )
    r_temp = np.broadcast_to(r_expanded, f.shape)
    if delta_r is None:  # if single fixed point
        f_zeroed = np.where(
            r_temp <= r0_brod, f, 0
        )  # putting 0 everywhere r < r0 to compute constant integrals
        constant_int = scipyint.trapezoid(f_zeroed, r_expanded)[..., None]
        del f_zeroed
        ind_r0 = np.argmin(np.abs(r0_brod - r_expanded), axis=-1)
        constant_MSV = np.take_along_axis(f, ind_r0[..., None], axis=-1)
        B = constant_MSV / C0_expanded + 4 * Xi * constant_int
    else:  # if averaged within delta_r
        f_zeroed = np.where(
            r_temp <= r0_brod + delta_r, f, 0
        )  # putting 0 everywhere r < r0 + delta_r to compute constant integrals
        constant_int_pos = scipyint.trapezoid(f_zeroed, r_expanded)[..., None]
        f_zeroed = np.where(
            r_temp <= r0_brod - delta_r, f, 0
        )  # putting 0 everywhere r < r0 - delta_r to compute constant integrals
        constant_int_neg = scipyint.trapezoid(f_zeroed, r_expanded)[..., None]
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
