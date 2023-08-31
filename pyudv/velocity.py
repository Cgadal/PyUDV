import numpy as np
import numpy.typing as npt
import scipy.interpolate as scipyinterp

from pyudv.probes import Probe, compute_vertical_axis, probe_crossing_point


def reconstruct_velocity(
    u1: npt.ArrayLike,
    u2: npt.ArrayLike,
    probe1_pars: list[npt.ArrayLike, float, list[float, npt.ArrayLike]],
    probe2_pars: list[npt.ArrayLike, float, list[float, npt.ArrayLike]],
) -> tuple[
    npt.NDArray, npt.NDArray, npt.NDArray, npt.NDArray, npt.NDArray, npt.NDArray
]:
    """
    This function takes the velocities measured by two probes, and reconstruct the velocity field, assuming that it only depends on the vertical coordinate.

    Parameters
    ----------
    u1 : ArrayLike
        velocity vector measured by first probe
    u2 : ArrayLike
        velocity vector measured by second probe
    probe1_pars : list[ArrayLike, float, list[float, ArrayLike]]
        list containing parameters of the first probe
    probe2_pars : list[ArrayLike, float, list[float, ArrayLike]]
        list containing parameters of the second probe

    Returns
    -------
    U : NDArray
        reconstructed velocity components in the coordinate system corresponding to the one of the refernce points of the probes.
    z_interp : NDArray
        coordinate vector corresponding to `U`
    X : NDArray
        crossing point of the probes
    dx_1 : NDArray
        horizontal distance between the first probe beam, and the vertical axis passing by the crossing point `X`
    dx_2 : NDArray
        horizontal distance between the second probe beam, and the vertical axis passing by the crossing point `X`

    """
    # # init probes
    probe1 = Probe(*probe1_pars)
    probe2 = Probe(*probe2_pars)
    #
    # # Build common vertical axis
    z_interp, _, _ = compute_vertical_axis(probe1, probe2)
    # # interpolating velocity
    interp1 = scipyinterp.interp1d(probe1.z, u1)
    interp2 = scipyinterp.interp1d(probe2.z, u2)
    #
    u1_interp = interp1(z_interp)
    u2_interp = interp2(z_interp)
    #
    # velocity reconstruction
    M = np.array([probe1.unit_vec, probe2.unit_vec])
    U = np.linalg.inv(M) @ np.array([u1_interp, u2_interp])
    #
    # # other quantities
    # crossing point
    X = probe_crossing_point(probe1, probe2)
    # horizontal distances to the vertical axis passing by the crossing point
    dx_1 = np.tan(np.radians(probe1.alpha)) * (z_interp - X[1])
    dx_2 = np.tan(np.radians(probe2.alpha)) * (z_interp - X[1])
    return U, z_interp, X, dx_1, dx_2


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    def U(z):
        u = 5 * (5 - z) ** 2
        v = u / 10
        U = np.array([u, v])
        return U

    # define probes
    r = np.linspace(0, 5, 100)
    alpha1, alpha2 = -120, -70  # deg
    O1, O2 = np.array([1, 8]), np.array([-1, 7])
    probe1_pars = [r, alpha1, [0, O1]]
    probe2_pars = [r, alpha2, [0, O2]]
    #
    probe1 = Probe(*probe1_pars)
    probe2 = Probe(*probe2_pars)
    # make fake signals
    u1 = U(probe1.z).T @ probe1.unit_vec
    u2 = U(probe2.z).T @ probe2.unit_vec
    #
    # #### reconstruction
    U_rec, z_interp, X, dx_1, dx_2 = reconstruct_velocity(
        u1, u2, probe1_pars, probe2_pars
    )
    U_th = U(z_interp)
    #
    fig, axarr = plt.subplots(1, 2, layout="constrained", sharey=True)
    for ax, u_th, u_rec in zip(axarr, U_th, U_rec):
        ax.plot(u_th, z_interp, ".", label="base")
        ax.plot(u_rec, z_interp, ".", label="reconstructed")
        ax.legend()
    axarr[0].set_xlabel("u")
    axarr[1].set_xlabel("v")
    axarr[0].set_ylabel("z")
    plt.show()
    #
