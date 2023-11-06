"""
This module helps you deal with probe geometries.
"""

import matplotlib.axes as mplaxes
import matplotlib.colors as mpcolor
import matplotlib.patches as mptch
import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
import scipy.interpolate as scipyinterp

import pyudv.helpers as ph


class Transducer:
    """
     A class used to caracterize a transducer physical properties.

    Parameters
    ----------
    activediameter : float, optional
        transducer active diameter, by default None
    overalldiameter : float, optional
        transducer overall diameter, by default None
    overalllength : float, optional
        transducer overall length, by default None
    nearfield : float, optional
        transducer near-field distance, by default None
    halfangle : float, optional
        transducer divegence half-angle (deg.), by default None
    freq : float, optional
        transducer centre frequency, by default None
    """

    def __init__(self,
                 activediameter: float = None,
                 overalldiameter: float = None,
                 overalllength: float = None,
                 nearfield: float = None,
                 halfangle: float = None,
                 freq: float = None) -> None:
        """
        Parameters
        ----------
        activediameter : float, optional
            transducer active diameter, by default None
        overalldiameter : float, optional
            transducer overall diameter, by default None
        overalllength : float, optional
            transducer overall length, by default None
        nearfield : float, optional
            transducer near-field distance, by default None
        halfangle : float, optional
            transducer divegence half-angle (deg.), by default None
        freq : float, optional
            transducer centre frequency, by default None
        """

        self.activediameter = activediameter  #: transducer active diameter
        self.overalldiameter = overalldiameter  #: transducer overall diameter
        self.overalllength = overalllength  #: transducer overall length
        self.nearfield = nearfield  #: transducer near-field distance
        self.halfangle = halfangle  #: transducer beam half divergence angle
        self.freq = freq  #: transducer frequency
        #
        self.width_at_nearfield = self.beam_width(
            self.nearfield)  #: beam width at near-field distance

    def beam_width(self, distance: float) -> float:
        """
        Compute the beam width at a certain distance from the beam half divergence angle.

        Parameters
        ----------
        distance : float
            distance

        Returns
        -------
        float
            beam width at a distance
        """
        return self.activediameter + 2 * np.tan(np.radians(
            self.halfangle)) * distance


# registering metflow probes
_MetFlow05Mhz = Transducer(freq=0.5,
                           activediameter=19,
                           overalldiameter=23,
                           overalllength=60,
                           nearfield=30.5,
                           halfangle=4.6)

_MetFlow1Mhz = Transducer(freq=1,
                          activediameter=13,
                          overalldiameter=16,
                          overalllength=60,
                          nearfield=28.5,
                          halfangle=3.4)

_MetFlow2Mhz = Transducer(freq=2,
                          activediameter=10,
                          overalldiameter=13,
                          overalllength=60,
                          nearfield=33.7,
                          halfangle=2.2)

_MetFlow4Mhz = Transducer(freq=4,
                          activediameter=5,
                          overalldiameter=8,
                          overalllength=60,
                          nearfield=16.9,
                          halfangle=2.2)
_MetFlow8Mhz = Transducer(freq=4,
                          activediameter=2.5,
                          overalldiameter=8,
                          overalllength=60,
                          nearfield=8.43,
                          halfangle=2.2)

MetFlowTransducers = {
    0.5: _MetFlow05Mhz,
    1: _MetFlow1Mhz,
    2: _MetFlow2Mhz,
    4: _MetFlow4Mhz,
    8: _MetFlow8Mhz
}
"""A dictionnary containing MetFlow :class:`Transducer <pyudv.geometry.Transducer>` object. All distances are stored in cm."""

# actually putting distances in cm
_distance_attr = [
    'activediameter', 'overalldiameter', 'overalllength', 'nearfield',
    'width_at_nearfield'
]
for transducer in MetFlowTransducers.values():
    for attr in _distance_attr:
        setattr(transducer, attr, getattr(transducer, attr) * 1e-2)


class Probe:
    """
     A class used to represent a single probe.

    Parameters
    ----------
    r : ndarray, (N,)
        a numpy array containing the radial coordinates of the beam.
    alpha : float
        orientation of the probe, in degree in the trigonometric referential.
    Pref : list[float, ndarray]
        A list a size 2 defining a reference point by containing its radial coordinate, and a numpy array containing its coordinate in the "real" reference space.
    transducer : Transducer, optional
        :class:`Transducer <pyudv.geometry.Transducer>` object giving the physical probe properties, by default None
    """

    def __init__(self,
                 r: np.ndarray,
                 alpha: float,
                 Pref: list[float, np.ndarray],
                 transducer: Transducer = None) -> None:
        """
        Parameters
        ----------
        r : ndarray
            a numpy array containing the radial coordinates of the beam.
        alpha : float
            orientation of the probe, in degree in the trigonometric referential.
        Pref : list[float, ndarray]
            A list a size 2 defining a reference point by containing its radial coordinate, and a numpy array containing its coordinate in the "real" reference space.
        transducer : Transducer, optional
            :class:`Transducer <pyudv.geometry.Transducer>` object giving the physical probe properties, by default None
        """

        self.r = r  #: vector along beam axis
        self.alpha = alpha  #: beam axis inclination, trigo referential
        self.Pref = Pref  #: reference point [r_ref, (x_ref, y_ref)]
        #
        self.r_ref = self.Pref[0]  #: reference point radial coordinats
        self.X_ref = self.Pref[1]  #: reference point real coordinates
        self.unit_vec = np.array([
            np.cos(np.radians(alpha)),
            np.sin(np.radians(alpha))
        ])  #: unit vector defining the beam direction in real coordinates
        self.unit_vec_perp = np.array(
            [-np.sin(np.radians(alpha)),
             np.cos(np.radians(alpha))]
        )  #: unit vector defining the beam perpendicular direction in real coordinates
        self.O = (Pref[1] - Pref[0] * self.unit_vec
                  )  #: coordinates of the origin of the beam
        self.E = (self.X_ref + (r.max() - self.r_ref) * self.unit_vec
                  )  #: coordinates of the end of the beam
        self.transducer = transducer
        # Beam coordinates
        self.x = self.r_to_x(self.r)
        self.z = self.r_to_z(self.r)

    def r_to_z(self, r: np.ndarray | float) -> np.ndarray | float:
        """
        Convert radial coordinates of the probe into z coordinates

        Parameters
        ----------
        r : np.ndarray | float
            radial coordinate

        Returns
        -------
        np.ndarray | float
            corresponding z coordinate
        """
        return self.O[1] + r * self.unit_vec[1]

    def r_to_x(self, r: np.ndarray | float) -> np.ndarray | float:
        """
        Convert radial coordinates of the probe into x coordinates

        Parameters
        ----------
        r : np.ndarray | float
            corresponding x coordinate

        Returns
        -------
        np.ndarray | float
            corresponding z coordinate
        """
        return self.O[0] + r * self.unit_vec[0]

    def xz_to_r(self, xz: np.ndarray, axis=-1) -> np.ndarray:
        """
        Convert xz coordinates into probe radial coordinates

        Parameters
        ----------
        xz : np.ndarray
            xz coordinate
        axis : int, optional
            axis along which the norm is calculated, by default -1

        Returns
        -------
        np.ndarray
            radial coordinates
        """
        return np.linalg.norm(xz - self.O, axis=axis)

    def plot_probe(self,
                   ax: mplaxes.Axes,
                   color: str = None,
                   plot_transducer: bool = True):
        """
        Sketch the probe and its beam on a given matplotlib axe.

        Parameters
        ----------
        ax : mplaxes.Axes
            matplotlib axe on which the probe is sketched.
        color : str, optional
            probe color, by default None
        plot_transducer: bool, optional
            whether or not to plot the transducer physical properties, by default True. Only works if a :class:`Transducer <pyudv.geometry.Transducer>` object has been passed to the :class:`Probe <pyudv.geometry.Probe>` object during its creation.
        """
        a = ax.scatter(self.O[0], self.O[1], color=color)
        color = a.get_facecolor()
        if (self.transducer is not None) & plot_transducer:
            # ### plot transducer
            a.remove()
            height = self.transducer.overalldiameter
            width = self.transducer.overalllength
            angle = self.alpha
            xy = self.O - np.array([width, height / 2])
            rect = plt.Rectangle(xy,
                                 width,
                                 height,
                                 angle=angle,
                                 rotation_point=tuple(self.O),
                                 facecolor=mpcolor.to_rgba(color, 0.2),
                                 edgecolor=color)
            ax.add_patch(rect)
            # #### plot beam
            max_width = self.transducer.beam_width(
                np.linalg.norm(self.E - self.O))
            xy = np.array([
                self.O -
                self.unit_vec_perp * self.transducer.activediameter / 2,
                self.E - self.unit_vec_perp * max_width / 2,
                self.E + self.unit_vec_perp * max_width / 2,
                self.O +
                self.unit_vec_perp * self.transducer.activediameter / 2,
            ])
            poly = plt.Polygon(xy,
                               edgecolor=color,
                               facecolor=mpcolor.to_rgba(color, 0.1),
                               linewidth=1)
            ax.add_patch(poly)
            # #### plot nearfield
            p = self.O + self.unit_vec * self.transducer.nearfield
            xy = np.array([
                self.O -
                self.unit_vec_perp * self.transducer.activediameter / 2,
                p -
                self.unit_vec_perp * self.transducer.width_at_nearfield / 2,
                p +
                self.unit_vec_perp * self.transducer.width_at_nearfield / 2,
                self.O +
                self.unit_vec_perp * self.transducer.activediameter / 2,
            ])
            poly = plt.Polygon(xy, color=color)
            ax.add_patch(poly)
        else:
            ax.annotate(
                "",
                xytext=self.O,
                xy=self.O + np.linalg.norm(
                    (self.E - self.O)) / 10 * self.unit_vec,
                arrowprops=dict(arrowstyle="->", color=color),
            )
        # plot unused beam part
        (a, ) = ax.plot([self.O[0], self.x[0]], [self.O[1], self.z[0]],
                        color=color,
                        ls="--",
                        lw=1,
                        alpha=0.3)
        # plot used beam part
        (a, ) = ax.plot(
            [self.x[0], self.E[0]],
            [self.z[0], self.E[1]],
            color=color,
            ls="--",
            lw=1,
        )


def sketch_probes(
    probes: list[Probe],
    combinations: list[tuple[int, int]] = None,
    ax: mplaxes.Axes = None,
    probe_colors: list[str] = None,
    combination_colors: list[str] = None,
):
    """
    Sketch multiple probes, and optionaly their intersections, etc ...

    Parameters
    ----------
    probes : list[Probe]
        list of :class:`Probe <pyudv.geometry.Probe>` object to be sketched.
    combinations : list[tuple[int, int]], optional
        list of tuple of size 2 containing indexes of `probes` whose intersection should be sketched.
    ax : mplaxes.Axes, optional
        matplotlib axe on which the sketch is drawn. If None, a new figure is created (default is None).
    probe_colors : list[str], optional
        list of colors for the probes, by default None. It should be the same size as `probes`.
    combination_colors : list[str], optional
        list of colors for the intersetctions, by default None. It should be the same size as `combinations`.
    """
    if ax is None:
        fig, ax = plt.subplots(1, 1, layout="constrained")
    ax.set_aspect("equal")
    #
    if probe_colors is None:
        probe_colors = [None for i in probes]
    if (combinations is not None) and (combination_colors is None):
        combination_colors = [None for i in combinations]
    #
    # plot probes
    for probe, color in zip(probes, probe_colors):
        probe.plot_probe(ax, color=color)
    if combinations is not None:
        for combination, combination_color in zip(combinations,
                                                  combination_colors):
            probe1, probe2 = probes[combination[0]], probes[combination[1]]
            X = probe_crossing_point(probe1, probe2)
            ax.scatter(X[0], X[1], color=combination_color)
            #
            z_interp, z_min, z_max = compute_vertical_axis(probe1, probe2)
            ax.plot([X[0], X[0]], [z_min, z_max], color=combination_color)


def probe_crossing_point(probe1: Probe, probe2: Probe) -> np.ndarray:
    """
    calculate the crossing point of two probe beams.

    Parameters
    ----------
    probe1 : Probe
        First :class:`Probe <pyudv.geometry.Probe>` object.
    probe2 : Probe
        Second :class:`Probe <pyudv.geometry.Probe>` object.

    Returns
    -------
    ndarray, (2. )
        coordinates of the crossing point.
    """
    M = np.array([-probe1.unit_vec, probe2.unit_vec]).T
    S = probe1.O - probe2.O
    T = np.linalg.inv(M) @ S
    X = probe1.O + T[0] * probe1.unit_vec
    return X


def compute_vertical_axis(probe1: Probe, probe2: Probe):
    """
    compute the common vertical axis between two probe beams.

    Parameters
    ----------
    probe1 : Probe
        First :class:`Probe <pyudv.geometry.Probe>` object.
    probe2 : Probe
        Second :class:`Probe <pyudv.geometry.Probe>` object.

    Returns
    -------
    z_interp : ndarray
        numpy array containg the vertical coordinates of the common axis.
    z_min : float
        minimal vertical coordinate
    z_max : float
        maximal vertical coordinate
    """
    return ph.compute_common_part(probe1.z, probe2.z)


def reconstruct_velocity(
    u1: npt.ArrayLike, u2: npt.ArrayLike, probe1: Probe, probe2: Probe
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray,
           np.ndarray]:
    """
    This function takes the velocities measured by two probes, and reconstruct the velocity field, assuming that it only depends on the vertical coordinate.

    Parameters
    ----------
    u1 : ArrayLike
        velocity vector measured by first probe. The first dimension should be the same dimension as the radial coordinate contained in `probe1.z`.
    u2 : ArrayLike
        velocity vector measured by second probe. The first dimension should be the same dimension as the radial coordinate contained in `probe2.z`.
    probe1 : Probe
        First :class:`Probe <pyudv.geometry.Probe>` object.
    probe2 : Probe
        Second :class:`Probe <pyudv.geometry.Probe>` object.

    Returns
    -------
    U : ndarray
        reconstructed velocity components in the coordinate system corresponding to the one of the reference points of the probes.
    z_interp : ndarray
        coordinate vector corresponding to `U`
    X : ndarray
        crossing point of the probes
    dx_1 : ndarray
        horizontal distance between the first probe beam, and the vertical axis passing by the crossing point `X`
    dx_2 : ndarray
        horizontal distance between the second probe beam, and the vertical axis passing by the crossing point `X`

    """
    # # Build common vertical axis
    z_interp, _, _ = compute_vertical_axis(probe1, probe2)
    # # interpolating velocity
    interp1 = scipyinterp.interp1d(probe1.z, u1, axis=0)
    interp2 = scipyinterp.interp1d(probe2.z, u2, axis=0)
    #
    u1_interp = interp1(z_interp)
    u2_interp = interp2(z_interp)
    #
    # velocity reconstruction
    M = np.array([probe1.unit_vec, probe2.unit_vec])
    # U = np.linalg.inv(M) @ np.array([u1_interp, u2_interp])
    U = np.einsum("ij, jk... -> ik...", np.linalg.inv(M),
                  np.array([u1_interp, u2_interp]))
    #
    # # other quantities
    # crossing point
    X = probe_crossing_point(probe1, probe2)
    # horizontal distances to the vertical axis passing by the crossing point
    dx_1 = np.tan(np.radians(probe1.alpha)) * (z_interp - X[1])
    dx_2 = np.tan(np.radians(probe2.alpha)) * (z_interp - X[1])
    return U, z_interp, X, dx_1, dx_2


def average_amplitude(
    a1: npt.ArrayLike, a2: npt.ArrayLike, probe1: Probe, probe2: Probe
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray,
           np.ndarray]:
    """
    This function takes the amplitude measured by two probes and average them, assuming that it only depends on the vertical coordinate.

    Parameters
    ----------
    a1 : npt.ArrayLike
        amplitude measured by first probe. The first dimension should be the same dimension as the radial coordinate contained in `probe1.z`.
    a2 : npt.ArrayLike
        amplitude measured by second probe. The first dimension should be the same dimension as the radial coordinate contained in `probe2.z`.
    probe1 : Probe
        First :class:`Probe <pyudv.geometry.Probe>` object.
    probe2 : Probe
        Second :class:`Probe <pyudv.geometry.Probe>` object.

    Returns
    -------
    A : ndarray
        averaged amplitud in the coordinate system corresponding to the one of the reference points of the probes.
    z_interp : ndarray
        coordinate vector corresponding to `U`
    X : ndarray
        crossing point of the probes
    dx_1 : ndarray
        horizontal distance between the first probe beam, and the vertical axis passing by the crossing point `X`
    dx_2 : ndarray
        horizontal distance between the second probe beam, and the vertical axis passing by the crossing point `X`

    """
    # # Build common vertical axis
    z_interp, _, _ = compute_vertical_axis(probe1, probe2)
    # # interpolating velocity
    interp1 = scipyinterp.interp1d(probe1.z, a1, axis=0)
    interp2 = scipyinterp.interp1d(probe2.z, a2, axis=0)
    #
    a1_interp = interp1(z_interp)
    a2_interp = interp2(z_interp)
    #
    # # other quantities
    # crossing point
    X = probe_crossing_point(probe1, probe2)
    # horizontal distances to the vertical axis passing by the crossing point
    dx_1 = np.tan(np.radians(probe1.alpha)) * (z_interp - X[1])
    dx_2 = np.tan(np.radians(probe2.alpha)) * (z_interp - X[1])
    return (a1_interp+a2_interp) / 2, z_interp, X, dx_1, dx_2


if __name__ == "__main__":
    #
    # # Test Probe object definition
    r = np.linspace(0, 5, 100)
    alpha1, alpha2 = -120, -70  # deg
    O1, O2 = np.array([1, 8]), np.array([-1, 7])
    probe1_pars = [r, alpha1, [0, O1]]
    probe2_pars = [r, alpha2, [0, O2]]
    #
    probe1 = Probe(*probe1_pars, transducer=MetFlowTransducers[4])
    probe2 = Probe(*probe2_pars)
    #
    sketch_probes(
        [probe1, probe2],
        combinations=[[0, 1]],
        probe_colors=[None, None],
        combination_colors=["k"],
    )

    # #### Test velocity reconstruction
    def U(z):
        u = 5 * (5 - z)**2
        v = u / 10
        U = np.array([u, v])
        return U

    #
    # make fake signals
    u1 = U(probe1.z).T @ probe1.unit_vec
    u2 = U(probe2.z).T @ probe2.unit_vec
    #
    # #### reconstruction
    U_rec, z_interp, X, dx_1, dx_2 = reconstruct_velocity(
        u1, u2, probe1, probe2)
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
