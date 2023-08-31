import numpy as np
from scipy.integrate import cumulative_trapezoid


def alpha_w(f, T):
    """Calculate the fresh water attenuation in m-1.

    Parameters
    ----------
    f : scalar, array
        Sound frequency in Hertz.
    T : scalar, array
        Temperature in Celsius degrees.

    Returns
    -------
    scalar, array
        return the fresh water attenuation in m-1.

    Examples
    --------
    >>> f = 2e6
    >>> T = 25
    >>> print(alpha_w(f, T))

    References
    --------
    [1] Fisher, F. H., & Simmons, V. P. (1977). Sound absorption in sea water. The Journal of the Acoustical Society of America, 62(3), 558-564.

    """
    return (
        10 ** (-15)
        * (55.9 - 2.37 * T + 4.77 * 10 ** (-2) * T**2 - 3.48 * 10 ** (-4) * T**3)
        * f**2
    )


def sound_velocity(T):
    """Compute sound velocity in fresh water in [m/s]

    Parameters
    ----------
    T : scalar, array
        Temperature in Celsius degrees.

    Returns
    -------
    scalar, array
        sound velocity in fresh water in [m/s]

    Examples
    --------
    >>> T = 25
    >>> print(sound_velocity(T))

    References
    --------
    [1] Lubbers, J., & Graaff, R. (1998). A simple and accurate formula for the sound velocity in water. Ultrasound in medicine & biology, 24(7), 1065-1068.


    """
    return 1405.03 + 4.624 * T - 3.83 * 10 ** (-2) * T**2


def near_field_theoretical(r, rn):
    """Return the near field function.

    Parameters
    ----------
    r : scalar, array
        Radial coordinate.
    rn : scalar, array
        Near field distance

    Returns
    -------
    scalar, array
        return the near field function.

    Examples
    --------
    >>> r = np.linspace(0.1, 1, 100)
    >>> rn = 0.2
    >>> Psi = near_field_theoretical(r, rn)

    References
    --------
    [1] Downing, A., Thorne, P. D., & Vincent, C. E. (1995). Backscattering from a suspension in the near field of a piston transducer. The Journal of the Acoustical Society of America, 97(3), 1614-1620.
    [2] Pedocchi, F., & García, M. H. (2012). Acoustic measurement of suspended sediment concentration profiles in an oscillatory boundary layer. Continental Shelf Research, 46, 87-95.


    """
    # return 1 + rn**3.2/(0.43*r*rn**2.2 + 0.48*r**3.2)
    return 1 + rn**3.2 / (1.35 * r * rn**2.2 + (2.5 * r) ** 3.2)


def correction_factor(r, alpha_w, Ks, Kt, psi):
    r"""Compute the part of the mean squared voltage independet of the sediment concentration: :math:`\frac{K_{\textup{s}} K_{\textup{t}}}{\psi r}^{2}e^{-4r\alpha_{\textup{w}}}`

    Parameters
    ----------
    r : scalar, array
        radial coordinate
    alpha_w : scalar, array
        water attenuation coefficient
    Ks : scalar, array
        Sediment constant.
    Kt : scalar, array
        Transducer constant.
    psi : scalar, array
        near field function.

    Returns
    -------
    scalar, array
        return the part of the mean squared voltage that is independent of the sediment concentration.

    """
    return np.exp(-4 * r * alpha_w) * ((Ks * Kt) / (r * psi)) ** 2


def create_MSvoltage(C, r, Xi, alpha_w, Ks, Kt, psi):
    r"""Compute the mean squared voltage from a concentration profile, given a set of acoustic coefficients: :math:`\langle V^{2} \rangle = C \frac{K_{\textup{s}} K_{\textup{t}}}{\psi r}^{2}e^{-4r\alpha_{\textup{w}} - \int_{0}^{r}\xi C dr}`

    Parameters
    ----------
    C : scalar, array
        Sediment concentration.
    r : scalar, array
        Radial coordinate.
    Xi : scalar, array
        sediment attenuation constant
    alpha_w : scalar, array
        water attenuation coefficient
    Ks : scalar, array
        Sediment constant.
    Kt : scalar, array
        Transducer constant.
    psi : scalar, array
        near field function.


    Returns
    -------
    scalar, array
        Mean squared voltage.
    >>>

    """
    if Xi == 0:
        integral = 0
    else:
        integral = cumulative_trapezoid(Xi * C + 0 * r, r + 0 * C, initial=0)
    factor = correction_factor(r, alpha_w, psi, Ks, Kt)
    return C * np.exp(-4 * integral) * factor
