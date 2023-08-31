import numpy as np
from scipy.integrate import trapezoid


def Chi_i(k, Rp):
    x = k * Rp
    phi = 1 - 0.4 * np.exp(((x - 5.5) / 0.7) ** 2)
    return (0.24 * phi * x**4) / (
        0.7 + 0.3 * x + 2.1 * x**2 - 0.7 * x**3 + 0.3 * x**4
    )


def Chi(k, Rp):
    if np.isscalar(Rp):
        return Chi_i(k, Rp)
    else:
        trapezoid(Rp[0] * Rp[1], Rp[0]) * trapezoid(
            Rp[0] ** 2 * Chi_i(k, Rp[0]) * Rp[1], Rp[0]
        ) / trapezoid(Rp[0] ** 3 * Rp[1], Rp[0])


def Xi(k, Rp):
    if np.isscalar(Rp):
        return 3 * Chi(k, Rp) / (4 * Rp)
    else:
        a0 = trapezoid(Rp[0] * Rp[1], Rp[0])  # mean value
        return 3 * Chi(k, Rp) / (4 * a0)
