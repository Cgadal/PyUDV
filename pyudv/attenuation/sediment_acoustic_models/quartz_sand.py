import numpy as np
import scipy.integrate as scipyint

# def Chi_i(k, Rp):
#     x = k*Rp
#     return (0.09*x**4/(1380 + 560*x**2 + 150*x**4))


def Chi_i(k, Rp):
    x = k * Rp
    return (1.1 * (4/3) * 0.18 * x**4) / (1 + 1.3 * x**2 + (4/3) * 0.18 * x**4)


def Chi(k, Rp):
    if np.isscalar(Rp):
        return Chi_i(k, Rp)
    else:
        scipyint.trapezoid(Rp[0] * Rp[1], Rp[0]) * scipyint.trapezoid(
            Rp[0]**2 * Chi_i(k, Rp[0]) * Rp[1], Rp[0]) / scipyint.trapezoid(
                Rp[0]**3 * Rp[1], Rp[0])


def Xi(k, Rp, rho=2.65 * 1e3):
    if np.isscalar(Rp):
        return 3 * Chi(k, Rp) / (4*Rp*rho)
    else:
        a0 = scipyint.trapezoid(Rp[0] * Rp[1], Rp[0])  # mean value
        return 3 * Chi(k, Rp) / (4*a0*rho)
