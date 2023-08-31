import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
import matplotlib.axes as mplaxes


class Probe:
    def __init__(
        self,
        r: npt.ArrayLike,
        alpha: float,
        Pref: list[float, npt.ArrayLike],
    ):
        self.r = r  # vector along beam axis
        self.alpha = alpha  # beam axis inclination, trigo referential
        self.Pref = Pref  # refernce point [r_ref, (x_ref, y_ref)]
        #
        self.r_ref = self.Pref[0]
        self.X_ref = self.Pref[1]
        self.unit_vec = np.array(
            [np.cos(np.radians(alpha)), np.sin(np.radians(alpha))]
        )  # beam unit vector in real coordinates
        self.O = Pref[1] - Pref[0] * self.unit_vec  # Beam origin in real coordinates
        self.E = (
            self.X_ref + (r.max() - self.r_ref) * self.unit_vec
        )  # Beam end in real coordinates

        # Beam coordinates
        self.x = self.O[0] + self.r * self.unit_vec[0]
        self.z = self.O[1] + self.r * self.unit_vec[1]

    def plot_probe(self, ax: mplaxes.Axes, color=None):
        ax.scatter(self.O[0], self.O[1], color=color)
        (a,) = ax.plot(
            [self.O[0], self.E[0]],
            [self.O[1], self.E[1]],
            color=color,
            ls="--",
            lw=1,
        )
        ax.annotate(
            "",
            xytext=self.O,
            xy=self.O + (self.r.max() - self.r.min()) / 10 * self.unit_vec,
            arrowprops=dict(arrowstyle="->", color=a.get_color()),
        )


def sketch_probes(
    probes: list[Probe],
    combinations: list[tuple] = None,
    ax=None,
    probe_colors: list[str] = None,
    combination_colors: list[str] = None,
):
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
        for combination, combination_color in zip(combinations, combination_colors):
            probe1, probe2 = probes[combination[0]], probes[combination[1]]
            X = probe_crossing_point(probe1, probe2)
            ax.scatter(X[0], X[1], color=combination_color)
            #
            z_interp, z_min, z_max = compute_vertical_axis(probe1, probe2)
            ax.plot([X[0], X[0]], [z_min, z_max], color=combination_color)


def probe_crossing_point(probe1: Probe, probe2: Probe):
    M = np.array([-probe1.unit_vec, probe2.unit_vec]).T
    S = probe1.O - probe2.O
    T = np.linalg.inv(M) @ S
    X = probe1.O + T[0] * probe1.unit_vec
    return X


def compute_vertical_axis(probe1: Probe, probe2: Probe):
    z_min = max(probe1.z.min(), probe2.z.min())
    z_max = min(probe1.z.max(), probe2.z.max())
    n_pts = (
        ((probe1.z <= z_max) & (probe1.z >= z_min)).sum()
        + ((probe2.z <= z_max) & (probe2.z >= z_min)).sum()
    ) / 2
    z_interp = np.linspace(z_min, z_max, int(n_pts))
    return z_interp, z_min, z_max


if __name__ == "__main__":
    r = np.linspace(0, 5, 100)
    alpha1, alpha2 = -120, -70  # deg
    O1, O2 = np.array([1, 8]), np.array([-1, 7])
    probe1_pars = [r, alpha1, [0, O1]]
    probe2_pars = [r, alpha2, [0, O2]]
    #
    probe1 = Probe(*probe1_pars)
    probe2 = Probe(*probe2_pars)
    #
    sketch_probes(
        [probe1, probe2],
        combinations=[[0, 1]],
        probe_colors=[None, None],
        combination_colors=["k"],
    )
