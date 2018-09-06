"""Plot figures from the paper."""
import logging
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from .utils import load_yaml
from .conf import config

log = logging.getLogger(__name__)


def plot_contours_stat():
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    for dataset in config.datasets.values():
        plot_contours_for_dataset(axes, dataset.name, dataset.color)

    add_legend(
        axes, [{"color": _.color, "label": _.label} for _ in config.datasets.values()]
    )

    plt.tight_layout()

    filename = "results/figures/contours.png"
    log.info(f"Writing {filename}")
    fig.savefig(filename)

    filename = "results/figures/contours.pdf"
    log.info(f"Writing {filename}")
    fig.savefig(filename)


def plot_contours_systematics():
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    for name in config.all_datasets:
        plot_contours_for_dataset(axes, name, color="gray")

    plot_contours_for_dataset(axes, "joint", color=config.datasets["joint"].color)
    plot_contours_for_dataset(axes, "joint_systematics", color="#002E63")

    add_legend(
        axes,
        [
            {"color": "crimson", "label": "joint fit \n stat. only"},
            {"color": "#002E63", "label": "joint fit \n stat. + syst."},
            {"color": "lightgray", "label": "instr. fit \n stat. only"},
        ],
    )

    plt.tight_layout()

    filename = "results/figures/contours_systematics.png"
    log.info(f"Writing {filename}")
    fig.savefig(filename)

    filename = "results/figures/contours_systematics.pdf"
    log.info(f"Writing {filename}")
    fig.savefig(filename)


def plot_contours_for_dataset(axes, name, color="lightgray"):
    fs = config.plot.fontsize_contours

    fit = load_yaml(f"results/fit/fit_{name}.yaml")["parameters"]
    fit = {_["name"]: _["value"] for _ in fit}
    contours = load_yaml(f"results/fit/contours_{name}.yaml")

    pars = {
        "phi": {
            "label": r"$\phi_0 \,/\,(10^{-11}\,{\rm TeV}^{-1} \, {\rm cm}^{-2} {\rm s}^{-1})$",
            "lim": [2.8, 5.2],
            "ticks": [3, 4, 5],
            "pos": np.array(1e11) * fit["amplitude"],
        },
        "gamma": {
            "label": r"$\Gamma$",
            "lim": [2.1, 2.9],
            "ticks": [2.2, 2.4, 2.6, 2.8],
            "pos": fit["alpha"],
        },
        "beta": {
            "label": r"$\beta$",
            "lim": [0.0, 0.8],
            "ticks": [0.2, 0.4, 0.6],
            "pos": fit["beta"],
        },
    }

    panels = [
        {
            "x": "phi",
            "y": "gamma",
            "cx": np.array(1e11) * contours["contour_amplitude_alpha"]["amplitude"],
            "cy": contours["contour_amplitude_alpha"]["alpha"],
        },
        {
            "x": "phi",
            "y": "beta",
            "cx": np.array(1e11) * contours["contour_amplitude_beta"]["amplitude"],
            "cy": contours["contour_amplitude_beta"]["beta"],
        },
        {
            "x": "gamma",
            "y": "beta",
            "cx": contours["contour_alpha_beta"]["alpha"],
            "cy": contours["contour_alpha_beta"]["beta"],
        },
    ]

    for p, ax in zip(panels, axes):
        x = pars[p["x"]]
        y = pars[p["y"]]
        plot_contour_line(ax, p["cx"], p["cy"], ls="-", lw=2.5, color=color)
        ax.plot(x["pos"], y["pos"], marker="X", markersize=7, color=color, lw=2.5)
        ax.set_xlabel(x["label"], size=fs)
        ax.set_ylabel(y["label"], size=fs)
        ax.set_xlim(x["lim"])
        ax.set_ylim(y["lim"])
        ax.set_xticks(x["ticks"])
        ax.set_yticks(y["ticks"])
        ax.tick_params("both", length=7, width=1.6, which="major", labelsize=fs)
        ax.tick_params("both", length=4, width=1.6, which="minor", labelsize=fs)
        for axis in ["top", "bottom", "left", "right"]:
            ax.spines[axis].set_linewidth(2.5)


def add_legend(axes, items):
    handles = [
        Line2D([], [], color=item["color"], ls="-", lw=2.5, label=item["label"])
        for item in items
    ]

    box = axes[2].get_position()
    axes[2].set_position([box.x0, box.y0, box.width * 0.97, box.height])
    axes[2].legend(
        handles=handles,
        loc="center left",
        fontsize=config.plot.fontsize,
        bbox_to_anchor=(1.0, 0.5),
    )


def plot_contour_line(ax, x, y, **kwargs):
    """Plot smooth curve from points.

    This doesn't seem to exist (or at least isn't easily accessible)
    in matplotlib. So we compute a spline using scipy following
    https://stackoverflow.com/questions/29837854
    """
    from scipy.interpolate import interp1d

    # ax.plot(x, y, "+")

    # Extend values so that spline doesn't have edge effects
    x = np.hstack([x, x, x])
    y = np.hstack([y, y, y])

    t = np.linspace(-1, 2, len(x), endpoint=False)
    tp = np.linspace(0, 1, 50)

    # kind="cubic" doesn't come out nicely.
    # Just use linear for now.
    xp = interp1d(t, x, kind="cubic")(tp)
    yp = interp1d(t, y, kind="cubic")(tp)

    ax.plot(xp, yp, **kwargs)
