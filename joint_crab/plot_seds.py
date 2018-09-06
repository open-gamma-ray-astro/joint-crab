"""Plot SEDs (Figure 2 in the paper)."""
import logging
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.table import Table
from gammapy.spectrum import CrabSpectrum
from .conf import config

log = logging.getLogger(__name__)


def main():
    fig, ax = plt.subplots()

    plot_meyer_2010()

    for name in config.all_datasets_plus_joint:
        plot_dataset(name)

    ax.legend(fontsize=config.plot.fontsize)
    ax.set_ylim([3e-13, 2e-10])
    ax.set_xlabel(config.plot.label_energy, size=config.plot.fontsize)
    ax.set_ylabel(config.plot.e2dnde_label, size=config.plot.fontsize)
    for axis in ["top", "bottom", "left", "right"]:
        ax.spines[axis].set_linewidth(1.6)
    opts = dict(axis="both", width=1.6, labelsize=config.plot.fontsize)
    ax.tick_params(which="major", length=7, **opts)
    ax.tick_params(which="minor", length=4, **opts)
    plt.tight_layout()

    filename = "results/figures/crab_sed_fit.png"
    log.info(f"Writing {filename}")
    fig.savefig(filename)

    filename = "results/figures/crab_sed_fit.pdf"
    log.info(f"Writing {filename}")
    fig.savefig(filename)


def plot_meyer_2010():
    model_meyer_ref = CrabSpectrum("meyer").model
    model_meyer_ref.plot(
        [10 * u.GeV, 100 * u.TeV],
        energy_power=2,
        flux_unit="erg-1 cm-2 s-1",
        ls=":",
        lw=2.2,
        color="#555555",
        label="Meyer et al. (2010)",
    )


def plot_dataset(name):
    dataset = config.datasets[name]
    ls = "-" if name == "joint" else "--"
    lw = 3 if name == "joint" else 2.2
    alpha = 0.38 if name == "joint" else 0.28

    table = Table.read(f"results/errorbands/sed_{name}.ecsv")
    mask = dataset.get_energy_mask(table["energy"].quantity)
    table = table[mask]

    x = table["energy"]
    y = table["e2dnde_fit"]
    y_lo = table["e2dnde_lo"]
    y_hi = table["e2dnde_hi"]

    plt.plot(x, y, ls=ls, lw=lw, color=dataset.color, label=dataset.label)

    plt.fill_between(x, y_lo, y_hi, color=dataset.color, alpha=alpha, label="")
