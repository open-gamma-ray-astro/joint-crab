"""Plot SED error bands (Figure 3 in the paper)."""
import logging
import matplotlib.pyplot as plt
from astropy.table import Table
from .conf import config

log = logging.getLogger(__name__)


def main():
    for dataset in config.datasets.values():
        log.info(f"Plotting SED error band for {dataset.name}")
        plot_dataset(dataset)


def plot_dataset(dataset):
    fig, ax = plt.subplots()

    # Plot SED for sampled models
    table = Table.read(f"results/errorbands/samples_{dataset.name}.fits.gz")
    for idx in range(100):
        energy = table["energy"].quantity[idx]
        dnde = table["dnde"].quantity[idx]
        e2dnde = ((energy ** 2) * dnde).to("erg cm-2 s-1")
        mask = dataset.get_energy_mask(energy)
        plt.plot(
            energy.value[mask],
            e2dnde.value[mask],
            ls="-",
            lw=1,
            color="gray",
            alpha=0.3,
        )

    # Plot best-fit model and error band
    table = Table.read(f"results/errorbands/sed_{dataset.name}.ecsv")
    mask = dataset.get_energy_mask(table["energy"].quantity)
    table = table[mask]
    plt.plot(
        table["energy"],
        table["e2dnde_fit"],
        ls="--",
        lw=2,
        color="k",
        label="best fit model",
    )
    plt.plot(
        table["energy"],
        table["e2dnde_lo"],
        lw=2,
        color="k",
        label="68% containment\nmultivariate sampling",
    )
    plt.plot(
        table["energy"], table["e2dnde_hi"], lw=2, color="k", label="best fit model"
    )

    # Plot formatting
    ax.legend(loc="lower left", fontsize=config.plot.fontsize)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(config.plot.label_energy, size=config.plot.fontsize)
    ax.set_ylabel(config.plot.e2dnde_label, size=config.plot.fontsize)
    for axis in ["top", "bottom", "left", "right"]:
        ax.spines[axis].set_linewidth(1.6)
    opts = dict(axis="both", width=1.6, labelsize=config.plot.fontsize)
    ax.tick_params(which="major", length=7, **opts)
    ax.tick_params(which="minor", length=4, **opts)
    plt.tight_layout()

    filename = f"results/figures/errorband_sed_{dataset.name}.png"
    log.info(f"Writing {filename}")
    fig.savefig(filename)

    filename = f"results/figures/errorband_sed_{dataset.name}.pdf"
    log.info(f"Writing {filename}")
    fig.savefig(filename)
