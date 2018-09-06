"""Compute and plot counts spectrum (Fig 1 in the paper)."""
import logging
import matplotlib.pyplot as plt
from .conf import config


log = logging.getLogger(__name__)


def main():
    fig, ax = plt.subplots()

    for name in config.all_datasets:
        dataset = config.datasets[name]
        obs = dataset.get_SpectrumObservationList().stack()
        ex_obs = obs.excess_vector
        energy = ex_obs.energy.nodes.to("TeV")
        mask = dataset.get_energy_mask(energy)
        bins = ex_obs.energy.bins.to("TeV")
        weights = ex_obs.data.data.value

        ax.hist(
            energy.value[mask],
            bins=bins.value,
            weights=weights[mask],
            histtype="step",
            ls="-",
            lw=2.2,
            color=dataset.color,
            label=dataset.label,
        )

    ax.legend(fontsize=config.plot.fontsize)

    ax.set_xscale("log")
    ax.set_ylabel("Excess counts", size=config.plot.fontsize)
    ax.set_xlabel(r"$E'\,/\,\mathrm{TeV}$", size=config.plot.fontsize)

    for axis in ["top", "bottom", "left", "right"]:
        ax.spines[axis].set_linewidth(1.6)

    opts = dict(axis="both", width=1.6, labelsize=config.plot.fontsize)
    ax.tick_params(which="major", length=7, **opts)
    ax.tick_params(which="minor", length=4, **opts)

    plt.tight_layout()

    filename = "results/figures/counts_spectra.png"
    log.info(f"Writing {filename}")
    fig.savefig(filename)

    filename = "results/figures/counts_spectra.pdf"
    log.info(f"Writing {filename}")
    fig.savefig(filename)
