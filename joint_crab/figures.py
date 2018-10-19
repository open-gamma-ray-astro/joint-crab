"""Plot figures from the paper."""
import logging
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from astropy import units as u
from astropy.table import Table
from gammapy.spectrum import CrabSpectrum, CountsPredictor
from gammapy.spectrum.models import LogParabola
import matplotlib.lines as mlines
from matplotlib import gridspec
from .models import Log10Parabola
from .utils import load_yaml
from .conf import config
from .errors import stat_errorband

log = logging.getLogger(__name__)

FONTSIZE = 15
FONTSIZE_CONTOURS = 18
E_UNIT_LABEL = r"$E\,/\,\mathrm{TeV}$"
SED_UNIT_LABEL = (
    r"$E^2 \cdot {\rm d}\phi/{\rm d}E\,/\,({\rm erg}\,{\rm cm}^{-2} {\rm s}^{-1})$"
)
COLORS = ["#21ABCD", "#FF9933", "#893F45", "#3EB489", "#002E63", "crimson"]


def plot_crab():
    """Plot Crab pulsar and nebula SED."""
    log.info("Executing plot_crab ...")

    fig, ax = plt.subplots()

    # Plot flux points
    for component in ["pulsar", "nebula"]:
        table = Table.read("data/other/crab_mwl.fits.gz")
        table = table[table["component"] == component]
        x = table["energy"].data
        y = table["energy_flux"].data
        yerr_lo = table["energy_flux_err_lo"].data
        yerr_hi = table["energy_flux_err_hi"].data
        ax.errorbar(x, y, yerr=(yerr_lo, yerr_hi), fmt="o", label=component)

    # Plot SED model
    energy = np.logspace(2, 8, 100) * u.MeV

    crab = CrabSpectrum(reference="meyer")

    flux = crab.model(energy)
    energy_flux = (energy ** 2 * flux).to("erg cm^-2 s^-1")
    ax.plot(energy.value, energy_flux.value, label="Meyer (2010) model", lw=3)

    ax.set_xlim((3e-1, 3e8))
    ax.set_ylim((3e-12, 3e-8))
    ax.set_xlabel("Energy (MeV)")
    ax.set_ylabel("E^2 dN/dE (erg cm^-2 s^-1)")
    fig.legend(loc="upper center", ncol=3)
    ax.grid()
    ax.loglog()

    path = Path("results/figures")
    path.mkdir(parents=True, exist_ok=True)
    filename = "results/figures/crab_mwl.png"
    log.info(f"Writing {filename}")
    fig.savefig(filename)


def plot_fit_results(tool):
    """plot the SEDs result of the gammapy / sherpa fit
    for comparison with the literature we choose only the Mayer spectrum

    Here we plot the butterfly as a result of the multivariate sampling
    with the 68% containment in flux
    """
    fig, ax = plt.subplots()

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

    # where to take the results, configurations for the individual butterflies
    instruments = ["fermi", "magic", "veritas", "fact", "hess", "joint"]
    labels = ["Fermi-LAT", "MAGIC", "VERITAS", "FACT", "H.E.S.S.", "joint fit"]
    lss = ["--", "--", "--", "--", "--", "-"]
    colors = COLORS
    # with one loop we realize all the butterfly plots
    for instrument, label, color, ls in zip(instruments, labels, colors, lss):

        path = (
            config.repo_path
            / f"results/fit/{tool}/{instrument}/fit_results_logparabola.yaml"
        )

        if not path.exists():
            log.warning(f"Missing: {path} . Skipping.")
            continue

        results = load_yaml(path)
        parameters = results["parameters"]

        model_lp = LogParabola.from_log10(
            amplitude=parameters[0]["value"] * u.Unit(parameters[0]["unit"]),
            reference=parameters[1]["value"] * u.Unit(parameters[1]["unit"]),
            alpha=parameters[2]["value"] * u.Unit(parameters[2]["unit"]),
            beta=parameters[3]["value"] * u.Unit(parameters[3]["unit"]),
        )

        # energy range for the plot
        dataset = config.get_dataset(instrument)
        energy_range = dataset.energy_range

        # just in case of the joint fit put a thicker line and a less transparent butterfly
        if instrument == "joint":
            model_lp.plot(
                energy_range,
                energy_power=2,
                flux_unit="erg-1 cm-2 s-1",
                ls=ls,
                lw=3,
                color=color,
                label=label,
            )
        else:
            model_lp.plot(
                energy_range,
                energy_power=2,
                flux_unit="erg-1 cm-2 s-1",
                ls=ls,
                lw=2.2,
                color=color,
                label=label,
            )

        # read the butterfly from the multivariate sampling results
        table_path = Path(
            f"{config.repo_path}/results/figures/stat_err/{instrument}_flux_errorband.dat"
        )
        log.info(f"reading butterfly values from {table_path}")
        t = Table.read(table_path, format="ascii.ecsv")
        energies = t["energies"].data * t["energies"].unit
        flux_lo = t["flux_lo"].data * t["flux_lo"].unit
        flux_hi = t["flux_hi"].data * t["flux_hi"].unit

        if instrument == "joint":
            alpha = 0.38
        else:
            alpha = 0.28

        plt.fill_between(
            energies.to("TeV"),
            (energies ** 2 * flux_lo).to("erg cm-2 s-1"),
            (energies ** 2 * flux_hi).to("erg cm-2 s-1"),
            color=color,
            alpha=alpha,
            label="",
        )

    ax.legend(fontsize=FONTSIZE)
    ax.set_ylim([1e-12, 2e-10])

    ax.set_xlabel(E_UNIT_LABEL, size=FONTSIZE)
    ax.set_ylabel(SED_UNIT_LABEL, size=FONTSIZE)
    # make axis thicker
    for axis in ["top", "bottom", "left", "right"]:
        ax.spines[axis].set_linewidth(1.6)
    ax.tick_params("both", length=7, width=1.6, which="major", labelsize=FONTSIZE)
    ax.tick_params("both", length=4, width=1.6, which="minor", labelsize=FONTSIZE)

    plt.tight_layout()

    filename = f"results/figures/crab_sed_{tool}_fit.png"
    filename_pdf = f"results/figures/crab_sed_{tool}_fit.pdf"
    log.info(f"Writing {filename}")
    fig.savefig(filename)
    fig.savefig(filename_pdf)


def plot_sherpa_contours():
    """plot the confidence contours obtained from the sherpa fit
    """
    log.info("plotting parameters contours obtained from sherpa")
    # where to take the results, configurations for the individual butterflies
    instruments = ["fermi", "magic", "hess", "fact", "veritas", "joint"]
    labels = ["Fermi-LAT", "MAGIC", "H.E.S.S.", "FACT", "VERITAS", "joint fit"]
    colors = ["#21ABCD", "#FF9933", "#5A4FCF", "#5CC184", "#702963", "crimson"]
    lss = ["--", "--", "--", "--", "--", "-"]

    fig, axarr = plt.subplots(1, 3, figsize=(18, 6))

    # with one loop we realize all the contour plots
    for instrument, label, color, ls in zip(instruments, labels, colors, lss):

        path = config.repo_path / f"results/fit/sherpa/{instrument}"

        contours_path = path / "fit_contours_logparabola.npy"
        results_path = path / "fit_results_logparabola.yaml"

        if not path.exists():
            log.warning(f"Missing: {path} . Skipping.")
            continue

        # load the contours and the results of the fit
        contours = np.load(contours_path).tolist()
        results = load_yaml(results_path)

        # define a 2 x 2 matrix to visualise the plot
        # we will delete one of the subplots and make something like a corner plot
        # useful variables for the plot
        ampl_range = contours["contour_ampl_c1"]["x0_range"]
        c1_range = contours["contour_ampl_c1"]["x1_range"]
        c2_range = contours["contour_ampl_c2"]["x1_range"]
        # actual values output of the fit
        # remember in sherpa notation: (amplitude->ampl, alpha->c1, beta->c2)
        ampl = results["parameters"][0]["value"]
        c1 = results["parameters"][2]["value"]
        c2 = results["parameters"][3]["value"]

        # axarr[0,0]
        extent = [ampl_range[0] * 1e9, ampl_range[1] * 1e9, c1_range[0], c1_range[1]]

        axarr[0].contour(
            contours["contour_ampl_c1"]["like_values"],
            contours["contour_ampl_c1"]["levels"],
            origin="lower",
            extent=extent,
            colors=color,
            linewidths=(2., 1.5, 1.3),
            linestyles=("-", "--", ":"),
        )

        # print actual value
        axarr[0].plot(ampl, c1, marker="X", markersize=7, color=color)
        axarr[0].set_xlabel(
            r"$f_0 / (\mathrm{TeV} \, \mathrm{cm}^{-2} \mathrm{s}^{-1})$"
        )
        axarr[0].set_ylabel(r"$\Gamma$")

        extent = [ampl_range[0] * 1e9, ampl_range[1] * 1e9, c2_range[0], c2_range[1]]

        axarr[1].contour(
            contours["contour_ampl_c2"]["like_values"],
            contours["contour_ampl_c2"]["levels"],
            origin="lower",
            extent=extent,
            colors=color,
            linewidths=(2., 1.5, 1.3),
            linestyles=("-", "--", ":"),
        )

        # print actual value
        axarr[1].plot(ampl, c2, marker="X", markersize=7, color=color)
        axarr[1].set_ylabel(r"$\beta$")
        axarr[1].set_xlabel(
            r"$f_0 / (\mathrm{TeV} \, \mathrm{cm}^{-2} \, \mathrm{s}^{-1})$"
        )

        extent = [c1_range[0], c1_range[1], c2_range[0], c2_range[1]]

        axarr[2].contour(
            contours["contour_c1_c2"]["like_values"],
            contours["contour_c1_c2"]["levels"],
            origin="lower",
            extent=extent,
            colors=color,
            linewidths=(2., 1.5, 1.3),
            linestyles=("-", "--", ":"),
        )

        # print actual value
        axarr[2].plot(c1, c2, marker="X", markersize=7, color=color)
        axarr[2].set_ylabel(r"$\beta$")
        axarr[2].set_xlabel(r"$\Gamma$")

    # axarr[0,1] is for the legend
    import matplotlib.lines as mlines

    sigma_1 = mlines.Line2D(
        [], [], color="k", marker="", ls="-", lw=2., label=r"1 $\sigma$ contour"
    )
    sigma_2 = mlines.Line2D(
        [], [], color="k", marker="", ls="--", lw=1.5, label=r"2 $\sigma$ contour"
    )
    sigma_3 = mlines.Line2D(
        [], [], color="k", marker="", ls=":", lw=1.3, label=r"3 $\sigma$ contour"
    )
    fermi = mlines.Line2D(
        [], [], color="#21ABCD", marker="", ls="-", lw=2., label="Fermi-LAT"
    )
    magic = mlines.Line2D(
        [], [], color="#FF9933", marker="", ls="-", lw=2., label="MAGIC"
    )
    hess = mlines.Line2D(
        [], [], color="#5A4FCF", marker="", ls="-", lw=2., label="H.E.S.S."
    )
    fact = mlines.Line2D(
        [], [], color="#5CC184", marker="", ls="-", lw=2., label="FACT"
    )
    veritas = mlines.Line2D(
        [], [], color="#702963", marker="", ls="-", lw=2., label="VERITAS"
    )
    joint = mlines.Line2D(
        [], [], color="crimson", marker="", ls="-", lw=2., label="joint fit"
    )
    axarr[2].legend(
        handles=[sigma_1, sigma_2, sigma_3, fermi, magic, hess, fact, veritas, joint],
        loc=3,
        fontsize=12,
    )
    # axarr[2].set_axis_off()

    plt.tight_layout()
    filename = "results/figures/sherpa_logparabola_contour.png"
    fig.savefig(filename)
    log.info(f"Writing {filename}")
    fig.savefig(filename)


def plot_iminuit_contours():
    """plot the confidence contours obtained from the sherpa fit
    """
    log.info("plotting parameters contours obtained from iminuit")
    # where to take the results, configurations for the individual butterflies
    instruments = ["fermi", "magic", "veritas", "fact", "hess", "joint"]
    labels = ["Fermi-LAT", "MAGIC", "VERITAS", "FACT", "H.E.S.S.", "joint fit"]
    colors = COLORS
    lss = ["--", "--", "--", "--", "--", "-"]

    fig, axarr = plt.subplots(1, 3, figsize=(16, 5))

    # with one loop we realize all the contour plots
    for instrument, label, color, ls in zip(instruments, labels, colors, lss):
        path = config.repo_path / f"results/fit/gammapy/{instrument}"

        contours_path = path / "fit_1.0_sigma_contours_logparabola.npy"
        results_path = path / "fit_results_logparabola.yaml"

        if not path.exists():
            log.warning(f"Missing: {path} . Skipping.")
            continue

        # load the contours and the results of the fit
        contours = np.load(contours_path).tolist()
        results = load_yaml(results_path)
        # true values to be plotted
        amplitude = float(results["parameters"][0]["value"])
        alpha = float(results["parameters"][2]["value"])
        beta = float(results["parameters"][3]["value"])

        # amplitude vs alpha
        amplitude_alpha = contours["contour_amplitude_alpha"]
        axarr[0].plot(
            amplitude_alpha["amplitude"] * 10,
            amplitude_alpha["alpha"],
            marker="",
            ls="-",
            lw=2.5,
            color=color,
        )
        # plot actual value
        axarr[0].plot(
            amplitude * 1e11, alpha, marker="X", markersize=7, color=color, lw=2.5
        )
        axarr[0].set_xlabel(
            r"$\phi_0 \,/\,(10^{-11}\,{\rm TeV} \, {\rm cm}^{-2} {\rm s}^{-1})$",
            size=FONTSIZE_CONTOURS,
        )
        axarr[0].set_ylabel(r"$\Gamma$", size=FONTSIZE_CONTOURS)
        # make axis thicker
        for axis in ["top", "bottom", "left", "right"]:
            axarr[0].spines[axis].set_linewidth(2.5)
        axarr[0].set_yticks([2.2, 2.4, 2.6, 2.8])
        axarr[0].set_ylim([2.1, 2.9])
        axarr[0].set_xticks([3, 4, 5])
        axarr[0].set_xlim([2.8, 5.2])
        axarr[0].tick_params(
            "both", length=7, width=1.6, which="major", labelsize=FONTSIZE_CONTOURS
        )
        axarr[0].tick_params(
            "both", length=4, width=1.6, which="minor", labelsize=FONTSIZE_CONTOURS
        )

        # amplitude vs beta
        amplitude_beta = contours["contour_amplitude_beta"]
        axarr[1].plot(
            amplitude_beta["amplitude"] * 10,
            # contour have a scale factor of 1e-10, parameters are in units of 1e-11
            amplitude_beta["beta"],
            marker="",
            ls="-",
            lw=2.5,
            color=color,
        )
        # plot actual value
        axarr[1].plot(amplitude * 1e11, beta, marker="X", markersize=7, color=color)
        axarr[1].set_xlabel(
            r"$\phi_0 \,/\,(10^{-11}\,{\rm TeV} \, {\rm cm}^{-2} {\rm s}^{-1})$",
            size=FONTSIZE_CONTOURS,
        )
        axarr[1].set_ylabel(r"$\beta$", size=FONTSIZE_CONTOURS)
        axarr[1].set_xticks([3, 4, 5])
        axarr[1].set_xlim([2.8, 5.2])
        axarr[1].set_yticks([0.2, 0.4, 0.6])
        axarr[1].set_ylim([0.0, 0.8])
        # make axis thicker
        for axis in ["top", "bottom", "left", "right"]:
            axarr[1].spines[axis].set_linewidth(2.5)
        axarr[1].tick_params(
            "both", length=7, width=1.6, which="major", labelsize=FONTSIZE_CONTOURS
        )
        axarr[1].tick_params(
            "both", length=4, width=1.6, which="minor", labelsize=FONTSIZE_CONTOURS
        )

        # alpha vs beta
        alpha_beta = contours["contour_alpha_beta"]
        axarr[2].plot(
            alpha_beta["alpha"],
            alpha_beta["beta"],
            marker="",
            ls="-",
            lw=2.5,
            color=color,
        )
        # plot actual value
        axarr[2].plot(alpha, beta, marker="X", markersize=7, color=color)
        axarr[2].set_xlabel(r"$\Gamma$", size=FONTSIZE_CONTOURS)
        axarr[2].set_ylabel(r"$\beta$", size=FONTSIZE_CONTOURS)
        axarr[2].set_xticks([2.2, 2.4, 2.6, 2.8])
        axarr[2].set_xlim([2.1, 2.9])
        axarr[2].set_yticks([0.2, 0.4, 0.6])
        axarr[2].set_ylim([0.0, 0.8])
        # make axis thicker
        for axis in ["top", "bottom", "left", "right"]:
            axarr[2].spines[axis].set_linewidth(2.5)
        axarr[2].tick_params(
            "both", length=7, width=1.6, which="major", labelsize=FONTSIZE_CONTOURS
        )
        axarr[2].tick_params(
            "both", length=4, width=1.6, which="minor", labelsize=FONTSIZE_CONTOURS
        )

    # legend
    import matplotlib.lines as mlines

    fermi = mlines.Line2D(
        [], [], color=COLORS[0], marker="", ls="-", lw=2.5, label="Fermi-LAT"
    )
    magic = mlines.Line2D(
        [], [], color=COLORS[1], marker="", ls="-", lw=2.5, label="MAGIC"
    )
    veritas = mlines.Line2D(
        [], [], color=COLORS[2], marker="", ls="-", lw=2.5, label="VERITAS"
    )
    fact = mlines.Line2D(
        [], [], color=COLORS[3], marker="", ls="-", lw=2.5, label="FACT"
    )
    hess = mlines.Line2D(
        [], [], color=COLORS[4], marker="", ls="-", lw=2.5, label="H.E.S.S."
    )
    joint = mlines.Line2D(
        [], [], color=COLORS[5], marker="", ls="-", lw=2.5, label="joint fit"
    )

    box = axarr[2].get_position()
    axarr[2].set_position([box.x0, box.y0, box.width * 0.97, box.height])
    # plot the legend on top of the central plot
    axarr[2].legend(
        handles=[fermi, magic, veritas, fact, hess, joint],
        loc="center left",
        fontsize=FONTSIZE_CONTOURS,
        bbox_to_anchor=(1., 0.5),
    )

    plt.tight_layout()
    filename = "results/figures/iminuit_logparabola_contour.png"
    filename_pdf = "results/figures/iminuit_logparabola_contour.pdf"
    fig.savefig(filename)
    log.info(f"Writing {filename}")
    fig.savefig(filename)
    fig.savefig(filename_pdf)


def run_butterfly_stat(which):
    """run the estimation of the butterfly"""
    if which in ["fermi", "all"]:
        butterfly_stat("fermi")
    if which in ["magic", "all"]:
        butterfly_stat("magic")
    if which in ["veritas", "all"]:
        butterfly_stat("veritas")
    if which in ["fact", "all"]:
        butterfly_stat("fact")
    if which in ["hess", "all"]:
        butterfly_stat("hess")
    if which in ["joint", "all"]:
        butterfly_stat("joint")


def butterfly_stat(which, tool="gammapy"):
    """plot an example figure showing how to correctly estimate a butterfly per each dataset
    We do it by default with the gammapy (iminuit) results
    """
    # the sampling in energy shall be consistent with the nbins of the spectrum extraction
    if which == "fermi":
        num_energy_points = 20
    else:
        num_energy_points = 60
    energies, flux_min, flux_max = stat_errorband(
        which=which, tool=tool, dim_sample=500, energy_points=num_energy_points, sigma=1
    )

    # plot best fit model
    fig, ax = plt.subplots()
    # best fit parameters for this fit
    result_file = f"results/fit/{tool}/{which}/fit_results_logparabola.yaml"
    results = load_yaml(result_file)

    parameters = results["parameters"]
    covariance = results["covariance"]

    # best fit parameters
    amplitude = parameters[0]["value"] * u.Unit(parameters[0]["unit"])
    reference = parameters[1]["value"] * u.Unit(parameters[1]["unit"])
    alpha = parameters[2]["value"] * u.Unit(parameters[2]["unit"])
    beta = parameters[3]["value"] * u.Unit(parameters[3]["unit"])

    # gammapy model with best-fit parameters
    model_lp = Log10Parabola(
        amplitude=amplitude, reference=reference, alpha=alpha, beta=beta
    )

    # set the covariance matrix from the output of the fit, this is needed to propagate the error
    model_lp.parameters.covariance = np.asarray(covariance)

    # plot some model representing the sampling
    # dictionary from the multivariate sampling
    path = Path(
        f"{config.repo_path}/results/debug/stat-err/{tool}/{which}/multivariate_sampling_fluxes.yaml"
    )
    sampled_dict = load_yaml(path)
    sampled_amplitude = sampled_dict["sampled_amplitude"]
    sampled_alpha = sampled_dict["sampled_alpha"]
    sampled_beta = sampled_dict["sampled_beta"]

    for (_ampl, _alpha, _beta) in zip(
        sampled_amplitude[:100], sampled_alpha[:100], sampled_beta[:100]
    ):
        _amplitude = _ampl * u.Unit(parameters[0]["unit"])
        _alpha = _alpha * u.Unit(parameters[2]["unit"])
        _beta = _beta * u.Unit(parameters[3]["unit"])

        _model_lp = Log10Parabola(
            amplitude=_amplitude, reference=reference, alpha=_alpha, beta=_beta
        )
        _model_lp.plot(
            energy_range=[energies[0], energies[-1]],
            energy_unit="TeV",
            flux_unit="cm-2 s-1 erg-1",
            energy_power=2,
            color="gray",
            lw=0.8,
            alpha=0.8,
            ax=ax,
        )

    # plot the 68 % containment correction
    ax.plot(energies, (energies ** 2 * flux_min).to("erg cm-2 s-1"), color="k", lw=2.5)
    ax.plot(
        energies,
        (energies ** 2 * flux_max).to("erg cm-2 s-1"),
        color="k",
        lw=2.5,
        label="68% containment" + "\n" + "multivariate sampling",
    )

    model_lp.plot(
        energy_range=[energies[0], energies[-1]],
        energy_unit="TeV",
        flux_unit="cm-2 s-1 erg-1",
        energy_power=2,
        color="k",
        lw=2.2,
        ls="--",
        label="best fit model",
        ax=ax,
    )

    # plt.title(f'{which} dataset')
    for axis in ["top", "bottom", "left", "right"]:
        ax.spines[axis].set_linewidth(2.)
    ax.tick_params("both", length=7, width=1.6, which="major", labelsize=FONTSIZE)
    ax.tick_params("both", length=4, width=1.6, which="minor", labelsize=FONTSIZE)
    ax.set_ylabel(SED_UNIT_LABEL, size=FONTSIZE)
    ax.set_xlabel(E_UNIT_LABEL, size=FONTSIZE)
    ax.set_xscale("log")
    ax.set_yscale("log")
    # ax.set_ylim([1e-12, 2e-10])
    ax.legend(fontsize=FONTSIZE)
    plt.tight_layout()

    figname = f"{config.repo_path}/results/figures/stat_err/{tool}_{which}_butterfly_comparison.png"
    figname_pdf = f"{config.repo_path}/results/figures/stat_err/{tool}_{which}_butterfly_comparison.pdf"
    logging.info(f"saving {figname}")
    fig.savefig(figname)
    fig.savefig(figname_pdf)

    log.info("making debug plots with histograms of the sampled parameters and fluxes")
    fig2, ax2 = plt.subplots(1, 3, sharey=True)
    nbins_params = int(len(sampled_amplitude) / 20)
    ax2[0].hist(sampled_amplitude, bins=nbins_params, color="lightgray")
    ax2[0].axvline(amplitude.value, lw=1.5, color="k")
    ax2[0].axvline(
        amplitude.value - np.sqrt(covariance[0][0]), lw=1.5, ls="--", color="k"
    )
    ax2[0].axvline(
        amplitude.value + np.sqrt(covariance[0][0]), lw=1.5, ls="--", color="k"
    )
    ax2[0].set_xlabel(
        r"$f_0 (\mathrm{TeV}^{-1} \, \mathrm{cm}^{-2} \, \mathrm{s}^{-1})$", labelpad=9
    )

    ax2[1].hist(sampled_alpha, bins=int(len(sampled_alpha) / 10), color="lightgray")
    ax2[1].axvline(alpha.value, lw=1.5, color="k")
    ax2[1].axvline(alpha.value - np.sqrt(covariance[2][2]), lw=1.5, ls="--", color="k")
    ax2[1].axvline(alpha.value + np.sqrt(covariance[2][2]), lw=1.5, ls="--", color="k")
    ax2[1].set_xlabel(r"$\Gamma$")

    ax2[2].hist(
        sampled_beta,
        bins=int(len(sampled_beta) / 10),
        color="lightgray",
        label="sampled",
    )
    ax2[2].axvline(beta.value, lw=1.5, color="k", label="best fit")
    ax2[2].axvline(
        beta.value - np.sqrt(covariance[3][3]),
        lw=1.5,
        ls="--",
        color="k",
        label=r"best fit $\pm$ 1 $\sigma$",
    )
    ax2[2].axvline(beta.value + np.sqrt(covariance[3][3]), lw=1.5, ls="--", color="k")
    ax2[2].legend()
    ax2[2].set_xlabel(r"$\beta$")

    figname2 = f"{config.repo_path}/results/figures/stat_err/{tool}_{which}_sampled_parameters.png"
    logging.info(f"saving {figname2}")
    fig2.savefig(figname2)

    # sampled fluxes
    fig3, ax3 = plt.subplots(1, 3, sharey=True, figsize=(12, 8))
    for i, ene in enumerate(["emin", "emid", "emax"]):
        flux_dict = sampled_dict["sampled_fluxes"][ene]
        energy = sampled_dict["sampled_fluxes"][ene]["value"] * u.Unit(
            sampled_dict["sampled_fluxes"]["energy_unit"]
        )
        log.info(f"building histogram of sampled fluxes for {energy}")
        nbins_fluxes = int(len(flux_dict["fluxes"]) / 20)
        bins_fluxes = np.logspace(
            np.log10(np.min(flux_dict["fluxes"])),
            np.log10(np.max(flux_dict["fluxes"])),
            nbins_fluxes,
        )
        ax3[i].hist(
            flux_dict["fluxes"],
            bins=bins_fluxes,
            color="lightgray",
            label="sampled fluxes",
        )
        ax3[i].axvline(
            np.mean(flux_dict["fluxes"]),
            lw=1.5,
            ls="-",
            color="k",
            label="mean sampled",
        )
        ax3[i].axvline(
            flux_dict["flux_quantiles"][0],
            lw=1.5,
            ls="--",
            color="k",
            label="68% containment sampled",
        )
        ax3[i].axvline(flux_dict["flux_quantiles"][1], lw=1.5, ls="--", color="k")
        ax3[i].set_xlabel(
            r"$F (\mathrm{TeV}^{-1} \, \mathrm{cm}^{-2} \, \mathrm{s}^{-1})$",
            labelpad=9,
        )
        # evaluate normal model and error
        fit_value = model_lp(energy).to("TeV-1 cm-2 s-1").value
        fit_value_err = model_lp.evaluate_error(energy).to("TeV-1 cm-2 s-1").value
        ax3[i].axvline(
            fit_value, lw=2, ls="-", color="crimson", label="best fit result"
        )
        ax3[i].axvline(
            fit_value + fit_value_err[1],
            lw=1.5,
            ls="--",
            color="crimson",
            label="error propagation",
        )
        ax3[i].axvline(fit_value - fit_value_err[1], lw=1.5, ls="--", color="crimson")
        ax3[i].set_title("E = {:.2f}".format(energy))
        ax3[i].legend()
        ax3[i].set_xscale("log")

    plt.tight_layout()
    figname3 = (
        f"{config.repo_path}/results/figures/stat_err/{tool}_{which}_sampled_fluxes.png"
    )
    logging.info(f"saving {figname3}")
    fig3.savefig(figname3)


def counts_histogram(predicted=False):
    """function to plot the excesses per dataset and compare them with the predicted counts
    if predicted == True will shwo the predicted counts from the results of the fit (for debug purpose)
    """
    log.info("loading the results from the joint fit to predict the counts")
    results = load_yaml(
        f"{config.repo_path}/results/fit/gammapy/joint/fit_results_logparabola.yaml"
    )
    parameters = results["parameters"]

    model_lp = LogParabola.from_log10(
        amplitude=parameters[0]["value"] * u.Unit(parameters[0]["unit"]),
        reference=parameters[1]["value"] * u.Unit(parameters[1]["unit"]),
        alpha=parameters[2]["value"] * u.Unit(parameters[2]["unit"]),
        beta=parameters[3]["value"] * u.Unit(parameters[3]["unit"]),
    )

    # defining the figure
    dict_color = {
        "fermi": COLORS[0],
        "magic": COLORS[1],
        "veritas": COLORS[2],
        "fact": COLORS[3],
        "hess": COLORS[4],
    }
    fig, ax = plt.subplots()

    for which in config.all_datasets:
        log.info(f"predicting counts for {which} dataset")
        dataset = config.get_dataset(which)
        obs = dataset.get_SpectrumObservationList().stack()
        cts_pred = CountsPredictor(
            model=model_lp, aeff=obs.aeff, edisp=obs.edisp, livetime=obs.livetime
        )
        cts_pred.run()

        e_max = dataset.energy_range[1].to("TeV").value
        e_min = dataset.energy_range[0].to("TeV").value

        kwargs_mdl = dict(ls=":", range=(e_min, e_max), lw=2.2, color=dict_color[which])
        kwargs_data = dict(
            ls="-", range=(e_min, e_max), lw=2.2, color=dict_color[which]
        )

        # CountsSpectrum with observed and predicted excesses
        ex_pred = cts_pred.npred
        ex_obs = obs.excess_vector
        # if it is an IACT rebin the counts before plotting
        if which != "fermi":
            ex_pred = ex_pred.rebin(2)
            ex_obs = ex_obs.rebin(2)

        if predicted:  # if you want to display the predicted counts
            ex_pred.plot_hist(ax, **kwargs_mdl)
        ex_obs.plot_hist(ax, **kwargs_data)

    # custom legend
    legend_observed = mlines.Line2D(
        [], [], color="gray", marker="", ls="-", lw=2, label="observed"
    )
    legend_expected = mlines.Line2D(
        [], [], color="gray", marker="", ls=":", lw=2, label="expected"
    )
    legend_fermi = mlines.Line2D(
        [], [], color=COLORS[0], marker="", ls="-", lw=2, label="Fermi-LAT"
    )
    legend_magic = mlines.Line2D(
        [], [], color=COLORS[1], marker="", ls="-", lw=2, label="MAGIC"
    )
    legend_veritas = mlines.Line2D(
        [], [], color=COLORS[2], marker="", ls="-", lw=2, label="VERITAS"
    )
    legend_fact = mlines.Line2D(
        [], [], color=COLORS[3], marker="", ls="-", lw=2, label="FACT"
    )
    legend_hess = mlines.Line2D(
        [], [], color=COLORS[4], marker="", ls="-", lw=2, label="H.E.S.S."
    )
    legend_handles = [
        legend_fermi,
        legend_magic,
        legend_veritas,
        legend_fact,
        legend_hess,
    ]
    if predicted:  # if you want to display the predicted counts
        legend_handles = [legend_observed, legend_expected] + legend_handles

    ax.legend(handles=legend_handles, fontsize=FONTSIZE)

    ax.set_xscale("log")
    ax.set_ylabel("Excess counts", size=FONTSIZE)
    ax.set_xlabel(E_UNIT_LABEL, size=FONTSIZE)

    # make axis thicker
    for axis in ["top", "bottom", "left", "right"]:
        ax.spines[axis].set_linewidth(1.6)
    ax.tick_params("both", length=7, width=1.6, which="major", labelsize=FONTSIZE)
    ax.tick_params("both", length=4, width=1.6, which="minor", labelsize=FONTSIZE)

    plt.tight_layout()

    filename = f"{config.repo_path}/results/figures/counts_spectra.png"
    filename_pdf = f"{config.repo_path}/results/figures/counts_spectra.pdf"
    log.info(f"saving figure in {filename}")
    fig.savefig(filename)
    fig.savefig(filename_pdf)


def butterfly_syst():
    """representation of the joint fit result w/ systematics uncertainty"""
    log.info("plotting syst. + stat. butterfly")
    # reading the joint stat butterfly values
    stat_table_path = Path(
        f"{config.repo_path}/results/figures/stat_err/joint_flux_errorband.dat"
    )
    log.info(f"reading butterfly values from {stat_table_path}")
    stat_table = Table.read(stat_table_path, format="ascii.ecsv")
    stat_energies = stat_table["energies"].data * stat_table["energies"].unit
    stat_flux_min = stat_table["flux_lo"].data * stat_table["flux_lo"].unit
    stat_flux_max = stat_table["flux_hi"].data * stat_table["flux_hi"].unit
    # load result of stat fit
    stat_result_file = f"results/fit/gammapy/joint/fit_results_logparabola.yaml"
    stat_results = load_yaml(stat_result_file)

    stat_parameters = stat_results["parameters"]

    # best fit parameters
    stat_amplitude = stat_parameters[0]["value"] * u.Unit(stat_parameters[0]["unit"])
    stat_reference = stat_parameters[1]["value"] * u.Unit(stat_parameters[1]["unit"])
    stat_alpha = stat_parameters[2]["value"] * u.Unit(stat_parameters[2]["unit"])
    stat_beta = stat_parameters[3]["value"] * u.Unit(stat_parameters[3]["unit"])

    # reading the joint stat butterfly values
    syst_table_path = Path(
        f"{config.repo_path}/results/figures/syst_err/joint_flux_errorband.dat"
    )
    log.info(f"reading butterfly values from {syst_table_path}")
    syst_table = Table.read(syst_table_path, format="ascii.ecsv")
    syst_energies = syst_table["energies"].data * syst_table["energies"].unit
    syst_flux_min = syst_table["flux_lo"].data * syst_table["flux_lo"].unit
    syst_flux_max = syst_table["flux_hi"].data * syst_table["flux_hi"].unit

    # load result of syst fit
    syst_result_file = (
        f"results/fit/gammapy/joint/fit_results_logparabola_energy_scale.yaml"
    )
    syst_results = load_yaml(syst_result_file)

    syst_parameters = syst_results["parameters"]

    # best fit parameters
    syst_amplitude = syst_parameters[0]["value"] * u.Unit(syst_parameters[0]["unit"])
    syst_reference = syst_parameters[1]["value"] * u.Unit(syst_parameters[1]["unit"])
    syst_alpha = syst_parameters[2]["value"] * u.Unit(syst_parameters[2]["unit"])
    syst_beta = syst_parameters[3]["value"] * u.Unit(syst_parameters[3]["unit"])

    fig = plt.figure()
    gs = gridspec.GridSpec(2, 1, height_ratios=[2.5, 1], wspace=0.)
    ax0 = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1], sharex=ax0)

    ax0.fill_between(
        syst_energies.to("TeV"),
        (syst_energies ** 2 * syst_flux_min).to("erg cm-2 s-1"),
        (syst_energies ** 2 * syst_flux_max).to("erg cm-2 s-1"),
        color="#002E63",
        alpha=0.4,
        label="stat. + syst.",
    )

    ax0.fill_between(
        stat_energies.to("TeV"),
        (stat_energies ** 2 * stat_flux_min).to("erg cm-2 s-1"),
        (stat_energies ** 2 * stat_flux_max).to("erg cm-2 s-1"),
        color="crimson",
        alpha=0.4,
        label="stat. only",
    )

    stat_flux = Log10Parabola.evaluate(
        stat_energies, stat_amplitude, stat_reference, stat_alpha, stat_beta
    )

    ax0.plot(
        stat_energies.to("TeV"),
        (stat_energies ** 2 * stat_flux).to("erg cm-2 s-1"),
        color="crimson",
        lw=2.2,
        ls="-",
    )

    syst_flux = Log10Parabola.evaluate(
        syst_energies, syst_amplitude, syst_reference, syst_alpha, syst_beta
    )

    ax0.plot(
        syst_energies.to("TeV"),
        (syst_energies ** 2 * syst_flux).to("erg cm-2 s-1"),
        color="#002E63",
        lw=2.2,
        ls="-",
    )

    ax0.set_xscale("log")
    ax0.set_yscale("log")
    ax0.set_ylabel(SED_UNIT_LABEL, size=FONTSIZE)
    ax0.legend(fontsize=FONTSIZE, loc=3)
    # make axis thicker
    for axis in ["top", "bottom", "left", "right"]:
        ax0.spines[axis].set_linewidth(1.6)
    ax0.tick_params("both", length=7, width=1.6, which="major", labelsize=FONTSIZE)
    ax0.tick_params("both", length=4, width=1.6, which="minor", labelsize=FONTSIZE)

    # plot the flux ratios
    ax1.fill_between(
        syst_energies.to("TeV"),
        (syst_flux_min - stat_flux) / stat_flux,
        (syst_flux_max - stat_flux) / stat_flux,
        color="#002E63",
        alpha=0.4,
        label="stat. + syst.",
    )

    ax1.fill_between(
        stat_energies.to("TeV"),
        (stat_flux_min - stat_flux) / stat_flux,
        (stat_flux_max - stat_flux) / stat_flux,
        color="crimson",
        alpha=0.4,
        label="stat. only",
    )

    ax1.plot(
        stat_energies.to("TeV"),
        (stat_flux - stat_flux) / stat_flux,
        color="crimson",
        lw=2.2,
        ls="-",
    )

    ax1.plot(
        syst_energies.to("TeV"),
        (syst_flux - stat_flux) / stat_flux,
        color="#002E63",
        lw=2.2,
        ls="-",
    )

    ax1.set_xscale("log")
    ax1.set_ylim([-0.35, 0.35])
    ax1.set_xlabel(E_UNIT_LABEL, size=FONTSIZE)
    ax1.set_ylabel("fractional diff. \n to stat.", size=FONTSIZE)
    # make axis thicker
    for axis in ["top", "bottom", "left", "right"]:
        ax1.spines[axis].set_linewidth(1.6)
    ax1.tick_params("both", length=7, width=1.6, which="major", labelsize=FONTSIZE)
    ax1.tick_params("both", length=4, width=1.6, which="minor", labelsize=FONTSIZE)

    plt.tight_layout()

    filename = "results/figures/crab_sed_joint_fit_syst.png"
    filename_pdf = "results/figures/crab_sed_joint_fit_syst.pdf"
    log.info(f"Writing {filename}")
    fig.savefig(filename)
    fig.savefig(filename_pdf)


def syst_contour():
    """plot a comparison of the contours obtained with the wstat statistics
    and with the stat + syst likelihood"""

    fig, axarr = plt.subplots(1, 3, figsize=(16, 5))

    # first plot in light gray the stat contours
    instruments = ["fermi", "magic", "veritas", "fact", "hess", "joint"]

    # with one loop we realize all the contour plots
    color = "lightgray"
    for instrument in instruments:
        path = f"{config.repo_path}/results/fit/gammapy/{instrument}"

        contours_path = f"{path}/fit_1.0_sigma_contours_logparabola.npy"
        results_path = f"{path}/fit_results_logparabola.yaml"

        # load the contours and the results of the fit
        contours = np.load(contours_path).tolist()
        results = load_yaml(results_path)
        # true values to be plotted
        amplitude = float(results["parameters"][0]["value"])
        alpha = float(results["parameters"][2]["value"])
        beta = float(results["parameters"][3]["value"])

        # amplitude vs alpha
        amplitude_alpha = contours["contour_amplitude_alpha"]
        axarr[0].plot(
            amplitude_alpha["amplitude"] * 10,
            amplitude_alpha["alpha"],
            marker="",
            ls="-",
            lw=2.5,
            color=color,
        )
        # plot actual value
        axarr[0].plot(
            amplitude * 1e11, alpha, marker="X", markersize=7, color=color, lw=2.5
        )

        # amplitude vs beta
        amplitude_beta = contours["contour_amplitude_beta"]
        axarr[1].plot(
            amplitude_beta["amplitude"] * 10,
            # contour have a scale factor of 1e-10, parameters are in units of 1e-11
            amplitude_beta["beta"],
            marker="",
            ls="-",
            lw=2.5,
            color=color,
        )
        # plot actual value
        axarr[1].plot(amplitude * 1e11, beta, marker="X", markersize=7, color=color)

        # alpha vs beta
        alpha_beta = contours["contour_alpha_beta"]
        axarr[2].plot(
            alpha_beta["alpha"],
            alpha_beta["beta"],
            marker="",
            ls="-",
            lw=2.5,
            color=color,
        )
        # plot actual value
        axarr[2].plot(alpha, beta, marker="X", markersize=7, color=color)

    # plotting stat vs syst result on top of instrument-wise contours in gray
    path = config.repo_path / f"results/fit/gammapy/joint"

    stat_contours_path = f"{path}/fit_1.0_sigma_contours_logparabola.npy"
    stat_results_path = f"{path}/fit_results_logparabola.yaml"

    syst_contours_path = f"{path}/fit_1.0_sigma_contours_logparabola_energy_scale.npy"
    syst_results_path = f"{path}/fit_results_logparabola_energy_scale.yaml"

    for (contours_path, results_path, color, label) in zip(
        [stat_contours_path, syst_contours_path],
        [stat_results_path, syst_results_path],
        ["crimson", "#002E63"],
        ["stat. only", "stat. + syst."],
    ):
        # load the contours and the results of the fit
        contours = np.load(contours_path).tolist()
        results = load_yaml(results_path)
        # true values to be plotted
        amplitude = float(results["parameters"][0]["value"])
        alpha = float(results["parameters"][2]["value"])
        beta = float(results["parameters"][3]["value"])

        # amplitude vs alpha
        amplitude_alpha = contours["contour_amplitude_alpha"]
        if label == "stat. only":
            axarr[0].plot(
                amplitude_alpha["amplitude"] * 10,
                # gammapy contour have a scale factor of 1e-10, parameters are in units of 1e-11
                amplitude_alpha["alpha"],
                marker="",
                ls="-",
                lw=2.5,
                color=color,
            )
        else:
            axarr[0].plot(
                amplitude_alpha["amplitude"] * 1e11,
                amplitude_alpha["alpha"],
                marker="",
                ls="-",
                lw=2.5,
                color=color,
            )

        # plot actual value
        axarr[0].plot(
            amplitude * 1e11, alpha, marker="X", markersize=7, color=color, lw=2.5
        )

        # amplitude vs beta
        amplitude_beta = contours["contour_amplitude_beta"]
        if label == "stat. only":
            axarr[1].plot(
                amplitude_beta["amplitude"] * 10,
                # gammapy contour have a scale factor of 1e-10, parameters are in units of 1e-11
                amplitude_beta["beta"],
                marker="",
                ls="-",
                lw=2.5,
                color=color,
            )
        else:
            axarr[1].plot(
                amplitude_beta["amplitude"] * 1e11,
                amplitude_beta["beta"],
                marker="",
                ls="-",
                lw=2.5,
                color=color,
            )

        # plot actual value
        axarr[1].plot(amplitude * 1e11, beta, marker="X", markersize=7, color=color)

        # alpha vs beta
        alpha_beta = contours["contour_alpha_beta"]
        axarr[2].plot(
            alpha_beta["alpha"],
            alpha_beta["beta"],
            marker="",
            ls="-",
            lw=2.5,
            color=color,
        )
        # plot actual value
        axarr[2].plot(alpha, beta, marker="X", markersize=7, color=color)

    # legend
    import matplotlib.lines as mlines

    stat_label = mlines.Line2D(
        [],
        [],
        color="crimson",
        marker="",
        ls="-",
        lw=2.5,
        label="joint fit \n stat. only",
    )
    syst_label = mlines.Line2D(
        [],
        [],
        color="#002E63",
        marker="",
        ls="-",
        lw=2.5,
        label="joint fit \n stat. + syst.",
    )
    single_label = mlines.Line2D(
        [],
        [],
        color="lightgray",
        marker="",
        ls="-",
        lw=2.5,
        label="instr. fit \n stat. only",
    )

    axarr[0].set_xlabel(
        r"$\phi_0 \,/\,(10^{-11}\,{\rm TeV} \, {\rm cm}^{-2} {\rm s}^{-1})$",
        size=FONTSIZE_CONTOURS,
    )
    axarr[0].set_ylabel(r"$\Gamma$", size=FONTSIZE_CONTOURS)
    # make axis thicker
    axarr[0].set_yticks([2.2, 2.4, 2.6, 2.8])
    axarr[0].set_ylim([2.1, 2.9])
    axarr[0].set_xticks([3, 4, 5])
    axarr[0].set_xlim([2.8, 5.2])
    for axis in ["top", "bottom", "left", "right"]:
        axarr[0].spines[axis].set_linewidth(2.5)
    axarr[0].tick_params(
        "both", length=7, width=1.6, which="major", labelsize=FONTSIZE_CONTOURS
    )
    axarr[0].tick_params(
        "both", length=4, width=1.6, which="minor", labelsize=FONTSIZE_CONTOURS
    )

    axarr[1].set_xlabel(
        r"$\phi_0 \,/\,(10^{-11}\,{\rm TeV} \, {\rm cm}^{-2} {\rm s}^{-1})$",
        size=FONTSIZE_CONTOURS,
    )
    axarr[1].set_ylabel(r"$\beta$", size=FONTSIZE_CONTOURS)
    # make axis thicker
    axarr[1].set_xticks([3, 4, 5])
    axarr[1].set_xlim([2.8, 5.2])
    axarr[1].set_yticks([0.2, 0.4, 0.6])
    axarr[1].set_ylim([0.0, 0.8])
    for axis in ["top", "bottom", "left", "right"]:
        axarr[1].spines[axis].set_linewidth(2.5)
    axarr[1].tick_params(
        "both", length=7, width=1.6, which="major", labelsize=FONTSIZE_CONTOURS
    )
    axarr[1].tick_params(
        "both", length=4, width=1.6, which="minor", labelsize=FONTSIZE_CONTOURS
    )

    axarr[2].set_xlabel(r"$\Gamma$", size=FONTSIZE_CONTOURS)
    axarr[2].set_ylabel(r"$\beta$", size=FONTSIZE_CONTOURS)
    axarr[2].set_xticks([2.2, 2.4, 2.6, 2.8])
    axarr[2].set_xlim([2.1, 2.9])
    axarr[2].set_yticks([0.2, 0.4, 0.6])
    axarr[2].set_ylim([0.0, 0.8])
    # make axis thicker
    for axis in ["top", "bottom", "left", "right"]:
        axarr[2].spines[axis].set_linewidth(2.5)
    axarr[2].tick_params(
        "both", length=7, width=1.6, which="major", labelsize=FONTSIZE_CONTOURS
    )
    axarr[2].tick_params(
        "both", length=4, width=1.6, which="minor", labelsize=FONTSIZE_CONTOURS
    )

    box = axarr[2].get_position()
    axarr[2].set_position([box.x0, box.y0, box.width * 0.97, box.height])
    # plot the legend on top of the central plot
    axarr[2].legend(
        handles=[stat_label, syst_label, single_label],
        loc="center left",
        fontsize=FONTSIZE_CONTOURS,
        bbox_to_anchor=(1., 0.5),
    )

    plt.tight_layout()
    filename = "results/figures/iminuit_logparabola_energy_scale_contour.png"
    filename_pdf = "results/figures/iminuit_logparabola_energy_scale_contour.pdf"
    fig.savefig(filename)
    log.info(f"Writing {filename}")
    fig.savefig(filename)
    fig.savefig(filename_pdf)
