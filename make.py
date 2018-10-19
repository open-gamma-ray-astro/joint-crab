#!/usr/bin/env python
import os
import logging
from distutils.dir_util import copy_tree
from pathlib import Path
import warnings
import yaml
import click
import matplotlib
from astropy.io import fits

# force a non-interactive backend to avoid windows popping up
matplotlib.use("agg")
click.disable_unicode_literals_warning = True

log = logging.getLogger(__name__)


@click.group()
@click.option(
    "--log-level",
    default="INFO",
    help="Logging verbosity level",
    type=click.Choice(["DEBUG", "INFO", "WARNING", "ERROR"]),
)
@click.option("--show-warnings", default=True, help="Show warnings?")
def cli(log_level, show_warnings):
    """Command line tool for this paper.

    To run all analyses, type `./make.py all`
    """
    logging.basicConfig(level=log_level)

    if not show_warnings:
        warnings.simplefilter("ignore")

    path = Path("results")
    path.mkdir(parents=True, exist_ok=True)


@cli.command("provenance")
def cli_provenance():
    """Write `results/provenance.yaml`."""
    from joint_crab.provenance import get_provenance

    data = get_provenance()

    filename = "results/provenance.yaml"
    log.info(f"Writing {filename}")

    with Path(filename).open("w") as fh:
        yaml.dump(data, fh, default_flow_style=False)


@cli.command("extract-spectra")
@click.option(
    "--dataset",
    default="all",
    help="Which dataset to process? (default: all)",
    type=click.Choice(["fermi", "fact", "magic", "hess", "veritas", "all"]),
)
def cli_extract_spectra(dataset):
    """Extract 1d spectra"""
    from joint_crab.extract_ogip_spectra import extract_spectra

    extract_spectra(dataset)


@cli.command("sky-images")
def cli_sky_images():
    """Compute and plot sky images."""
    from joint_crab.sky_image import main

    main()


@cli.command("fit-spec")
@click.option(
    "--dataset",
    default="all",
    help="Which dataset to process? (default: all)",
    type=click.Choice(["fermi", "fact", "magic", "hess", "veritas", "joint", "all"]),
)
@click.option(
    "--tool",
    default="all",
    help="Which tool to use? (default: all)",
    type=click.Choice(["sherpa", "gammapy", "all"]),
)
def cli_fit_spec(dataset, tool):
    """Execute spectrum fit"""
    from joint_crab.fit_spec import main

    main(dataset, tool)


@cli.command("plot-crab")
def cli_plot_crab():
    """Plot Crab SED"""
    from joint_crab.figures import plot_crab

    plot_crab()


@cli.command("plot-fit")
def cli_plot_fit():
    """Plot SEDs result of SpectrumFit and sherpa fit"""
    log.info("Executing plot-fit")
    from joint_crab.figures import (
        plot_fit_results,
        plot_sherpa_contours,
        plot_iminuit_contours,
        butterfly_syst,
        syst_contour,
    )

    plot_fit_results("gammapy")
    plot_fit_results("sherpa")
    plot_iminuit_contours()
    plot_sherpa_contours()
    butterfly_syst()
    syst_contour()


@cli.command("counts-histogram")
def cli_check_predicted_counts():
    """Plot counts histogram and check predicted counts"""
    from joint_crab.figures import counts_histogram

    counts_histogram(predicted=False)


@cli.command("stat-err")
@click.option(
    "--dataset",
    help="Which dataset to process?",
    default="all",
    type=click.Choice(["fermi", "fact", "magic", "hess", "veritas", "joint", "all"]),
)
def cli_butterfly_comparison(dataset):
    """estimate flux error band per each dataset"""
    log.info("Computing flux stat. error band")
    from joint_crab.figures import run_butterfly_stat

    run_butterfly_stat(dataset)


@cli.command("syst-err")
def cli_syst_err():
    """run the fit with the modified syst likelihood and plot the results agains the stat. only fit"""
    from joint_crab.errors import systematic, syst_errorband

    systematic()
    syst_errorband()


@cli.command("summary-data")
def cli_summary_data():
    """Write results/summary_data.md"""
    from joint_crab.summary_data import make_summary_data

    make_summary_data()


@cli.command("check-data")
def cli_check_data():
    """Check DL3 files for spec compliance.

    Write report in results/data_checks.
    """
    from joint_crab.test_dl3_fits import check_files
    path = Path("results/data_checks")
    path.mkdir(parents=True, exist_ok=True)

    for name in ["magic", "hess", "fact", "veritas"]:
        logfile = f"results/data_checks/check-data-{name}.txt"
        log.info(f"Writing {logfile}")

        # See https://docs.python.org/3.7/howto/logging-cookbook.html#logging-to-multiple-destinations
        logger = logging.getLogger("checkdl3")
        # See https://stackoverflow.com/a/2267567/498873
        logger.propagate = False

        fh = logging.FileHandler(logfile, mode="w")
        fh.setLevel(logging.INFO)
        logger.addHandler(fh)
        fh.setFormatter(logging.Formatter("%(levelname)8s %(message)s"))

        path = Path("data") / name
        filenames = sorted(path.glob("**/*.fits")) + sorted(path.glob("**/*.fits.gz"))

        # Note: we don't use the `errors` return value from `check_files` here because it's
        # too hard to work with (a nested dict of dicts of whatever)
        # It's not needed, the log contains the information.
        check_files(filenames)

        logger.removeHandler(fh)


@cli.command("check-data2")
@click.option(
    "--dataset",
    default="all",
    help="Which dataset to process? (default: all)",
    type=click.Choice(["fact", "magic", "hess", "veritas", "all"]),
)
def cli_check_data2(dataset):
    """Check DL3 with DataStore.check

    Write report in results/data_checks.
    """
    from gammapy.data import DataStore
    from joint_crab.utils import write_yaml
    path = Path("results/data_checks")
    path.mkdir(parents=True, exist_ok=True)

    if dataset == "all":
        names = ["magic", "hess", "fact", "veritas"]
    else:
        names = [dataset]

    for name in names:
        log.info(f"Running data checks for: {name}")
        data_store = DataStore.from_dir(f"data/{name}")

        results = data_store.check()
        results = (_ for _ in results if _["level"] not in {"debug", "info"})

        if name in ["magic", "fact", "veritas"]:
            results = (_ for _ in results if _["msg"] != "Loading psf failed")

        results = list(results)

        filename = f"results/data_checks/check-data-{name}2.yaml"
        write_yaml(results, filename)


@cli.command("check-data3")
@click.option(
    "--dataset",
    default="all",
    help="Which dataset to process? (default: all)",
    type=click.Choice(["fermi", "fact", "magic", "hess", "veritas", "all"]),
)
def cli_check_data3(dataset):
    """Check FITS headers of data files.

    Write report in results/data_checks.
    """
    path = Path("results/data_checks")
    path.mkdir(parents=True, exist_ok=True)

    if dataset == "all":
        names = ["fermi", "magic", "hess", "fact", "veritas"]
    else:
        names = [dataset]

    for name in names:
        log.info(f"Running data check 3 for: {name}")

        path = Path("data") / name
        paths = sorted(path.glob("**/*.fits")) + sorted(path.glob("**/*.fits.gz"))

        txt = ""
        for path in paths:
            hdu_list = fits.open(str(path))
            for hdu in hdu_list:
                txt += hdu.header.tostring(sep="\n", padding=False)
                txt += "\n\n**********\n\n"

        filename = f"results/data_checks/check-data-{name}3.txt"
        Path(filename).write_text(txt)


@cli.command("check-data4")
@click.option(
    "--dataset",
    default="all",
    help="Which dataset to process? (default: all)",
    type=click.Choice(["fermi", "fact", "magic", "hess", "veritas", "all"]),
)
def cli_check_data4(dataset):
    """Make standard peek plots.

    Write to results/data_checks/plots.
    """
    import matplotlib.pyplot as plt
    from gammapy.data import DataStore
    path = Path("results/data_checks/plots")
    path.mkdir(parents=True, exist_ok=True)

    if dataset == "all":
        names = ["magic", "hess", "fact", "veritas"]
    else:
        names = [dataset]

    for name in names:
        log.info(f"Running data check 4 for: {name}")

        path = Path("data") / name
        ds = DataStore.from_dir(path)
        for obs_id in ds.obs_table["OBS_ID"]:
            log.info(f"Processing {obs_id}")
            obs = ds.obs(obs_id)
            for hdu_type in ["events", "aeff", "edisp", "psf"]:
                try:
                    thing = getattr(obs, hdu_type)
                except Exception:

                    continue

                # if hdu_type == 'psf':
                #     plt.gcf().axes[2].set_ylim(0, 0.6)

                try:
                    thing.peek()
                except Exception:
                    log.exception("Oh no!")

                path = "results/data_checks/plots/"
                path += f"{name}_{hdu_type}_{obs_id:06d}.png"
                log.info(f"Writing {path}")
                plt.savefig(path)


@cli.command("summary-results")
def cli_summary_results():
    """Write results/summary_results.md"""
    from joint_crab.summary_results import make_summary_results

    make_summary_results()


@cli.command(name="all")
@click.pass_context
def cli_all(ctx):
    """Run all steps"""
    log.info("Run all steps ...")
    ctx.invoke(cli_provenance)
    ctx.invoke(cli_sky_images)
    ctx.invoke(cli_extract_spectra)
    ctx.invoke(cli_fit_spec)
    ctx.invoke(cli_plot_crab)
    ctx.invoke(cli_syst_err)
    ctx.invoke(cli_butterfly_comparison)
    ctx.invoke(cli_plot_fit)
    ctx.invoke(cli_check_predicted_counts)
    ctx.invoke(cli_check_data)
    ctx.invoke(cli_check_data2)
    ctx.invoke(cli_check_data3)
    ctx.invoke(cli_summary_data)
    ctx.invoke(cli_summary_results)


@cli.command(name="all-docker")
@click.pass_context
def make_all_docker(ctx):
    """Run all steps inside the docker container and copy results folder to host"""
    if "DOCKER_INSIDE" in os.environ:
        ctx.invoke(cli_all)
        copy_tree("results", "host/joint-crab-results")
        log.info("Copying results to host ...")
    else:
        log.info("The all-docker option is only allowed from a docker container shell")


if __name__ == "__main__":
    cli()
