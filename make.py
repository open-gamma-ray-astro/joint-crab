#!/usr/bin/env python
"""Command line script to run the analysis.

All results and plot are produced like this:

$ ./make.py all
"""
import os
import logging
from distutils import dir_util
from pathlib import Path
import warnings
import yaml
import click
import matplotlib

# force a non-interactive backend to avoid windows popping up
matplotlib.use("agg")
click.disable_unicode_literals_warning = True

import joint_crab

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


@cli.command(name="all")
@click.pass_context
def cli_all(ctx):
    """Run all steps."""
    log.info("Executing task: all")

    ctx.invoke(cli_clean)
    ctx.invoke(cli_provenance)

    ctx.invoke(cli_maps)
    ctx.invoke(cli_extract_spectra)
    ctx.invoke(cli_fit_spec)
    ctx.invoke(cli_syst_err)
    ctx.invoke(cli_errorbands)

    ctx.invoke(cli_summary_data)
    ctx.invoke(cli_summary_results)

    ctx.invoke(cli_plot_counts)
    ctx.invoke(cli_plot_seds)
    ctx.invoke(cli_plot_errorbands)
    ctx.invoke(cli_plot_contours)

    # copy results from docker container to host
    if "DOCKER_INSIDE" in os.environ:
        dir_util.copy_tree("results", "host/joint-crab-results")
        log.info("Copying results to host.")


@cli.command("clean")
def cli_clean():
    """Clean out results folder."""
    log.info("Executing task: clean")

    dir_util.remove_tree("results")

    Path("results").mkdir()
    txt = "# Results\n\nThis folder contains the output of make.py"
    Path("results/README.md").write_text(txt)

    Path("results/maps").mkdir()
    Path("results/spectra").mkdir()
    Path("results/fit").mkdir()
    Path("results/figures").mkdir()
    Path("results/summary").mkdir()


@cli.command("provenance")
def cli_provenance():
    """Write `results/provenance.yaml`."""
    log.info("Executing task: provenance")
    data = joint_crab.provenance.get_provenance()

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
    log.info("Executing task: extract-spectra.")
    joint_crab.extract_ogip_spectra.extract_spectra(dataset)


@cli.command("maps")
def cli_maps():
    """Make and plot sky maps."""
    log.info("Executing task: maps")
    joint_crab.maps.main()


@cli.command("fit-spec")
@click.option(
    "--dataset",
    default="all",
    help="Which dataset to process? (default: all)",
    type=click.Choice(["fermi", "fact", "magic", "hess", "veritas", "joint", "all"]),
)
def cli_fit_spec(dataset):
    """Execute spectrum fit."""
    log.info("Executing task: fit-spec")
    joint_crab.fit_spec.main(dataset)


@cli.command("syst-err")
def cli_syst_err():
    """Fit that includes systematics."""
    log.info("Executing task: syst-err")
    joint_crab.systematics.main()


@cli.command("errorbands")
def cli_errorbands():
    """Compute flux error bands."""
    log.info("Executing task: errorbands")
    joint_crab.errorbands.main()


@cli.command("summary-data")
def cli_summary_data():
    """Write summary for data."""
    log.info("Executing task: summary-data")
    joint_crab.summary_data.make_summary_data()
    joint_crab.summary_data.make_summary_latex()


@cli.command("summary-results")
def cli_summary_results():
    """Write summary for results."""
    log.info("Executing task: summary-results")
    joint_crab.summary_results.make_summary_results()
    joint_crab.summary_results.make_summary_latex()


@cli.command("plot-counts")
def cli_plot_counts():
    """Plot counts spectra."""
    log.info("Executing task: plot-counts")
    joint_crab.plot_counts.main()


@cli.command("plot-seds")
def cli_plot_seds():
    """Plot SEDs."""
    log.info("Executing task: plot-seds")
    joint_crab.plot_seds.main()


@cli.command("plot-errorbands")
def cli_plot_errorbands():
    """Plot SED error bands."""
    log.info("Executing task: plot-errorbands")
    joint_crab.plot_errorbands.main()


@cli.command("plot-contours")
def cli_plot_contours():
    """Plot contours."""
    joint_crab.plot_contours.plot_contours_stat()
    joint_crab.plot_contours.plot_contours_systematics()


if __name__ == "__main__":
    cli()
