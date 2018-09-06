"""Print a results summary text file."""
import logging
from pathlib import Path
from . import utils

log = logging.getLogger(__name__)


def make_summary_results():
    lines = make_result_summary_lines()
    filename = "results/summary/results.md"
    log.info(f"Writing {filename}")
    Path(filename).write_text("\n".join(lines))


def make_result_summary_lines():
    yield "# Result summary\n"

    yield "## Spectral model fit results\n"
    yield "Log parabola; reference energy for amplitude parameter is at 1 TeV"
    yield ""
    for dataset in ["fermi", "magic", "hess", "fact", "veritas", "joint"]:
        yield f"### {dataset}"
        yield ""
        results = utils.load_yaml(f"results/fit/fit_{dataset}.yaml")
        for line in results_to_text(results):
            yield line
        yield ""


def results_to_text(results):
    p = results["parameters"]
    val, err = 1e12 * p[0]["value"], 1e12 * p[0]["error"]
    yield f"* amplitude: {val:.5f} +/- {err:.5f} 1e-12 cm-2 s-1 TeV-1"
    val, err = p[2]["value"], p[2]["error"]
    yield f"* alpha: {val:.5f} +/- {err:.5f}"
    val, err = p[3]["value"], p[3]["error"]
    yield f"* beta: {val:.5f} +/- {err:.5f}"


def make_summary_latex():
    """Make summary LaTeX table for the paper."""
    lines = []  # TODO: implement
    filename = "results/summary/results.tex"
    log.info(f"Writing {filename}")
    Path(filename).write_text("\n".join(lines))
