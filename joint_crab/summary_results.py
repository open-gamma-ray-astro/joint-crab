"""Print a results summary text file."""
import logging
from pathlib import Path
from . import utils

log = logging.getLogger(__name__)


def make_summary_results():
    lines = make_result_summary_lines()
    filename = "results/summary_results.md"
    log.info(f"Writing {filename}")
    Path(filename).write_text("\n".join(lines))


def read_results():
    results = {}
    for dataset in ["fermi", "magic", "hess", "fact", "veritas", "joint"]:
        for tool in ["gammapy", "sherpa"]:
            path = f"results/fit/{tool}/{dataset}/fit_results_logparabola.yaml"
            data = utils.load_yaml(path)
            data["tool"] = tool
            data["dataset"] = dataset
            key = (tool, dataset)
            results[key] = data

    return results


def make_result_summary_lines():
    yield "# Result summary\n"

    yield "## Spectral model fit results\n"
    yield "Log parabola; reference energy for amplitude parameter is at 1 TeV"
    yield ""
    results = read_results()
    for result in results.values():
        for line in result_to_text(result):
            yield line


def result_to_text(r):
    p = r["parameters"]
    yield f'### {r["tool"]}, {r["dataset"]}'
    yield ""
    val, err = 1e12 * p[0]["value"], 1e12 * p[0]["error"]
    yield f"* amplitude: {val:.5f} +/- {err:.5f} 1e-12 cm-2 s-1 TeV"
    val, err = p[2]["value"], p[2]["error"]
    yield f"* alpha: {val:.5f} +/- {err:.5f}"
    val, err = p[3]["value"], p[3]["error"]
    yield f"* beta: {val:.5f} +/- {err:.5f}"
    yield ""
