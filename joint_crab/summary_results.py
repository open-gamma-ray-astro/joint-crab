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
    path = "results/summary/results.tex"
    log.info(f"Writing {path}")

    table = open(path, "w")
    table.write(r"\begin{tabular}{lccc}")
    table.write("\n")
    table.write(r"\hline")
    table.write("\n")
    table.write(r"Dataset &  $\phi_0$ & $\Gamma$ & $\beta$ \\")
    table.write("\n")
    table.write("\hline")
    table.write("\n")

    for dataset in ["fermi", "magic", "veritas", "fact", "hess", "joint"]:
        results = utils.load_yaml(f"results/fit/fit_{dataset}.yaml")

        p = results["parameters"]
        phi_0_val, phi_0_err = 1e11 * p[0]["value"], 1e11 * p[0]["error"]
        gamma_val, gamma_err = p[2]["value"], p[2]["error"]
        beta_val, beta_err = p[3]["value"], p[3]["error"]

        row_name = rf"\{dataset} & " if dataset != "joint" else f"{dataset} & "

        phi_0 = f"{phi_0_val:.2f} $\pm$ {phi_0_err:.2f} & "
        gamma = f"{gamma_val:.2f} $\pm$ {gamma_err:.2f} & "
        beta = rf"{beta_val:.2f} $\pm$ {beta_err:.2f} \\"

        table.write(row_name + phi_0 + gamma + beta)
        table.write("\n")

    table.write("\hline")
    table.write("\n")
    table.write("\end{tabular}")
    table.close()
