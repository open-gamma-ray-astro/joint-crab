"""Run spectral fits with Gammapy and sherpa."""
import logging
from gammapy.spectrum import SpectrumFit
from . import utils
from .conf import config
from .fit_models import Log10Parabola

log = logging.getLogger(__name__)


def main(dataset):
    if dataset in {"fermi", "all"}:
        run_analysis("fermi")
    if dataset in {"magic", "all"}:
        run_analysis("magic")
    if dataset in {"hess", "all"}:
        run_analysis("hess")
    if dataset in {"fact", "all"}:
        run_analysis("fact")
    if dataset in {"veritas", "all"}:
        run_analysis("veritas")
    if dataset in {"joint", "all"}:
        run_analysis("joint")


def run_analysis(which):
    log.info(f"Fitting dataset: {which}")

    dataset = config.datasets[which]
    obs_list = dataset.get_SpectrumObservationList()
    fit_range = dataset.energy_range
    log.info(f"obs_list: {obs_list}")

    model = Log10Parabola(
        amplitude="3.8e-11 cm-2 s-1 TeV-1", reference="1 TeV", alpha=2.47, beta=0.24
    )

    fit = SpectrumFit(obs_list=obs_list, model=model, fit_range=fit_range)

    log.info("Running fit")
    fit.run()

    fit_results = make_results_dict(fit)
    utils.write_yaml(fit_results, f"results/fit/fit_{which}.yaml")

    contour_results = compute_contours(fit)
    utils.write_yaml(contour_results, f"results/fit/contours_{which}.yaml")


# TODO: Use parameters.to_dict instead!
def make_results_dict(fit):
    pars = fit.result[0].model.parameters
    results = {}
    results["parameters"] = []
    for par in pars:
        info = par.to_dict()
        info["error"] = float(pars.error(par))
        results["parameters"].append(info)
    results["covariance"] = pars.covariance.tolist()
    results["statname"] = fit.result[0].statname
    results["statval"] = float(fit.result[0].statval)
    results["fit_range"] = {
        "min": float(fit.result[0].fit_range[0].to("keV").value),
        "max": float(fit.result[0].fit_range[1].to("keV").value),
        "unit": "keV",
    }
    return results


def compute_contours(fit, numpoints=10):
    contours = {}

    log.info("Computing contour: (amplitude, alpha)")
    result = fit.minos_contour("amplitude", "alpha", numpoints=numpoints)
    contours["contour_amplitude_alpha"] = {
        "amplitude": result["x"].tolist(),
        "alpha": result["y"].tolist(),
    }

    log.info("Computing contour: (amplitude, beta)")
    result = fit.minos_contour("amplitude", "beta", numpoints=numpoints)
    contours["contour_amplitude_beta"] = {
        "amplitude": result["x"].tolist(),
        "beta": result["y"].tolist(),
    }

    log.info("Computing contour: (alpha, beta)")
    result = fit.minos_contour("alpha", "beta", numpoints=numpoints)
    contours["contour_alpha_beta"] = {
        "alpha": result["x"].tolist(),
        "beta": result["y"].tolist(),
    }

    return contours
