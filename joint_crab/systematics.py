"""statistics for the joint-crab fit"""
import logging
from gammapy.spectrum import SpectrumFit
from gammapy.spectrum.models import SpectralModel
from gammapy.utils.fitting import Parameter, Parameters
from .models import Log10ParabolaEnergyScale
from .fit_spec import compute_contours
from .conf import config
from . import utils

log = logging.getLogger(__name__)


def main():
    """Fit data with likelihood that includes energy scale uncertainty."""
    log.info("Joint fit including systematics")

    log.info("Running fit")
    fit = SystematicsSpectrumFit()
    result = fit.run()

    fit_results = make_results_dict(result)
    utils.write_yaml(fit_results, f"results/fit/fit_joint_systematics.yaml")

    contour_results = compute_contours(fit)
    utils.write_yaml(contour_results, f"results/fit/contours_joint_systematics.yaml")


class SystematicsSpetrumModel(SpectralModel):
    def __init__(self):
        self.parameters = Parameters(
            [
                Parameter("amplitude", 3e-12, unit="cm-2 s-1 TeV-1"),
                Parameter("reference", 1, unit="TeV", frozen=True),
                Parameter("alpha", 2.4, min=1, max=5),
                Parameter("beta", 0.2, min=0.001, max=1),
                Parameter("z_fermi", 0),
                Parameter("z_magic", 0),
                Parameter("z_veritas", 0),
                Parameter("z_fact", 0),
                Parameter("z_hess", 0),
            ]
        )

    def make_model(self):
        model = Log10ParabolaEnergyScale()
        z = Parameter("z", 0, min=-0.3, max=0.3)
        parameters = self.parameters.parameters[:4] + [z]
        model.parameters = Parameters(parameters)
        return model


class SystematicsSpectrumFit(SpectrumFit):
    delta = {"fermi": 0.05, "magic": 0.15, "hess": 0.15, "veritas": 0.15, "fact": 0.15}

    def __init__(self):
        self._model = SystematicsSpetrumModel()
        self.fits = {name: self.make_fit(name) for name in config.all_datasets}

    def make_fit(self, name):
        dataset = config.datasets[name]
        observations = dataset.get_SpectrumObservationList()
        model = self._model.make_model()
        fit_range = dataset.energy_range
        return SpectrumFit(obs_list=observations, model=model, fit_range=fit_range)

    def total_stat(self, parameters):
        stat = 0

        for name, fit in self.fits.items():
            # Parameter update
            pars = fit._model.parameters
            pars["amplitude"].value = parameters["amplitude"].value
            pars["alpha"].value = parameters["alpha"].value
            pars["beta"].value = parameters["beta"].value
            pars["z"].value = parameters["z_" + name].value

            # Poisson likelihood
            stat += fit.total_stat(pars)

            # Gaussian likelihood (energy scale systematics)
            stat += (pars["z"].value / self.delta[name]) ** 2

        return stat


def make_results_dict(fit_result):
    pars = fit_result.model.parameters

    results = {"parameters": []}
    for par in pars:
        par_info = par.to_dict()
        par_info["error"] = float(pars.error(par))
        results["parameters"].append(par_info)

    results["covariance"] = pars.covariance.tolist()

    return results
