"""Calculate the errors on the fitted quantities."""
import logging
from pathlib import Path
import numpy as np
import scipy.stats
import astropy.units as u
from astropy.table import Table
from .utils import load_yaml
from .models import Log10Parabola
from .conf import config

log = logging.getLogger(__name__)


def main():
    Path("results/errorbands").mkdir(exist_ok=True, parents=True)

    datasets = list(config.all_datasets) + ["joint"]
    for dataset in datasets:
        errorband_for_dataset(dataset)


def errorband_for_dataset(dataset, size=500, seed=0):
    log.info(f"Computing errorband for {dataset}")

    fit = load_yaml(f"results/fit/fit_{dataset}.yaml")

    mean = np.array([_["value"] for _ in fit["parameters"]])
    names = [_["name"] for _ in fit["parameters"]]
    units = [u.Unit(_["unit"]) for _ in fit["parameters"]]
    cov = np.array(fit["covariance"])

    # Sample parameters
    rng = np.random.RandomState(seed=seed)
    parameter_samples = rng.multivariate_normal(mean, cov, size)

    table = Table(parameter_samples, names=names)

    # Compute fluxes
    energy = config.energy_bins
    dnde = []
    for parameter_sample in parameter_samples:
        model = Log10Parabola(
            parameter_sample[0] * units[0],
            parameter_sample[1] * units[1],
            parameter_sample[2] * units[2],
            parameter_sample[3] * units[3],
        )
        dnde.append(model(energy))
    dnde = np.array(dnde) * u.Unit("cm-2 s-1 TeV-1")

    table["energy"] = energy[None, :]
    table["dnde"] = dnde

    path = f"results/errorbands/samples_{dataset}.fits.gz"
    log.info(f"Writing {path}")
    table.write(path, overwrite=True)

    # Compute error band as flux quantiles
    sigma = 1
    sed = Table()
    sed["energy"] = energy
    sed["dnde_mean"] = np.mean(dnde, axis=0)
    sed["dnde_median"] = np.percentile(dnde.value, 50, axis=0) * dnde.unit
    q = 100 * scipy.stats.norm.cdf(-sigma)
    sed["dnde_lo"] = np.percentile(dnde.value, q, axis=0) * dnde.unit
    q = 100 * scipy.stats.norm.cdf(+sigma)
    sed["dnde_hi"] = np.percentile(dnde.value, q, axis=0) * dnde.unit

    # Compute and add best-fit SED curve
    model = Log10Parabola(
        mean[0] * units[0], mean[1] * units[1], mean[2] * units[2], mean[3] * units[3]
    )
    model.parameters.covariance = cov
    sed["dnde_fit"] = model(energy)

    for name in ["fit", "mean", "median", "lo", "hi"]:
        value = (energy ** 2) * sed[f"dnde_{name}"]
        sed[f"e2dnde_{name}"] = value.to("erg cm-2 s-1")

    for name in sed.columns:
        sed[name].format = "3g"

    path = f"results/errorbands/sed_{dataset}.ecsv"
    log.info(f"Writing {path}")
    sed.write(path, overwrite=True)
