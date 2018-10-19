"""Calculate theerrors on the fitted quantities."""
import numpy as np
import astropy.units as u
from pathlib import Path
from scipy.stats import norm
from astropy.table import Table
from .utils import load_yaml, write_yaml
from .models import Log10Parabola
from .conf import config
import logging
from .models import Log10ParabolaEnergyScale
from gammapy.spectrum import CountsPredictor, SpectrumObservationList
from gammapy.stats.fit_statistics import wstat
from iminuit import Minuit

log = logging.getLogger(__name__)


def stat_errorband(which, tool, dim_sample, energy_points, sigma):
    """compute the statistical error on the flux for a dataset

    Parameters
    ----------
    which : string
        dataset whose error on the fit has to be calculated
    tool : {'gammapy', 'sherpa'}
        which results have to be used, gammapy (iMinuit) or sherpa
    sigma : int
        number of sigma to consider for the confidence contour
    dim_sample : int
        dimension of the sample, i.e. how many times to sample from
        the multivariate distribution
    energy_points : int
        number of energy points we want to generate

    Returns
    -------
    `~astropy.units.Quantity` array of energies and
    corresponding lower and upper error on the flux
    """
    # read the output of the fitter
    result_file = f"results/fit/{tool}/{which}/fit_results_logparabola.yaml"
    log.info(
        f"estimating statistical error for {which} dataset, using {tool} fit results"
    )

    results = load_yaml(result_file)

    parameters = results["parameters"]
    amplitude = parameters[0]["value"] * u.Unit(parameters[0]["unit"])
    reference = parameters[1]["value"] * u.Unit(parameters[1]["unit"])
    alpha = parameters[2]["value"] * u.Unit(parameters[2]["unit"])
    beta = parameters[3]["value"] * u.Unit(parameters[3]["unit"])

    covariance = np.asarray(results["covariance"])

    # sample from a multivariate having as mean values the fitted results and covariance terms defined
    # by the covariance matrix of the fit
    pars = np.asarray([amplitude.value, reference.value, alpha.value, beta.value])
    sampled_amplitude, sampled_reference, sampled_alpha, sampled_beta = np.random.multivariate_normal(
        pars, covariance, dim_sample
    ).T
    # now evaluate the model on a series of energy points
    # first fetch the energy range defined for this dataset
    dataset = config.get_dataset(which)
    energy_range = dataset.energy_range
    energy_unit = energy_range[0].unit
    energies = (
        np.logspace(
            np.log10(energy_range[0].value),
            np.log10(energy_range[1].value),
            energy_points,
        )
        * energy_unit
    )

    # empty lists to determine the quantiles
    flux_unit = parameters[0]["unit"]
    flux_min = []
    flux_max = []

    # this dictionary will be used for debugging purpose
    # we will store the fluxes at the extreme and at a medium energy
    sampled_fluxes = {
        "emin": {"value": energies[0].value, "fluxes": [], "flux_quantiles": []},
        "emid": {
            "value": energies[int(energy_points / 2)].value,
            "fluxes": [],
            "flux_quantiles": [],
        },
        "emax": {"value": energies[-1].value, "fluxes": [], "flux_quantiles": []},
        "energy_unit": str(energy_unit),
        "flux_unit": flux_unit,
    }
    for i, ene in enumerate(energies):
        # loop through the sampled values and estimate flux at this energies
        _flux = []
        for amp, ref, a, b in zip(
            sampled_amplitude, sampled_reference, sampled_alpha, sampled_beta
        ):
            # give them the units of the sampled quantities
            amp *= u.Unit(parameters[0]["unit"])
            ref *= u.Unit(parameters[1]["unit"])
            a *= u.Unit(parameters[2]["unit"])
            b *= u.Unit(parameters[3]["unit"])

            _flux.append(
                Log10Parabola.evaluate(ene, amp, ref, a, b).to(flux_unit).value
            )
        # now to define the minimum and maximum flux let's take the n-sigma containment
        low_quantile = norm.cdf(-sigma)
        high_quantile = norm.cdf(sigma)
        _flux_min = np.percentile(np.asarray(_flux), 100 * low_quantile)
        _flux_max = np.percentile(np.asarray(_flux), 100 * high_quantile)

        # store the sampled fluxes at the extremes and at a medium energy
        if i == 0:
            sampled_fluxes["emin"]["fluxes"] = [float(_) for _ in _flux]
            sampled_fluxes["emin"]["flux_quantiles"] = [
                float(_flux_min),
                float(_flux_max),
            ]

        if i == int(energy_points / 2):
            sampled_fluxes["emid"]["fluxes"] = [float(_) for _ in _flux]
            sampled_fluxes["emid"]["flux_quantiles"] = [
                float(_flux_min),
                float(_flux_max),
            ]

        if i == energy_points - 1:
            sampled_fluxes["emax"]["fluxes"] = [float(_) for _ in _flux]
            sampled_fluxes["emax"]["flux_quantiles"] = [
                float(_flux_min),
                float(_flux_max),
            ]

        # list with the same length of the energies that will define the butterfly
        flux_min.append(_flux_min)
        flux_max.append(_flux_max)

    # to obtain the proper units
    flux_min = np.asarray(flux_min) * u.Unit(flux_unit)
    flux_max = np.asarray(flux_max) * u.Unit(flux_unit)

    sampled_dict = {
        "sampled_amplitude": [float(_) for _ in sampled_amplitude.tolist()],
        "sampled_alpha": [float(_) for _ in sampled_alpha.tolist()],
        "sampled_reference": [float(_) for _ in sampled_reference.tolist()],
        "sampled_beta": [float(_) for _ in sampled_beta.tolist()],
        "sampled_fluxes": sampled_fluxes,
    }

    path = Path(
        f"{config.repo_path}/results/debug/stat-err/{tool}/{which}/multivariate_sampling_fluxes.yaml"
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    write_yaml(sampled_dict, path)

    # we also save the upper and lower limit in an astropy table in order to just read them when we do the
    # sed plot, instead of recomputing them each time (very time consuming due to multivariate sampling)
    t = Table(
        [energies.value, flux_min.value, flux_max.value],
        names=("energies", "flux_lo", "flux_hi"),
        meta={"name": "flux error band"},
    )
    t["energies"].unit = energy_unit
    t["flux_lo"].unit = u.Unit(flux_unit)
    t["flux_hi"].unit = u.Unit(flux_unit)
    table_path = Path(
        f"{config.repo_path}/results/figures/stat_err/{which}_flux_errorband.dat"
    )
    table_path.parent.mkdir(parents=True, exist_ok=True)
    t.write(table_path, format="ascii.ecsv")

    # return them
    return energies, flux_min, flux_max


def systematic():
    """compute the systematic error only on the joint fit
    """
    log.info("running systematic error estimation")
    # let us read the SpectrumObservationList from the ogips, we will use them for the ON and OFF counts
    # and to provide the IRFs to compute the predicted counts

    # fitter dict will contain the SpectrumObservationLists and the bin in fit range per each dataset
    fitter = {"spec_obs_lists": {}, "bins_in_fit_range": {}}
    fitter["spec_obs_lists"]["fermi"] = SpectrumObservationList.read(
        f"{config.repo_path}/results/spectra/fermi"
    )
    fitter["spec_obs_lists"]["magic"] = SpectrumObservationList.read(
        f"{config.repo_path}/results/spectra/magic"
    )
    fitter["spec_obs_lists"]["hess"] = SpectrumObservationList.read(
        f"{config.repo_path}/results/spectra/hess"
    )
    fitter["spec_obs_lists"]["veritas"] = SpectrumObservationList.read(
        f"{config.repo_path}/results/spectra/veritas"
    )
    fitter["spec_obs_lists"]["fact"] = SpectrumObservationList.read(
        f"{config.repo_path}/results/spectra/fact"
    )

    for which in ["fermi", "magic", "hess", "veritas", "fact"]:
        log.info(f"checking for energy range of {which} dataset")

        obs_list = fitter["spec_obs_lists"][which]

        dataset = config.get_dataset(which)
        fit_range = dataset.energy_range

        _bins_in_fit_range = []

        for obs in obs_list:
            # Take into account fit range, copied from SpectrumFit class
            energy = obs.e_reco
            valid_range = np.zeros(energy.nbins)

            if fit_range is not None:

                precision = 1e-3  # to avoid floating round precision
                idx_lo = np.where(energy * (1 + precision) < fit_range[0])[0]
                valid_range[idx_lo] = 1

                idx_hi = np.where(energy[:-1] * (1 - precision) > fit_range[1])[0]
                if len(idx_hi) != 0:
                    idx_hi = np.insert(idx_hi, 0, idx_hi[0] - 1)
                valid_range[idx_hi] = 1

            # Take into account thresholds
            try:
                quality = obs.on_vector.quality
            except AttributeError:
                quality = np.zeros(obs.e_reco.nbins)

            intersection = np.logical_and(1 - quality, 1 - valid_range)

            _bins_in_fit_range.append(intersection)
        # add it to the dictionary
        fitter["bins_in_fit_range"][which] = _bins_in_fit_range

    def wstat_energy_scale(amplitude, reference, alpha, beta, z, which):
        """single instrument log likelihood"""
        # note there is a bins_in_fit range element per each obs in teh fitter dict
        spec_obs_list = fitter["spec_obs_lists"][which]
        bins_in_fit_range = fitter["bins_in_fit_range"][which]

        amplitude *= u.Unit("TeV-1 cm-2 s-1")
        reference *= u.TeV

        model = Log10ParabolaEnergyScale(
            amplitude=amplitude, reference=reference, alpha=alpha, beta=beta, z=z
        )

        total_wstat = 0

        for spec_obs, _bins_in_fit_range in zip(spec_obs_list, bins_in_fit_range):
            # compute predicted counts
            predictor = CountsPredictor(
                model=model,
                aeff=spec_obs.aeff,
                edisp=spec_obs.edisp,
                livetime=spec_obs.livetime,
            )
            predictor.run()

            # calculate wstat statistics per each osbervation
            _wstat = wstat(
                n_on=spec_obs.on_vector.data.data,
                n_off=spec_obs.off_vector.data.data,
                alpha=spec_obs.alpha,
                mu_sig=predictor.npred.data.data,
            )
            # restrict the wstat to the bins in fit range
            restricted_wstat = _wstat[_bins_in_fit_range]
            total_wstat += np.sum(restricted_wstat)

        return total_wstat

    def all_instrument_wstat(
        amplitude, reference, alpha, beta, z_fermi, z_magic, z_hess, z_veritas, z_fact
    ):
        """sum of all the instrument likelihoods with log-parabola w/ energy scale correction"""
        all_instrument_wstat = (
            wstat_energy_scale(amplitude, reference, alpha, beta, z_fermi, "fermi")
            + wstat_energy_scale(amplitude, reference, alpha, beta, z_magic, "magic")
            + wstat_energy_scale(amplitude, reference, alpha, beta, z_hess, "hess")
            + wstat_energy_scale(
                amplitude, reference, alpha, beta, z_veritas, "veritas"
            )
            + wstat_energy_scale(amplitude, reference, alpha, beta, z_fact, "fact")
            + (z_fermi / 0.10) ** 2
            + (z_magic / 0.15) ** 2
            + (z_hess / 0.15) ** 2
            + (z_veritas / 0.15) ** 2
            + (z_fact / 0.15) ** 2
        )
        return all_instrument_wstat

    # set the limits of the parameters with a dictionary, so we can reuse it when we write the yaml
    # accessing the limits by name
    limits_dict = {
        "amplitude": (1e-15, 1e-9),
        "alpha": (1, 5),
        "beta": (0.001, 1),
        "reference": (0.3, 30),
        "z_fermi": (-0.3, 0.3),
        "z_magic": (-0.3, 0.3),
        "z_hess": (-0.3, 0.3),
        "z_veritas": (-0.3, 0.3),
        "z_fact": (-0.3, 0.3),
    }

    # minimize the likelihood
    m = Minuit(
        all_instrument_wstat,
        amplitude=3e-12,
        reference=1,
        alpha=2.4,
        beta=0.2,
        z_fermi=0.0,
        z_magic=0.0,
        z_hess=0.0,
        z_veritas=0.0,
        z_fact=0.0,
        fix_reference=True,
        # fix_z_fermi = True,
        # fix_z_magic=True,
        # fix_z_hess=True,
        # fix_z_veritas=True,
        # fix_z_fact=True,
        error_amplitude=1e-13,
        error_alpha=0.01,
        error_beta=0.001,
        error_reference=0.,
        error_z_fermi=0.001,
        error_z_magic=0.001,
        error_z_hess=0.001,
        error_z_veritas=0.001,
        error_z_fact=0.001,
        limit_amplitude=limits_dict["amplitude"],
        limit_alpha=limits_dict["alpha"],
        limit_beta=limits_dict["beta"],
        limit_reference=limits_dict["reference"],
        limit_z_fermi=limits_dict["z_fermi"],
        limit_z_magic=limits_dict["z_magic"],
        limit_z_hess=limits_dict["z_hess"],
        limit_z_veritas=limits_dict["z_veritas"],
        limit_z_fact=limits_dict["z_fact"],
    )

    m.migrad()
    m.hesse()
    m.minos()

    # write the output in a yaml format similar to the one used for the gammapy - sherpa fit
    results = {}  # dict to be written as a yaml
    parameters = []  # list of parameters
    for name, value, error, fixed in zip(
        m.values.keys(), m.values.values(), m.errors.values(), m.fixed.values()
    ):
        # the output with the modified likelihood does not have astropy units
        unit = ""
        if name == "amplitude":
            unit = "1 / (cm2 s TeV)"
        if name == "reference":
            unit = "TeV"
        # dictionary to be appended to parameters
        param_dict = {
            "name": name,
            "value": value,
            "unit": unit,
            "frozen": fixed,
            "min": limits_dict[name][0],
            "max": limits_dict[name][1],
            "error": error,
        }
        parameters.append(param_dict)

    results["parameters"] = parameters

    # define the covariance matrix, the output of m.covariance is a dictionary, convert it to a numpy matrix
    n_params = len(m.values)
    covariance = np.zeros((n_params, n_params), float)

    for i, i_name in enumerate(m.values.keys()):
        for j, j_name in enumerate(m.values.keys()):
            # reference is not in iminuit covariance matrix being fixed
            # we add a 0 in the matrix manually
            if (i_name == "reference") or (j_name == "reference"):
                covariance[i, j] = 0.
            else:
                covariance[i, j] = m.covariance[(i_name, j_name)]

    results["covariance"] = covariance.tolist()
    # stat info
    results["statname"] = "wstat"
    results["staval"] = m.fval

    path = Path(
        f"{config.repo_path}/results/fit/gammapy/joint/fit_results_logparabola_energy_scale.yaml"
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    write_yaml(results, path)

    # store the contours
    like_points = 80
    sigma = 1.

    log.info(f"computing amplitude vs alpha {sigma} sigma contour")
    cont = m.mncontour("amplitude", "alpha", numpoints=like_points, sigma=sigma)
    # the third element of mncontour's returned object is a list of tuples with the contour coordinates
    # (x_1,y_1), ..., (x_n, y_n)]
    cont = np.asarray(cont[2])
    amplitude = cont.T[0]  # transpose and take the first row
    alpha = cont.T[1]  # transpose and take the
    # trick to make a close circle when plotting: just repeat the first coordinate
    amplitude = np.append(amplitude, amplitude[0])
    alpha = np.append(alpha, alpha[0])
    contour_amplitude_alpha = {"amplitude": amplitude, "alpha": alpha}

    log.info(f"computing amplitude vs beta {sigma} sigma contour")
    cont = m.mncontour("amplitude", "beta", numpoints=like_points, sigma=sigma)
    cont = np.asarray(cont[2])
    amplitude = cont.T[0]  # transpose and take the first row
    beta = cont.T[1]  # transpose and take the
    # trick to make a close circle when plotting: just repeat the first coordinate
    amplitude = np.append(amplitude, amplitude[0])
    beta = np.append(beta, beta[0])
    contour_amplitude_beta = {"amplitude": amplitude, "beta": beta}

    log.info(f"computing alpha vs beta {sigma} sigma contour")
    cont = m.mncontour("alpha", "beta", numpoints=like_points, sigma=sigma)
    cont = np.asarray(cont[2])
    alpha = cont.T[0]  # transpose and take the first row
    beta = cont.T[1]  # transpose and take the
    # trick to make a close circle when plotting: just repeat the first coordinate
    alpha = np.append(alpha, alpha[0])
    beta = np.append(beta, beta[0])
    contour_alpha_beta = {"alpha": alpha, "beta": beta}

    # define the general dictionary and dump it in a .npy object
    contours = {
        "contour_amplitude_alpha": contour_amplitude_alpha,
        "contour_amplitude_beta": contour_amplitude_beta,
        "contour_alpha_beta": contour_alpha_beta,
    }
    outpath = f"{config.repo_path}/results/fit/gammapy/joint/fit_{sigma}_sigma_contours_logparabola_energy_scale.npy"
    logging.info(f"storing .yaml with contours in {outpath}")
    np.save(outpath, contours)


def syst_errorband(dim_sample=500, energy_points=60, sigma=1):
    """compute the butterfly given the result of the joint fit w/ systematics"""
    results = load_yaml(
        f"{config.repo_path}/results/fit/gammapy/joint/fit_results_logparabola_energy_scale.yaml"
    )

    parameters = results["parameters"]
    amplitude = parameters[0]["value"] * u.Unit(parameters[0]["unit"])
    reference = parameters[1]["value"] * u.Unit(parameters[1]["unit"])
    alpha = parameters[2]["value"] * u.Unit(parameters[2]["unit"])
    beta = parameters[3]["value"] * u.Unit(parameters[3]["unit"])

    # we do not care about the nuisance parameters so we choose only the first 4 lines / rows of the matrix
    covariance = np.asarray(results["covariance"])[:4, :4]

    # sample from the multivariate
    pars = np.asarray([amplitude.value, reference.value, alpha.value, beta.value])
    sampled_amplitude, sampled_reference, sampled_alpha, sampled_beta = np.random.multivariate_normal(
        pars, covariance, dim_sample
    ).T
    # now evaluate the model on a series of energy points
    # fetch the energy range for this dataset
    dataset = config.get_dataset("joint")
    energy_range = dataset.energy_range
    energy_unit = energy_range[0].unit
    energies = (
        np.logspace(
            np.log10(energy_range[0].value),
            np.log10(energy_range[1].value),
            energy_points,
        )
        * energy_range[0].unit
    )

    # determine the quantiles
    flux_min = []
    flux_max = []

    for ene in energies:
        # loop through the sampled values and estimate flux at this energies
        _flux = []
        for amp, ref, a, b in zip(
            sampled_amplitude, sampled_reference, sampled_alpha, sampled_beta
        ):
            # give them the units of the sampled quantities
            amp *= u.Unit(parameters[0]["unit"])
            ref *= u.Unit(parameters[1]["unit"])
            a *= u.Unit(parameters[2]["unit"])
            b *= u.Unit(parameters[3]["unit"])

            _flux.append(Log10Parabola.evaluate(ene, amp, ref, a, b).value)

        # now to define the minimum and maximum flux let's take the n-sigma containment
        low_quantile = norm.cdf(-sigma)
        high_quantile = norm.cdf(sigma)
        _flux_min = np.percentile(np.asarray(_flux), 100 * low_quantile)
        _flux_max = np.percentile(np.asarray(_flux), 100 * high_quantile)

        flux_min.append(_flux_min)
        flux_max.append(_flux_max)

    # to obtain the proper units
    flux_unit = Log10Parabola.evaluate(
        energies[0], amplitude, reference, alpha, beta
    ).unit
    flux_min = np.asarray(flux_min) * u.Unit(flux_unit)
    flux_max = np.asarray(flux_max) * u.Unit(flux_unit)

    # store them in an astropy table to be read later when producing the figure
    t = Table(
        [energies.value, flux_min.value, flux_max.value],
        names=("energies", "flux_lo", "flux_hi"),
        meta={"name": "flux error band"},
    )
    t["energies"].unit = energy_unit
    t["flux_lo"].unit = u.Unit(flux_unit)
    t["flux_hi"].unit = u.Unit(flux_unit)
    table_path = Path(
        f"{config.repo_path}/results/figures/syst_err/joint_flux_errorband.dat"
    )
    table_path.parent.mkdir(parents=True, exist_ok=True)
    t.write(table_path, format="ascii.ecsv")

    # return them
    return energies, flux_min, flux_max
