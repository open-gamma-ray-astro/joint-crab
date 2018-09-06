"""statistics for the joint-crab fit"""
import logging
import numpy as np
import astropy.units as u
from iminuit import Minuit
from gammapy.spectrum import CountsPredictor, SpectrumObservationList
from gammapy.stats.fit_statistics import wstat
from .models import Log10ParabolaEnergyScale
from .conf import config
from .utils import write_yaml

log = logging.getLogger(__name__)

# FITTER dict will contain the SpectrumObservationList and the bin in fit range
# per each dataset, it is a dictionary we define here otherwise
# each time all_instruments_wstat_energy_scale is called from the
# fit optimizer would have to read and create the SpectrumObservationList
FITTER = {"spec_obs_lists": {}, "bins_in_fit_range": {}}
FITTER["spec_obs_lists"]["fermi"] = SpectrumObservationList.read(
    f"results/spectra/fermi"
)
FITTER["spec_obs_lists"]["magic"] = SpectrumObservationList.read(
    f"results/spectra/magic"
)
FITTER["spec_obs_lists"]["hess"] = SpectrumObservationList.read(f"results/spectra/hess")
FITTER["spec_obs_lists"]["veritas"] = SpectrumObservationList.read(
    f"results/spectra/veritas"
)
FITTER["spec_obs_lists"]["fact"] = SpectrumObservationList.read(f"results/spectra/fact")

for which in config.all_datasets:
    obs_list = FITTER["spec_obs_lists"][which]
    dataset = config.datasets[which]
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
    FITTER["bins_in_fit_range"][which] = _bins_in_fit_range


def wstat_energy_scale(
    amplitude, reference, alpha, beta, z_instr, spec_obs_list, bins_in_fit_range
):
    """single instrument -2 log likelihood for the joint-fit w/ systematic
    uncertainties.
    Parameters
    ----------
    amplitude : float
        amplitude of the log-parabola diff. flux in units of TeV-1 cm-2 s-1
    reference : float
        reference of the log-parabola diff. flux in TeV
    alpha : float
        index of the log-parabola diff. flux
    beta: float
        second index of the log-parabola diff. flux,
    z_instr : float
        energy scale nuisance parameter for a single instrument
    spec_obs_list . `~gammapy.spectrum.SpectrumObservationList`
        observation list for a single instrument, it will be used to obtain
        N_on and N_off for the likelihood and to predict the counts with the IRFs
    bins_in_fit_range : array of bools
        restrict the single instrument wstat to bins in fit range
    """
    amplitude *= u.Unit("TeV-1 cm-2 s-1")
    reference *= u.TeV

    model = Log10ParabolaEnergyScale(
        amplitude=amplitude, reference=reference, alpha=alpha, beta=beta, z=z_instr
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


def all_instruments_energy_scale_wstat(
    amplitude, reference, alpha, beta, z_fermi, z_magic, z_hess, z_veritas, z_fact
):
    """Global likelihood modified to account for the energy-scale systematic
    uncertainties on the single instruments. Same Parameters as wstat_energy_scale
    """
    # assumed values for standard deviation of the Gaussian constraining the nuisance parameters
    delta_fermi = 0.05
    delta_magic = 0.15
    delta_hess = 0.15
    delta_veritas = 0.15
    delta_fact = 0.15

    # sum of all the instrument likelihoods with log-parabola w/ energy scale correction
    all_instrument_wstat = (
        wstat_energy_scale(
            amplitude,
            reference,
            alpha,
            beta,
            z_fermi,
            FITTER["spec_obs_lists"]["fermi"],
            FITTER["bins_in_fit_range"]["fermi"],
        )
        + wstat_energy_scale(
            amplitude,
            reference,
            alpha,
            beta,
            z_magic,
            FITTER["spec_obs_lists"]["magic"],
            FITTER["bins_in_fit_range"]["magic"],
        )
        + wstat_energy_scale(
            amplitude,
            reference,
            alpha,
            beta,
            z_hess,
            FITTER["spec_obs_lists"]["hess"],
            FITTER["bins_in_fit_range"]["hess"],
        )
        + wstat_energy_scale(
            amplitude,
            reference,
            alpha,
            beta,
            z_veritas,
            FITTER["spec_obs_lists"]["veritas"],
            FITTER["bins_in_fit_range"]["veritas"],
        )
        + wstat_energy_scale(
            amplitude,
            reference,
            alpha,
            beta,
            z_fact,
            FITTER["spec_obs_lists"]["fact"],
            FITTER["bins_in_fit_range"]["fact"],
        )
        + (z_fermi / delta_fermi) ** 2
        + (z_magic / delta_magic) ** 2
        + (z_hess / delta_hess) ** 2
        + (z_veritas / delta_veritas) ** 2
        + (z_fact / delta_fact) ** 2
    )
    return all_instrument_wstat


def main():
    """Fit data with likelihood that includes energy scale uncertainty."""
    log.info("Joint fit including systematics")
    minuit = run_analysis()
    compute_fit_results(minuit)
    compute_contours(minuit)


def run_analysis():
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
    m = Minuit(
        all_instruments_energy_scale_wstat,
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
        error_amplitude=1e-13,
        error_alpha=0.01,
        error_beta=0.001,
        error_reference=0.0,
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
        print_level=0,
        errordef=1,
    )

    m.migrad()
    m.hesse()
    return m


def compute_fit_results(minuit):
    # write the output in a yaml format similar to the one used for the gammapy - sherpa fit
    results = {}  # dict to be written as a yaml
    parameters = []  # list of parameters
    for name, value, error, fixed in zip(
        minuit.values.keys(),
        minuit.values.values(),
        minuit.errors.values(),
        minuit.fixed.values(),
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
            # "min": limits_dict[name][0],
            # "max": limits_dict[name][1],
            "error": error,
        }
        parameters.append(param_dict)

    results["parameters"] = parameters

    # define the covariance matrix, the output of m.covariance is a dictionary, convert it to a numpy matrix
    n_params = len(minuit.values)
    covariance = np.zeros((n_params, n_params), float)

    for i, i_name in enumerate(minuit.values.keys()):
        for j, j_name in enumerate(minuit.values.keys()):
            # reference is not in iminuit covariance matrix being fixed
            # we add a 0 in the matrix manually
            if (i_name == "reference") or (j_name == "reference"):
                covariance[i, j] = 0.0
            else:
                covariance[i, j] = minuit.covariance[(i_name, j_name)]

    results["covariance"] = covariance.tolist()
    # stat info
    results["statname"] = "wstat"
    results["staval"] = minuit.fval
    path = "results/fit/fit_joint_systematics.yaml"
    write_yaml(results, path)


def compute_contours(minuit, numpoints=3):
    contours = {}

    log.info("Computing contour: (amplitude, alpha)")
    cont = minuit.mncontour("amplitude", "alpha", numpoints=numpoints)
    cont = np.asarray(cont[2])
    amplitude = cont.T[0]
    alpha = cont.T[1]
    contours["contour_amplitude_alpha"] = {
        "amplitude": amplitude.tolist(),
        "alpha": alpha.tolist(),
    }

    log.info("Computing contour: (amplitude, beta)")
    cont = minuit.mncontour("amplitude", "beta", numpoints=numpoints)
    cont = np.asarray(cont[2])
    amplitude = cont.T[0]
    beta = cont.T[1]
    contours["contour_amplitude_beta"] = {
        "amplitude": amplitude.tolist(),
        "beta": beta.tolist(),
    }

    log.info("Computing contour: (alpha, beta)")
    cont = minuit.mncontour("alpha", "beta", numpoints=numpoints)
    cont = np.asarray(cont[2])
    alpha = cont.T[0]
    beta = cont.T[1]
    contours["contour_alpha_beta"] = {"alpha": alpha.tolist(), "beta": beta.tolist()}

    path = "results/fit/contours_joint_systematics.yaml"
    write_yaml(contours, path)
