"""Run spectral fits with Gammapy and sherpa."""
import logging
from pathlib import Path
import numpy as np
from time import time
import copy
import astropy.units as u
from .models import Log10Parabola
from gammapy.spectrum import SpectrumFit
from . import utils
import sherpa.astro.datastack as sh
from .conf import config

log = logging.getLogger(__name__)


def main(dataset, tool):
    if tool in {"sherpa", "all"}:
        if dataset in {"fermi", "all"}:
            run_analysis_sherpa("fermi")
        if dataset in {"magic", "all"}:
            run_analysis_sherpa("magic")
        if dataset in {"hess", "all"}:
            run_analysis_sherpa("hess")
        if dataset in {"fact", "all"}:
            run_analysis_sherpa("fact")
        if dataset in {"veritas", "all"}:
            run_analysis_sherpa("veritas")
        if dataset in {"joint", "all"}:
            run_analysis_sherpa("joint")
    if tool in {"gammapy", "all"}:
        if dataset in {"fermi", "all"}:
            run_analysis_gammapy("fermi")
        if dataset in {"magic", "all"}:
            run_analysis_gammapy("magic")
        if dataset in {"hess", "all"}:
            run_analysis_gammapy("hess")
        if dataset in {"fact", "all"}:
            run_analysis_gammapy("fact")
        if dataset in {"veritas", "all"}:
            run_analysis_gammapy("veritas")
        if dataset in {"joint", "all"}:
            run_analysis_gammapy("joint")


def run_analysis_gammapy(which):
    log.info(f"Fitting with Gammapy the {which} dataset")
    # write the output of the fit in a yaml file
    path = Path(f"{config.repo_path}/results/fit/gammapy/{which}")
    path.parent.mkdir(parents=True, exist_ok=True)
    log.info(f"results will be stored in: {path}")

    dataset = config.get_dataset(which)
    obs_list = dataset.get_SpectrumObservationList()
    fit_range = dataset.energy_range
    log.info(f"obs_list: {obs_list}")

    results = run_gammapy_fit(obs_list, fit_range, eval_contours=path)

    # extract the yaml with the results and write it
    utils.write_yaml(results, f"{path}/fit_results_logparabola.yaml")


def extract_spectrum_results_gammapy(fit):
    results = {}
    results["parameters"] = []
    # fit.result[0].model.parameters.covariance.to_list_of_dict() removed in one of the last gammapy PRs
    # reimplemented manually here
    params = fit.result[0].model.parameters
    for param_index, parameter in enumerate(params):
        vals = parameter.to_dict()
        error = np.sqrt(
            fit.result[0].model.parameters.covariance[param_index][param_index]
        )
        vals["error"] = float(error)
        results["parameters"].append(vals)
    results["covariance"] = fit.result[0].model.parameters.covariance.tolist()
    results["statname"] = fit.result[0].statname
    results["statval"] = float(fit.result[0].statval)
    results["fit_range"] = {
        "min": float(fit.result[0].fit_range[0].to("keV").value),
        "max": float(fit.result[0].fit_range[1].to("keV").value),
        "unit": "keV",
    }
    return results


def run_gammapy_fit(obs_list, fit_range, eval_contours=None):
    """Run fit with Gammapy, using iminiuit"""
    model_lp = Log10Parabola(
        amplitude=3.80 * 1e-11 * u.Unit("cm-2 s-1 TeV-1"),
        reference=1 * u.Unit("TeV"),
        alpha=2.47 * u.Unit(""),
        beta=0.24 * u.Unit(""),
    )

    # note this step is very important to iminuit, to initialize the parameter error!
    fit = SpectrumFit(obs_list=obs_list, model=model_lp, fit_range=fit_range)

    log.info("Starting fit ...")
    t = time()
    fit.optimize(backend="minuit", pedantic=True)
    fit.run()
    t = time() - t
    log.info(f"Finished fit in {t} seconds.")

    results = extract_spectrum_results_gammapy(fit)
    print(fit.result[0])

    if eval_contours is not None:
        log.info("storing likelihood contours ...")
        # points along the contour
        like_points = 80
        sigma = 1.

        log.info(f"computing amplitude vs alpha {sigma} sigma contour")
        cont = fit.minuit.mncontour(
            "par_000_amplitude", "par_002_alpha", numpoints=like_points, sigma=sigma
        )
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
        cont = fit.minuit.mncontour(
            "par_000_amplitude", "par_003_beta", numpoints=like_points, sigma=sigma
        )
        cont = np.asarray(cont[2])
        amplitude = cont.T[0]  # transpose and take the first row
        beta = cont.T[1]  # transpose and take the
        # trick to make a close circle when plotting: just repeat the first coordinate
        amplitude = np.append(amplitude, amplitude[0])
        beta = np.append(beta, beta[0])
        contour_amplitude_beta = {"amplitude": amplitude, "beta": beta}

        log.info(f"computing alpha vs beta {sigma} sigma contour")
        cont = fit.minuit.mncontour(
            "par_002_alpha", "par_003_beta", numpoints=like_points, sigma=sigma
        )
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
        logging.info(f"storing .yaml with contours in {eval_contours}")
        path = Path(eval_contours)
        path.mkdir(parents=True, exist_ok=True)
        np.save(f"{eval_contours}/fit_{sigma}_sigma_contours_logparabola.npy", contours)

        return results


def run_analysis_sherpa(which):
    log.info(f"Fitting with sherpa the {which} dataset")
    path = f"{config.repo_path}/results/fit/sherpa/{which}"
    Path(path).mkdir(parents=True, exist_ok=True)

    dataset = config.get_dataset(which)
    datastack = dataset.get_sherpa_datastack()

    fit = run_sherpa_fit(datastack, plot_fit=path, eval_contours=path)

    # fit is the dict with the results of the fit
    results = extract_spectrum_results_sherpa(fit)

    # write the output of the fit in a yaml file
    path = Path(
        f"{config.repo_path}/results/fit/sherpa/{which}/fit_results_logparabola.yaml"
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    utils.write_yaml(results, path)

    return fit


def extract_spectrum_results_sherpa(fit):
    """function to produce a dictionary/yaml output comparable with the output
    of SpectrumFit"""
    # names compatible with gammapy
    names = ["alpha", "beta", "amplitude"]
    # units
    units = ["", "", "keV-1 cm-2 s-1"]
    # gammapy results are stored as a list of dictionaries
    parameters = []
    for idx, (name, value, unit, _min, _max) in enumerate(
        zip(names, fit["parvals"], units, fit["parmins"], fit["parmaxes"])
    ):
        _dict = {
            "name": name,
            "value": float(value),
            "unit": unit,
            "frozen": False,
            "min": float(_min),
            "max": float(_max),
            "error": float(np.sqrt(fit["extra_output"][idx][idx])),
        }
        parameters.append(_dict)

    # dictionary with info to be stored
    # we'll create a dictionary similar to the one produced by SpectrumFit
    results = dict()
    # and introduce 'reference' that is missing in the fit
    reference = {
        "name": "reference",
        "value": 1,
        "unit": "TeV",
        "frozen": True,
        "min": float("nan"),
        "max": float("nan"),
        "error": 0,
    }
    parameters.append(reference)
    # we store them in the same way they are stored in gammapy
    # convert the normalization
    parameters[2]["value"] *= 1e9  # (keV-1 m-2 s-1) -> (TeV-1 m-2 s-1)
    parameters[2]["min"] *= 1e9  # (keV-1 m-2 s-1) -> (TeV-1 m-2 s-1)
    parameters[2]["max"] *= 1e9  # (keV-1 m-2 s-1) -> (TeV-1 m-2 s-1)
    parameters[2]["unit"] = "1 / (cm2 s TeV)"
    parameters[2]["error"] *= 1e9
    results["parameters"] = [parameters[2], parameters[3], parameters[0], parameters[1]]
    # covariance matrix
    covariance = fit["extra_output"]
    # from sherpa: order of the parameters is alpha, beta, amplitude
    # all terms involving amplitude in the covariance have to be converted
    # to TeV-1 cm-2 s-1, note covariance has same units of quantities used for
    # calculation
    sigma_alpha_alpha = float(covariance[0][0])
    sigma_beta_beta = float(covariance[1][1])
    sigma_ampl_ampl = float(covariance[2][2] * (1e9 * 1e9))
    sigma_alpha_beta = float(covariance[0][1])
    sigma_alpha_ampl = float(covariance[0][2] * 1e9)
    sigma_beta_ampl = float(covariance[1][2] * 1e9)

    # write a covariance like the one in gammapy
    # oreder of the parameters is amplitude, reference, alpha, beta
    # reference is fixed so terms associated with its covariance will be null
    results["covariance"] = [
        [sigma_ampl_ampl, 0, sigma_alpha_ampl, sigma_beta_ampl],
        [0, 0, 0, 0],
        [sigma_alpha_ampl, 0, sigma_alpha_alpha, sigma_alpha_beta],
        [sigma_beta_ampl, 0, sigma_alpha_beta, sigma_beta_beta],
    ]
    results["statname"] = fit["statname"]
    results["statval"] = fit["statval"]
    # last keyword is the energy range
    results["fit_range"] = {
        "min": float(fit["fit_range"][0]),
        "max": float(fit["fit_range"][1]),
        "unit": "keV",
    }

    return results


def run_sherpa_fit(data, plot_fit=None, eval_contours=None):
    """Perform a spectrum fit using sherpa
    http://cxc.harvard.edu/sherpa/ahelp/fit.html
    """
    data.show_stack()
    # define the source model
    data.set_source("logparabola.p1")

    # Change reference energy of the model
    p1.ref = 1e9  # 1 TeV = 1e9 keV
    p1.c1 = 2.0
    p1.c2 = 0.5
    p1.ampl = 1e-20  # in cm**-2 s**-1 keV**-1
    # view parameters
    print(p1)

    # define the statistic
    sh.set_stat("wstat")

    # we retrieve the ids of the observations in the datastack
    # this is useful both to plot single run fit and also for the covariance
    # estimation
    data_dict = data.dataset_ids
    obs_ids = []
    for key in data_dict:
        obs_ids.append(data_dict[key]["id"])
    print("found observation ids", obs_ids)

    # ignore the bins above the energy threshold
    for _id in obs_ids:
        sh.ignore_bad(_id)

    # run the fit
    data.fit()

    # produce diagnostic plots, data and model per each run
    if plot_fit is not None:
        import matplotlib
        import matplotlib.pyplot as plt

        matplotlib.use("agg")
        for idx, _id in enumerate(obs_ids):
            sh.set_analysis(_id, "energy", "counts", factor=0)
            sh.get_data_plot_prefs()["xlog"] = True
            sh.get_data_plot_prefs()["ylog"] = True
            sh.plot_fit(id=_id)
            sh.plot_fit_resid(id=_id)
            # TODO: simplify obs_id handling!
            dl3_obs_id = data.datasets[idx]["data"].header["OBS_ID"]
            path_plot_fit = plot_fit + "/plot_fit_sherpa_obs_id_{}.png".format(
                dl3_obs_id
            )
            plt.savefig(path_plot_fit)

    # evaluate covariance and errors
    sh.covariance(*obs_ids)
    covar = sh.get_covar_results()

    # store the output in a dictionary
    results = dict()
    # we keep sherpa nomenclature
    results["parnames"] = covar.parnames  # parameter names, tuple
    results["parvals"] = covar.parvals  # parameter values, tuple
    results["parmins"] = covar.parmins  # parameter min error, tuple
    results["parmaxes"] = covar.parmaxes  # parameter max error, tuple
    results["extra_output"] = covar.extra_output  # covariance matrix, numpy array

    status = sh.get_fit_results()
    results["statname"] = status.statname
    results["statval"] = status.statval
    # for the energy range, take the x axis of the plot of the first id
    fit_range = sh.get_fit_plot(obs_ids[0]).dataplot.x
    results["fit_range"] = (fit_range[0], fit_range[-1])

    # dictionary for the contours, it will be empty
    # if eval_contours = False
    contours = dict()

    if eval_contours:
        # evaluate the confidence contours of two thawed parameters
        contour_ampl_c1 = dict()
        contour_ampl_c2 = dict()
        contour_c1_c2 = dict()

        # dimension of the grid to scan the confidence contours
        nloop = 50
        # ranges of the parameters
        c1_range = [1.5, 3.5]
        c2_range = [-0.2, 0.8]
        ampl_range = [2. * 1e-20, 6. * 1e-20]

        # amplitude vs c1
        sh.reg_proj(
            p1.ampl,
            p1.c1,
            id=obs_ids[0],
            otherids=obs_ids[1:],
            sigma=[1, 2, 3],
            min=[ampl_range[0], c1_range[0]],
            max=[ampl_range[1], c1_range[1]],
            nloop=(nloop, nloop),
        )

        # tmp_proj is an object we will use to store the get_reg_proj() method
        # info on the output of this method
        # http://cxc.harvard.edu/sherpa/ahelp/get_reg_proj.html
        tmp_proj = sh.get_reg_proj()
        # reshape in matrix form
        tmp_proj.y.shape = (nloop, nloop)
        # from now on we use the copy method because if we declare a variable
        # to point at the output of `sh.get_reg_proj()` this varaible
        # change when we rerun the `sh.get_reg_proj()` method on another
        # couple of parameters
        contour_ampl_c1["like_values"] = tmp_proj.y.copy()
        # the x0 array repeats its values every 20
        contour_ampl_c1["x0"] = tmp_proj.x0[:nloop].copy()
        # the x1 array repeats its values
        contour_ampl_c1["x1"] = tmp_proj.x1[:nloop].copy()
        contour_ampl_c1["levels"] = tmp_proj.levels.copy()
        # store also the parameter range we have investigated
        contour_ampl_c1["x0_range"] = (ampl_range[0], ampl_range[1])
        contour_ampl_c1["x1_range"] = (c1_range[0], c1_range[1])

        # amplitude vs c2
        sh.reg_proj(
            p1.ampl,
            p1.c2,
            id=obs_ids[0],
            otherids=obs_ids[1:],
            sigma=[1, 2, 3],
            min=[ampl_range[0], c2_range[0]],
            max=[ampl_range[1], c2_range[1]],
            nloop=(nloop, nloop),
        )

        tmp_proj = sh.get_reg_proj()
        # reshape in matrix form
        tmp_proj.y.shape = (nloop, nloop)
        contour_ampl_c2["like_values"] = tmp_proj.y.copy()
        contour_ampl_c2["x0"] = tmp_proj.x0[:nloop].copy()
        contour_ampl_c2["x1"] = tmp_proj.x1[:nloop].copy()
        contour_ampl_c2["levels"] = tmp_proj.levels.copy()
        contour_ampl_c2["x0_range"] = (ampl_range[0], ampl_range[1])
        contour_ampl_c2["x1_range"] = (c2_range[0], c2_range[1])

        # c1 vs c2
        sh.reg_proj(
            p1.c1,
            p1.c2,
            id=obs_ids[0],
            otherids=obs_ids[1:],
            sigma=[1, 2, 3],
            min=[c1_range[0], c2_range[0]],
            max=[c1_range[1], c2_range[1]],
            nloop=(nloop, nloop),
        )

        tmp_proj = sh.get_reg_proj()
        # reshape in matrix form
        tmp_proj.y.shape = (nloop, nloop)
        contour_c1_c2["like_values"] = tmp_proj.y.copy()
        contour_c1_c2["x0"] = tmp_proj.x0[:nloop].copy()
        contour_c1_c2["x1"] = tmp_proj.x1[:nloop].copy()
        contour_c1_c2["levels"] = tmp_proj.levels.copy()
        contour_c1_c2["x0_range"] = (c1_range[0], c1_range[1])
        contour_c1_c2["x1_range"] = (c2_range[0], c2_range[1])

        # add the dictionaries with the confidence contours to the final
        # output dictionary
        contours["contour_ampl_c1"] = contour_ampl_c1
        contours["contour_ampl_c2"] = contour_ampl_c2
        contours["contour_c1_c2"] = contour_c1_c2

    if eval_contours is not None:
        # write the contours in a .npy file
        path_contours = eval_contours + "/fit_contours_logparabola.npy"
        np.save(path_contours, contours)

    # return the dictionary with the results of the fit
    return results
