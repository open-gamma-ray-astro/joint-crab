"""
Investigate Gammapy vs Sherpa spectral results differences.

See https://github.com/gammasky/joint-crab/issues/48
"""
from pprint import pprint
import astropy.units as u
from regions import CircleSkyRegion
from gammapy.data import DataStore
from gammapy.spectrum import SpectrumObservation, SpectrumFit
from gammapy.spectrum.models import PowerLaw
from gammapy.scripts import SpectrumAnalysisIACT
import sherpa.astro.ui as sh
from . import conf
from .extract_ogip_spectra import get_exclusion_mask

obs_path = "results/spectra/hess/pha_obs23523.fits"
energy_range = [1, 30] * u.TeV


def fit_gammapy():
    """
    Current results

    Parameters:

           name     value     error         unit      min max frozen
        --------- --------- --------- --------------- --- --- ------
            index 2.602e+00 1.555e-01                 nan nan  False
        amplitude 2.441e-11 3.452e-12 1 / (cm2 s TeV) nan nan  False
        reference 1.000e+00 0.000e+00             TeV nan nan   True

    Covariance:

        name/name  index   amplitude
        --------- -------- ---------
            index   0.0242  3.79e-13
        amplitude 3.79e-13  1.19e-23

    Statistic: -157.719 (cash)
    Fit Range: [1.0000000e+09 2.7825594e+10] keV


    """
    obs = SpectrumObservation.read(obs_path)
    # obs.peek()
    # plt.show()
    model = PowerLaw(
        amplitude=1.23 * 1e-11 * u.Unit("cm-2 s-1 TeV-1"),
        reference=1 * u.Unit("TeV"),
        index=2.14 * u.Unit(""),
    )

    fit = SpectrumFit(obs_list=obs, model=model, fit_range=energy_range, stat="cash")
    fit.run()
    print(fit.result[0])
    pprint(fit.__dict__)

    obs = fit.obs_list[0]
    print(obs)
    print("This is fit_gammapy")
    obs.peek()
    import matplotlib.pyplot as plt

    plt.savefig("fit_gammapy.png")


def fit_sherpa():
    sh.load_data(obs_path)

    sh.set_source("powlaw1d.model")
    model.ref = 1e9  # 1 TeV = 1e9 keV
    model.ampl = 1.23e-20  # in cm**-2 s**-1 keV**-1
    model.gamma = 2.8

    sh.set_stat("cash")
    sh.notice(energy_range[0].to("keV").value, energy_range[1].to("keV").value)
    sh.fit()
    print(sh.get_fit_results())


def check_edisp_normalisation():
    import matplotlib.pyplot as plt

    obs = SpectrumObservation.read(obs_path)
    p = obs.edisp.data.data.sum(axis=1)
    e = obs.edisp.data.axis("e_true").nodes
    plt.plot(e, p)
    plt.semilogx()
    plt.show()
    # import IPython; IPython.embed()


def check_npred_vs_nobs():
    pass


def check_energy_binning_effects():
    """Check how spectral fit results change with energy binnings.

    Actually this is still using the default:

    In [14]: print(analysis.extraction.observations[0].edisp)
    EnergyDispersion
    NDDataArray summary info
    e_true         : size =   108, min =  0.010 TeV, max = 301.416 TeV
    e_reco         : size =    72, min =  0.011 TeV, max = 93.804 TeV
    Data           : size =  7776, min =  0.000, max =  1.000

    But now, the fit results from SpectrumAnalysisIACT, which is just
    driving Fitspectrum, are different, almost the same as Sherpa.


    Current results

    Parameters:

           name     value     error         unit      min max frozen
        --------- --------- --------- --------------- --- --- ------
            index 2.620e+00 1.540e-01                 nan nan  False
        amplitude 3.849e-11 5.407e-12 1 / (cm2 s TeV) nan nan  False
        reference 1.000e+00 0.000e+00             TeV nan nan   True

    Covariance:

        name/name  index   amplitude
        --------- -------- ---------
            index   0.0237  5.85e-13
        amplitude 5.85e-13  2.92e-23

    Statistic: -157.166 (cash)
    Fit Range: [ 1.         27.82559402] TeV


    ???
    """
    data_store = DataStore.from_dir("data/hess")
    obs_list = data_store.obs_list([23523])
    on_region = CircleSkyRegion(conf.crab_position, conf.on_radius["hess"])
    fp_binning = [1, 10, 30] * u.TeV
    exclusion_mask = get_exclusion_mask(conf.crab_position)
    model = PowerLaw(
        amplitude=1.23 * 1e-11 * u.Unit("cm-2 s-1 TeV-1"),
        reference=1 * u.Unit("TeV"),
        index=2.14 * u.Unit(""),
    )
    cfg = dict(
        outdir=None,
        background=dict(
            on_region=on_region,
            exclusion_mask=exclusion_mask,
            # min_distance=0.1 * u.rad,
        ),
        extraction=dict(containment_correction=True),
        fit=dict(
            model=model,
            stat="cash",
            # forward_folded=True,
            fit_range=energy_range,
        ),
        fp_binning=fp_binning,
    )
    analysis = SpectrumAnalysisIACT(obs_list, cfg)
    analysis.run()
    analysis.fit.est_errors()
    print("BEFORE", analysis.fit.result[0])
    # print(analysis.extraction.observations[0].edisp)
    # pprint(analysis.fit.__dict__)
    # obs = analysis.fit.obs_list[0]
    # print(obs)
    # # import IPython; IPython.embed()
    # obs.peek()
    # import matplotlib.pyplot as plt
    # plt.savefig('check_energy_binning_effects.png')
    # print('This is check_energy_binning_effects')

    # Try to see if the I/O causes a change in results.
    analysis.extraction.observations.write("temp123")
    analysis.extraction.observations.write("temp123", use_sherpa=True)
    obs = SpectrumObservation.read("temp123/pha_obs23523.fits")
    model = PowerLaw(
        amplitude=1.23 * 1e-11 * u.Unit("cm-2 s-1 TeV-1"),
        reference=1 * u.Unit("TeV"),
        index=2.14 * u.Unit(""),
    )

    fit = SpectrumFit(obs_list=obs, model=model, fit_range=energy_range, stat="cash")
    fit.run()
    print("AFTER", fit.result[0])


if __name__ == "__main__":
    # fit_gammapy()
    # fit_sherpa()
    # check_edisp_normalisation()
    check_energy_binning_effects()

"""

This is what we currently get.
The results are different, for unknown reasons.
Energy dispersion matrix looks to be normalised OK, i.e. to sum to 1 along the e_reco axis.

$ python joint_crab/check_gammapy_sherpa_spectrum_differences.py
WARNING: imaging routines will not be available, 
failed to import sherpa.image.ds9_backend due to 
'RuntimeErr: DS9Win unusable: Could not find ds9 on your PATH'
WARNING: failed to import WCS module; WCS routines will not be available
WARNING: failed to import sherpa.astro.xspec; XSPEC models will not be available
/Users/deil/code/gammapy/gammapy/stats/fit_statistics.py:55: RuntimeWarning: divide by zero encountered in log
  stat = 2 * (mu_on - n_on * np.log(mu_on))

Fit result info 
--------------- 
Model: PowerLaw

Parameters: 

	   name     value     error         unit      min max frozen
	--------- --------- --------- --------------- --- --- ------
	    index 2.602e+00 1.555e-01                 nan nan  False
	amplitude 2.441e-11 3.452e-12 1 / (cm2 s TeV) nan nan  False
	reference 1.000e+00 0.000e+00             TeV nan nan   True

Covariance: 

	name/name  index   amplitude
	--------- -------- ---------
	    index   0.0242  3.79e-13
	amplitude 3.79e-13  1.19e-23 

Statistic: -157.719 (cash)
Fit Range: [1.0000000e+09 2.7825594e+10] keV

read ARF file data/hess/ogip-spectra/arf_obs23523.fits
read RMF file data/hess/ogip-spectra/rmf_obs23523.fits
read background file data/hess/ogip-spectra/bkg_obs23523.fits
WARNING: data set 1 has associated backgrounds, but they have not been subtracted, nor have background models been set
Dataset               = 1
Method                = levmar
Statistic             = cash
Initial fit statistic = -48.8737
Final fit statistic   = -156.778 at function evaluation 202
Data points           = 27
Degrees of freedom    = 25
Change in statistic   = 107.905
   model.gamma    2.63201     
   model.ampl     3.86992e-20 
WARNING: parameter value model.ampl is at its minimum boundary 0.0
datasets       = (1,)
itermethodname = none
methodname     = levmar
statname       = cash
succeeded      = True
parnames       = ('model.gamma', 'model.ampl')
parvals        = (2.6320092037472103, 3.8699228650782024e-20)
statval        = -156.77839882881938
istatval       = -48.87368976960537
dstatval       = 107.90470905921401
numpoints      = 27
dof            = 25
qval           = None
rstat          = None
message        = successful termination
nfev           = 202


"""
