"""Extract 1D spectra for IACT datasets."""
import logging
import numpy as np
import astropy.units as u
from regions import CircleSkyRegion
from gammapy.data import EventList
from gammapy.maps import HpxNDMap
from gammapy.irf import EnergyDependentTablePSF
from gammapy.background import ReflectedRegionsBackgroundEstimator
from gammapy.spectrum import SpectrumExtraction
from .extract_lat import fermi_ring_background_extract, SpectrumExtractionFermi1D
from .conf import config

log = logging.getLogger(__name__)


def extract_spectra(which):
    if which in {"fermi", "all"}:
        dataset = config.get_dataset("fermi")
        extract_spectra_fermi(dataset.source_pos, dataset.on_radius)

    if which in {"magic", "all"}:
        dataset = config.get_dataset("magic")
        extract_spectra_IACT(dataset)

    if which in {"hess", "all"}:
        dataset = config.get_dataset("hess")
        extract_spectra_IACT(dataset)

    if which in {"fact", "all"}:
        dataset = config.get_dataset("fact")
        extract_spectra_IACT(dataset)

    if which in {"veritas", "all"}:
        dataset = config.get_dataset("veritas")
        extract_spectra_IACT(dataset)


def extract_spectra_fermi(target_position, on_radius):
    """Extract 1d spectra for Fermi-LAT"""
    log.info("Extracting 1d spectra for Fermi-LAT")
    events = EventList.read("data/fermi/events.fits.gz")
    exposure = HpxNDMap.read("data/fermi/exposure_cube.fits.gz")
    psf = EnergyDependentTablePSF.read("data/fermi/psf.fits.gz")

    emin, emax, dex = 0.03, 2, 0.1
    num = int(np.log10(emax / emin) / dex)
    energy = np.logspace(start=np.log10(emin), stop=np.log10(emax), num=num) * u.TeV

    bkg_estimate = fermi_ring_background_extract(events, target_position, on_radius)

    extract = SpectrumExtractionFermi1D(
        events=events,
        exposure=exposure,
        psf=psf,
        bkg_estimate=bkg_estimate,
        target_position=target_position,
        on_radius=on_radius,
        energy=energy,
        containment_correction=True,
    )
    obs = extract.run()

    path = f"{config.repo_path}/results/spectra/fermi"
    log.info(f"Writing to {path}")
    obs.write(path, use_sherpa=True, overwrite=True)


def extract_spectra_IACT(dataset):
    """Extract 1d spectra for IACT dataset"""
    log.info(f"Extracting 1d spectra for {dataset.name} dataset")
    # Dataset class has already the method to obtain the gammapy DataStore object
    datastore = dataset.get_DataStore()
    obs_list = datastore.obs_list(dataset.obs_ids)

    on_region = CircleSkyRegion(center=dataset.source_pos, radius=dataset.on_radius)

    exclusion_mask = config.get_exclusion_mask()

    bkg_estimate = ReflectedRegionsBackgroundEstimator(
        obs_list=obs_list, on_region=on_region, exclusion_mask=exclusion_mask
    )
    bkg_estimate.run()

    extract = SpectrumExtraction(
        obs_list=obs_list,
        bkg_estimate=bkg_estimate.result,
        containment_correction=dataset.containment_correction,
    )
    extract.run()

    path = f"{config.repo_path}/results/spectra/{dataset.name}/"
    log.info(f"Writing to {path}")

    if dataset.name == "fact":
        # For FACT the IRFs are the same for all observations
        # So we only store a stacked spectrum and response
        # plus we add a LO_THRESHOLD keyword was missing
        obs = extract.observations.stack()
        obs.lo_threshold = 0.4 * u.TeV
        # we are writing a single observation, as for Fermi
        obs.write(path, use_sherpa=True, overwrite=True)
    else:
        extract.write(path, ogipdir="", use_sherpa=True, overwrite=True)
