"""Extract 1D spectra for IACT datasets."""
import logging
import astropy.units as u
from regions import CircleSkyRegion
from gammapy.data import EventList, DataStore
from gammapy.maps import HpxNDMap
from gammapy.irf import EnergyDependentTablePSF
from gammapy.background import (
    ReflectedRegionsBackgroundEstimator,
    ring_background_estimate,
)
from gammapy.spectrum import SpectrumExtraction
from .extract_lat import SpectrumExtractionFermi1D
from .conf import config

log = logging.getLogger(__name__)


def extract_spectra(which):
    if which in {"fermi", "all"}:
        dataset = config.datasets["fermi"]
        extract_spectra_fermi(config.source_pos, dataset.on_radius)

    if which in {"magic", "all"}:
        dataset = config.datasets["magic"]
        extract_spectra_iact(dataset)

    if which in {"hess", "all"}:
        dataset = config.datasets["hess"]
        extract_spectra_iact(dataset)

    if which in {"fact", "all"}:
        dataset = config.datasets["fact"]
        extract_spectra_iact(dataset)

    if which in {"veritas", "all"}:
        dataset = config.datasets["veritas"]
        extract_spectra_iact(dataset)


def extract_spectra_fermi(target_position, on_radius):
    """Extract 1d spectra for Fermi-LAT"""
    log.info("Extracting 1d spectra for Fermi-LAT")
    events = EventList.read("data/fermi/events.fits.gz")
    exposure = HpxNDMap.read("data/fermi/exposure_cube.fits.gz")
    psf = EnergyDependentTablePSF.read("data/fermi/psf.fits.gz")

    valid_range = (config.energy_bins >= 30 * u.GeV) * (config.energy_bins <= 2 * u.TeV)
    energy = config.energy_bins[valid_range]

    bkg_estimate = ring_background_estimate(
        pos=target_position,
        on_radius=on_radius,
        inner_radius=1 * u.deg,
        outer_radius=2 * u.deg,
        events=events,
    )

    extract = SpectrumExtractionFermi1D(
        events=events,
        exposure=exposure,
        psf=psf,
        bkg_estimate=bkg_estimate,
        target_position=target_position,
        on_radius=on_radius,
        energy=energy,
    )
    obs = extract.run()

    path = "results/spectra/fermi"
    log.info(f"Writing to {path}")
    obs.write(path, use_sherpa=True, overwrite=True)


def extract_spectra_iact(dataset):
    """Extract 1d spectra for IACT dataset"""
    log.info(f"Extracting 1d spectra for {dataset.name} dataset")
    datastore = DataStore.from_dir(f"data/{dataset.name}")
    observations = datastore.get_observations(dataset.obs_ids)

    on_region = CircleSkyRegion(center=config.source_pos, radius=dataset.on_radius)

    exclusion_mask = config.get_exclusion_mask()

    bkg_estimate = ReflectedRegionsBackgroundEstimator(
        observations=observations, on_region=on_region, exclusion_mask=exclusion_mask
    )
    bkg_estimate.run()

    extract = SpectrumExtraction(
        observations=observations,
        bkg_estimate=bkg_estimate.result,
        e_true=config.energy_bins,
        e_reco=config.energy_bins,
        containment_correction=dataset.containment_correction,
    )
    extract.run()

    path = f"results/spectra/{dataset.name}/"
    log.info(f"Writing to {path}")

    if dataset.name == "fact":
        # For FACT the IRFs are the same for all observations
        # So we only store a stacked spectrum and response
        # plus we add a LO_THRESHOLD keyword was missing
        obs = extract.spectrum_observations.stack()
        obs.lo_threshold = 0.4 * u.TeV
        # we are writing a single observation, as for Fermi
        obs.write(path, use_sherpa=True, overwrite=True)
    else:
        extract.write(path, ogipdir="", use_sherpa=True, overwrite=True)
