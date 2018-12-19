"""Extract 1D spectrum information for Fermi-LAT"""
import logging
import numpy as np
import astropy.units as u
from astropy.table import Table
from gammapy.irf import EffectiveAreaTable, EnergyDispersion
from gammapy.spectrum import PHACountsSpectrum, SpectrumObservation

log = logging.getLogger(__name__)


class SpectrumExtractionFermi1D(object):
    def __init__(
        self, events, exposure, psf, bkg_estimate, target_position, on_radius, energy
    ):
        self.events = events
        self.exposure = exposure
        self.psf = psf
        self.bkg_estimate = bkg_estimate
        self.target_position = target_position
        self.on_radius = on_radius
        self.energy = energy.to("MeV")

    def make_empty_vectors(self, bkg_estimate):
        self._on_vector = PHACountsSpectrum(
            energy_lo=self.energy[:-1],
            energy_hi=self.energy[1:],
            backscal=bkg_estimate.a_on,
            obs_id=0,
        )

        self._off_vector = self._on_vector.copy()
        self._off_vector.is_bkg = True
        self._off_vector.backscal = bkg_estimate.a_off
        # here we set the livetime of 1s, because we are actually storing an
        # exposure rather than an effective area
        self._on_vector.livetime = 1.0 * u.s
        self._off_vector.livetime = 1.0 * u.s

    def extract_counts(self, bkg_estimate):
        self._on_vector.fill(bkg_estimate.on_events)
        self._off_vector.fill(bkg_estimate.off_events)

    def run(self):
        self.make_empty_vectors(self.bkg_estimate)
        self.extract_counts(self.bkg_estimate)

        aeff = extract_aeff(self.exposure, self.target_position, self.energy)

        containment_factor = extract_psf_containment(
            self.psf, self.on_radius, self.energy
        )
        aeff.data.data *= containment_factor

        # Not the real Fermi-LAT EDISP
        # Use 5% energy resolution as approximation
        edisp = EnergyDispersion.from_gauss(e_true=self.energy, e_reco=self.energy, sigma=0.05, bias=0)

        return SpectrumObservation(
            on_vector=self._on_vector,
            aeff=aeff,
            off_vector=self._off_vector,
            edisp=edisp,
        )


def extract_aeff(exposure, target_position, energy):
    energy_log_ctr = np.sqrt(energy[1:] * energy[:-1])
    lon = target_position.galactic.l
    lat = target_position.galactic.b
    expo_values = exposure.get_by_coord((lon.value, lat.value, energy_log_ctr.value))

    table = Table(
        [energy[:-1], energy[1:], expo_values],
        names=("ENERG_LO", "ENERG_HI", "SPECRESP"),
        dtype=("float64", "float64", "float64"),
        meta={"name": "Fermi-LAT exposure"},
    )

    table["ENERG_LO"].unit = str(energy.unit)
    table["ENERG_HI"].unit = str(energy.unit)
    table["SPECRESP"].unit = "cm2"

    table.meta["EXTNAME"] = "SPECRESP"
    table.meta["TELESCOP"] = "Fermi"
    table.meta["INSTRUME"] = "LAT"
    table.meta["EXPOSURE"] = "1"
    table.meta["FILTER"] = ""
    table.meta["HDUCLASS"] = "OGIP"
    table.meta["HDUCLAS1"] = "RESPONSE"
    table.meta["HDUCLAS2"] = "SPECRESP"
    table.meta["HDUVERS"] = "1.1.0"

    return EffectiveAreaTable.from_table(table)


def extract_psf_containment(psf, on_radius, energy):
    energy_log_ctr = np.sqrt(energy[:-1] * energy[1:])
    containment_factor = np.asarray(
        [psf.integral(_, rad_min="0 deg", rad_max=on_radius) for _ in energy_log_ctr]
    )
    return containment_factor
