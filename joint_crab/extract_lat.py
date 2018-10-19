"""Extract 1D spectrum information for Fermi-LAT"""
import logging
import numpy as np
import astropy.units as u
from astropy.coordinates import Angle
from astropy.table import Table
from regions import CircleSkyRegion
from gammapy.irf import EffectiveAreaTable, EnergyDispersion
from gammapy.spectrum import PHACountsSpectrum, SpectrumObservation
from gammapy.background import ring_background_estimate
from gammapy.utils.scripts import make_path

__all__ = [
    "extract_aeff",
    "extract_psf_containment",
    "fermi_ring_background_extract",
    "extract_edisp",
    "SpectrumExtractionFermi1D",
]

log = logging.getLogger(__name__)


def extract_aeff(exposure, target_position, energy):
    """Function to extract the effective area of Fermi-LAT starting from a
    file containing the exposure.

    We rely on the following assumption:

    1) the exposure is flat in the ROI we are considering, so we assume the
    value at the source position is valid all over the ROI

    2) since the exposure is the convolution of the effective area and time
    what we do is simply to store the exposure as it is [cm2 s] and then
    put a mock livetime of 1 s in the PHACountsSpectrum definition later on

    Parameters
    ----------
    exposure : `~gammapy.maps.HpxMapND`
        Fermi-LAT exposure
    target_position : `~astropy.coordinates.SkyCoord`
        position of the source of interest
    energy : `~astropy.units.Quantity`
        energy binning to calculate the effective area

    Returns
    -------
    `~gammapy.irf.EffectiveAreaTable`
    table containing the effective area
    """
    energy_log_ctr = np.sqrt(energy[1:] * energy[:-1])
    # the exposure is defined pixel-wise, fetch the value at source position
    lon = target_position.galactic.l
    lat = target_position.galactic.b
    expo_values = exposure.get_by_coord((lon.value, lat.value, energy_log_ctr.value))
    # generate an ARF table for the effective area
    # according to the format defined in
    # https://heasarc.gsfc.nasa.gov/docs/heasarc/caldb/docs/memos/cal_gen_92_002/cal_gen_92_002.html#Sec:ARF-format
    # we will create an astropy table that can be read by EffectiveAreaTable
    # via from_table method

    # we need low and high bin edges
    table = Table(
        [energy[:-1], energy[1:], expo_values],
        names=("ENERG_LO", "ENERG_HI", "SPECRESP"),
        dtype=("float64", "float64", "float64"),
        meta={"name": "Fermi-LAT exposure"},
    )

    # declare units
    table["ENERG_LO"].unit = str(energy.unit)
    table["ENERG_HI"].unit = str(energy.unit)
    table["SPECRESP"].unit = "cm2"

    # add mandatory keywords in header
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
    """Function the extract the containment fraction
    from the PSF of Fermi-LAT.

    note : we are assuming the PSF was calculated at
    the source postion

    Parameters
    ----------
    psf : `~gammapy.irf.EnergyDependentTablePSF`
        Fermi-LAT PSF
    on_radius : `~astropy.units.Quantity`
        size in degree of the region from which we take the ON counts
    energy : `~astropy.units.Quantity`
        energy binning to calculate the the psf
    Returns
    -------
    `~np.ndarray`
    containtment fraction in each energy bin
    """
    energy_log_ctr = np.sqrt(energy[:-1] * energy[1:])
    # integrate up to the on_radius radius
    containment_factor = np.asarray(
        [psf.integral(_, rad_min=0 * u.deg, rad_max=on_radius) for _ in energy_log_ctr]
    )
    return containment_factor


def fermi_ring_background_extract(
    events, target_position, on_radius, inner_radius=1 * u.deg, outer_radius=2 * u.deg
):
    """Function to extract the on and off event list through a ring background
    estimation

    Parameters
    ----------
    events : `~gammapy.data.EventList`
        Fermi-LAT event list
    on_radius : `~astropy.units.Quantity`
        size in degree of the region from which we take the ON counts

    Returns
    -------
    `~gammapy.background.BackgroundEstimate`
    """
    # we start from with a default ROI of 3 degree, should be enough for our purposes
    roi_radius = Angle("3 deg")
    roi_region = CircleSkyRegion(center=target_position, radius=roi_radius)
    evtlist = events.select_circular_region(roi_region)

    # estimate background
    ring_background = ring_background_estimate(
        pos=target_position,
        on_radius=on_radius,
        inner_radius=inner_radius,
        outer_radius=outer_radius,
        events=evtlist,
    )

    return ring_background


def extract_edisp(energy):
    """Dummy migration matrix for the fermi filename

    Its presence is needed for the sherpa fit, otherwise will trigger the issue
    commented in
    https://github.com/sherpa/sherpa/blob/7f9d31714eccfc7851fe2c439e5f186781109211/sherpa/astro/instrument.py#L556

    Parameters
    ----------
    energy : `~astropy.units.Quantity`
        energy binning to calculate the the migration matrix,
        since it's a mock energy dispersion we assume same binning
        for `e_true` and `e_reco`

    Returns
    -------
    `~gammapy.irf.EnergyDispersion`
    energy dispersion ontained assuming bias 0 and 0.05 resolution
    """
    edisp = EnergyDispersion.from_gauss(
        e_true=energy, e_reco=energy, sigma=0.05, bias=0
    )
    return edisp


class SpectrumExtractionFermi1D(object):
    """Creating input data to 1D spectrum fitting using Fermi-LAT data

    Class based on
    http://docs.gammapy.org/en/latest/_modules/gammapy/spectrum/extract.html#SpectrumExtraction

    Parameters
    ----------
    events : `~gammapy.data.EventList`
        Fermi-LAT event list
    exposure : `~gammapy.maps.HpxMapND`
        Fermi-LAT exposure
    psf : `~gammapy.irf.EnergyDependentTablePSF`
        Fermi-LAT PSF
    bkg_estimate : `~gammapy.background.BackgroundEstimate`
        Background estimate
    target_position : `~astropy.coordinates.SkyCoord`
        postion of the source of interest
    on_radius : `~astropy.units.Quantity`
        size in degree of the region from which we take the ON counts
        this is important both in ON counts calculation than in PSF
    containment correction : bool
        apply the psf correction when extracting the effective area
    energy : `~astropy.units.Quantity`
        energy binning to calculate the the irfs and extract the counts

    Returns
    -------
    `~gammapy.spectrum.SpectrumObservation`
    SpectrumObservation object
    """

    # TODO: this is probably not a useful default for most people -> remove?
    DEFAULT_ENERGY = np.logspace(4, np.log10(2e6), 21) * u.MeV
    """default energy axis to calculate irfs and extract counts, energy
    range of the 3FHL divided in 20 bins"""

    def __init__(
        self,
        events,
        exposure,
        psf,
        bkg_estimate,
        target_position,
        on_radius,
        energy=None,
        containment_correction=True,
    ):
        self.events = events
        self.exposure = exposure
        self.psf = psf
        self.bkg_estimate = bkg_estimate
        self.target_position = target_position
        self.on_radius = on_radius
        self.containment_correction = containment_correction
        # TODO: at the moment, the energy axis for exposure doesn't have unit information.
        # So needs to be MeV as used in Fermi ST
        self.energy = energy.to("MeV") or self.DEFAULT_ENERGY

    def make_empty_vectors(self, bkg_estimate):
        """Create empty vectors.

        This method sets up the energy binning of the PHACountsSpectrum

        Parameters
        ----------
        bkg : `~gammapy.background.BackgroundEstimate`
            Background estimate
        """
        log.info("Update observation meta info")
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
        self._on_vector.livetime = 1. * u.s
        self._off_vector.livetime = 1. * u.s

    def extract_counts(self, bkg_estimate):
        """Fill on and off vector for one observation.

        Parameters
        ----------
        bkg : `~gammapy.background.BackgroundEstimate`
            Background estimate
        """
        log.info("Fill events")
        self._on_vector.fill(bkg_estimate.on_events)
        self._off_vector.fill(bkg_estimate.off_events)

    def run(self):
        """Performs the spectrum extraction

        Returns
        -------
        spectrum_observation : `~gammapy.spectrum.SpectrumObservation`
            Spectrum observation
        """
        self.make_empty_vectors(self.bkg_estimate)
        self.extract_counts(self.bkg_estimate)

        aeff = extract_aeff(self.exposure, self.target_position, self.energy)

        if self.containment_correction:
            containment_factor = extract_psf_containment(
                self.psf, self.on_radius, self.energy
            )
            aeff.data.data *= containment_factor

        edisp = extract_edisp(self.energy)

        self._aeff = aeff
        self._edisp = edisp

        return SpectrumObservation(
            on_vector=self._on_vector,
            aeff=self._aeff,
            off_vector=self._off_vector,
            edisp=self._edisp,
        )

    def write(self, outdir, ogipdir="ogip_data", use_sherpa=False):
        """Write results to disk.

        Parameters
        ----------
        outdir : `~gammapy.extern.pathlib.Path`
            Output folder
        ogipdir : str, optional
            Folder name for OGIP data, default: 'ogip_data'
        use_sherpa : bool, optional
            Write Sherpa compliant files, default: False
        """
        outdir = make_path(outdir)
        log.info("Writing OGIP files to {}".format(outdir / ogipdir))
        outdir.mkdir(exist_ok=True, parents=True)
        self.observations.write(outdir / ogipdir, use_sherpa=use_sherpa)
