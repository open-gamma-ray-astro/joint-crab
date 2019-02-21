"""Configuration and data access helper class"""
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
from gammapy.spectrum import SpectrumObservationList
from gammapy.maps import Map

DATASETS = [
    {
        "name": "fermi",
        "label": "Fermi-LAT",
        "obs_ids": [0],
        "on_radius": "0.3 deg",
        "containment_correction": True,
        "energy_range": {"min": "0.03 TeV", "max": "2 TeV"},
        "color": "#21ABCD",
    },
    {
        "name": "magic",
        "label": "MAGIC",
        "on_radius": "0.141 deg",
        "containment_correction": False,
        "energy_range": {"min": "0.08 TeV", "max": "30 TeV"},
        "color": "#FF9933",
    },
    {
        "name": "veritas",
        "label": "VERITAS",
        "on_radius": "0.1 deg",
        "containment_correction": False,
        "energy_range": {"min": "0.15 TeV", "max": "30 TeV"},
        "color": "#893F45",
    },
    {
        "name": "fact",
        "label": "FACT",
        "on_radius": "0.1732 deg",
        "containment_correction": False,
        "energy_range": {"min": "0.4 TeV", "max": "30 TeV"},
        "color": "#3EB489",
    },
    {
        "name": "hess",
        "label": "H.E.S.S.",
        "on_radius": "0.11 deg",
        "containment_correction": True,
        "energy_range": {"min": "0.66 TeV", "max": "30 TeV"},
        "color": "#002E63",
    },
    {
        "name": "joint",
        "label": "joint fit",
        "energy_range": {"min": "0.03 TeV", "max": "30 TeV"},
        "color": "crimson",
    },
]


class PlotConfig:
    fontsize = 15
    fontsize_contours = 18
    label_energy = r"$E\,/\,\mathrm{TeV}$"
    e2dnde_label = (
        r"$E^2 \cdot {\rm d}\phi/{\rm d}E\,/\,({\rm erg}\,{\rm cm}^{-2}\,{\rm s}^{-1})$"
    )


class Config:
    """Configuration options."""

    source_pos = SkyCoord("83d37m59.0988s", "22d00m52.2s")
    plot = PlotConfig()

    all_datasets = ["fermi", "magic", "veritas", "fact", "hess"]
    all_datasets_plus_joint = all_datasets + ["joint"]

    def __init__(self):
        self.datasets = {_["name"]: Dataset.from_dict(_) for _ in DATASETS}

        # energy bins for all the datasets 10 GeV to 100 TeV, 20 bins per decade
        self.energy_bins = np.logspace(np.log10(0.01), np.log10(100), 81) * u.TeV

    @staticmethod
    def get_exclusion_mask():
        return Map.read("results/maps/exclusion_mask.fits.gz")


class Dataset(object):
    """Information about each dataset."""

    def __init__(
        self,
        name,
        label,
        obs_ids,
        on_radius,
        containment_correction,
        energy_range,
        color,
    ):
        self.name = name
        self.label = label
        self.obs_ids = obs_ids
        self.on_radius = on_radius
        self.containment_correction = containment_correction
        self.energy_range = energy_range
        self.color = color

        # directories to read / store results per dataset
        self.data_path = f"data/{self.name}"
        self.spectra_path = f"results/spectra/{self.name}"

    @classmethod
    def from_dict(cls, d):
        """function to initialize the class starting from a dictionary with keys
        dictionary as defined per each source in the config.yaml
        """
        e_min = u.Quantity(d["energy_range"]["min"])
        e_max = u.Quantity(d["energy_range"]["max"])
        energy_range = [e_min.value, e_max.value] * e_min.unit

        if d["name"] in ["magic", "veritas", "fact", "hess"]:
            table = Table.read(f"data/{d['name']}/obs-index.fits.gz")
            obs_ids = np.array(table["OBS_ID"])
        else:
            obs_ids = None

        return cls(
            name=d["name"],
            label=d["label"],
            obs_ids=obs_ids,
            on_radius=u.Quantity(d.get("on_radius", "nan")),
            containment_correction=d.get("containment_correction", "nan"),
            energy_range=energy_range,
            color=d["color"],
        )

    def get_energy_mask(self, energy):
        return (self.energy_range[0] <= energy) & (energy < self.energy_range[1])

    def get_SpectrumObservationList(self):
        """load the OGIP files and return a SpectrumObservationList
        SpectrumObservationList has already a method to read from a directory
        http://docs.gammapy.org/dev/api/gammapy.spectrum.SpectrumObservationList.html
        """
        if self.name == "joint":
            spec_obs_list = SpectrumObservationList()
            # extend the list adding all the other SpectrumObservationList
            for name in {"fermi", "magic", "veritas", "fact", "hess"}:
                spec_obs = SpectrumObservationList.read(f"results/spectra/{name}")
                spec_obs_list.extend(spec_obs)
        else:
            spec_obs_list = SpectrumObservationList.read(self.spectra_path)

        return spec_obs_list


config = Config()
