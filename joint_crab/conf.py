"""Configuration and data access helper class"""
import astropy.units as u
from astropy.coordinates import SkyCoord
from gammapy.data import DataStore
from gammapy.spectrum import SpectrumObservationList
from gammapy.maps import Map
from pathlib import Path
import glob
from .utils import load_yaml


class Config:
    """class containing the configutations for the analysis, for now it contains only the
    list of Dataset objects (see below)"""

    def __init__(self):
        self.repo_path = Path(__file__).parent.parent
        self.all_datasets = ("fermi", "magic", "hess", "fact", "veritas")
        self.datasets = self.load_datasets_from_yaml()

    def load_datasets_from_yaml(self):
        """load a list of particular datasets"""
        datasets = []
        config_yaml = load_yaml(f"{self.repo_path}/config/datasets.yaml")
        for _dict in config_yaml["datasets"]:
            # load the single dataset with all his info
            datasets.append(Dataset.from_dict(_dict))
        return datasets

    def get_dataset(self, name):
        """get a particular dataset by name {'fermi', 'magic', 'hess', 'fact', 'veritas', 'joint'}"""
        dataset_dict = {
            "fermi": self.datasets[0],
            "magic": self.datasets[1],
            "hess": self.datasets[2],
            "fact": self.datasets[3],
            "veritas": self.datasets[4],
            "joint": self.datasets[5],
        }
        if name in dataset_dict.keys():
            return dataset_dict[name]
        else:
            raise ValueError(f"Unknown dataset: {name}")

    def get_exclusion_mask(self):
        return Map.read(f"{self.repo_path}/results/skyimage/exclusion_mask.fits.gz")


class Dataset(object):
    """Dataset is a class holding the property of a each dataset:
    obs IDs, on region size, energy range to be used in the fit
    It will have methods to return the proper gammapy and sherpa objects needed for the analysis
    (gammapy DataStore and SpectrumObservationList, sherpa datastack)

    Parameters:
    -----------
    name : string
        string with source name, check that the method
    obs_ids : list
        list of the observations IDs
    source_pos : `~astropy.coordinates.SkyCoord`
        coordinates of the source (Crab for all the datasets in this case)
    on_region : `~astropy.units.Quantity`
        size of the on region to be used for the signal extraction
    containment_correction : boolean
        use the containment correction or not
    e_min, e_max : `~astropy.units.Quantity`
        energies defining the range to be fitted
    bins_decade : int
        number of bins per decade for the counts plot
    """

    def __init__(
        self,
        name,
        obs_ids,
        source_pos,
        on_radius,
        containment_correction,
        energy_range,
        bins_decade,
    ):
        self.name = name
        self.obs_ids = obs_ids
        self.source_pos = source_pos
        self.on_radius = on_radius
        self.containment_correction = containment_correction
        self.energy_range = energy_range
        self.bins_decade = bins_decade
        # directories to read / store results per dataset
        self.main_repo_path = Path(__file__).parent.parent
        self.data_path = f"{self.main_repo_path}/data/{self.name}"
        self.spectra_path = f"{self.main_repo_path}/results/spectra/{self.name}"

    @classmethod
    def from_dict(cls, _dict):
        """function to initialize the class starting from a dictionary with keys
        dictionary as defined per each source in the config.yaml
        """
        name = _dict["name"]
        obs_ids = _dict["obs_ids"]
        ra = _dict["source_pos"]["ra"]
        dec = _dict["source_pos"]["dec"]
        source_pos = SkyCoord(ra=ra, dec=dec)
        containment_correction = _dict["containment_correction"]
        on_radius = u.Quantity(_dict["on_radius"])
        e_min = float(_dict["energy_range"]["min"]) * u.Unit(
            _dict["energy_range"]["unit"]
        )
        e_max = float(_dict["energy_range"]["max"]) * u.Unit(
            _dict["energy_range"]["unit"]
        )
        # better for energy_range to be a quantity list than a list of quantity
        energy_range = [e_min.value, e_max.value] * e_min.unit
        bins_decade = int(_dict["energy_range"]["bins_decade"])

        return cls(
            name,
            obs_ids,
            source_pos,
            on_radius,
            containment_correction,
            energy_range,
            bins_decade,
        )

    def get_DataStore(self):
        """get gammapy DataStore object"""
        if self.name in ["fermi", "joint"]:
            return None
        elif self.name in ["magic", "hess", "fact", "veritas"]:
            return DataStore.from_dir(self.data_path)
        else:
            raise ValueError(f"No datastore for: {self.name}")

    def get_SpectrumObservationList(self):
        """load the OGIP files and return a SpectrumObservationList
        SpectrumObservationList has already a method to read from a directory
        http://docs.gammapy.org/dev/api/gammapy.spectrum.SpectrumObservationList.html
        """
        if self.name == "joint":
            spec_obs_list = SpectrumObservationList()
            # extend the list adding all the other SpectrumObservationList
            for name in {"fermi", "magic", "hess", "fact", "veritas"}:
                spec_obs = SpectrumObservationList.read(
                    f"{self.main_repo_path}/results/spectra/{name}"
                )
                spec_obs_list.extend(spec_obs)
        else:
            spec_obs_list = SpectrumObservationList.read(self.spectra_path)

        return spec_obs_list

    def get_sherpa_datastack(self):
        """return a sherpa datastack object

         http://cxc.harvard.edu/sherpa/ahelp/datastack.html
         depending on the instrument selected,
         set also the appropriate fit range
         """
        import sherpa.astro.datastack as sh

        data = sh.DataStack()

        if self.name == "joint":
            pha_list = glob.glob(f"{self.main_repo_path}/results/spectra/*/pha_*.fits")
        else:
            pha_list = glob.glob(f"{self.spectra_path}/pha_*.fits")

        # load all the OGIP files in the DataStack object
        for pha in pha_list:
            sh.load_data(data, pha)

        # convert the fit range to keV and strip units
        fit_range = self.energy_range.to("keV").value
        data.notice(*fit_range)

        return data


config = Config()
