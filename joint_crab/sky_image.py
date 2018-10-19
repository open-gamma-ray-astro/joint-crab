"""Make a sky image for each dataset."""
import logging
from pathlib import Path
from astropy.table import Table
from gammapy.maps import WcsNDMap
from .conf import config

log = logging.getLogger(__name__)


def main():
    make_exclusion_mask()

    for name in config.all_datasets:
        sky_image_compute(name)
        sky_image_plot(name)


def make_map():
    # TODO: make source_pos a global setting, instead of per-dataset?
    skydir = config.get_dataset("hess").source_pos
    return WcsNDMap.create(skydir=skydir, proj="TAN", width=8, binsz=0.05)


def make_exclusion_mask():
    m = make_map()
    m.data += 1  # zero is excluded, one is not excluded
    path = Path("results/skyimage")
    path.mkdir(parents=True, exist_ok=True)
    filename = "results/skyimage/exclusion_mask.fits.gz"
    log.info(f"Writing {filename}")
    m.write(filename, overwrite=True)


def sky_image_compute(which):
    m = make_map()

    if which == "fermi":
        table = Table.read("data/fermi/events.fits.gz", hdu="EVENTS")
        m.fill_by_coord((table["RA"], table["DEC"]))
    else:
        dataset = config.get_dataset(which)
        datastore = dataset.get_DataStore()
        for obs_id in dataset.obs_ids:
            table = datastore.obs(obs_id).events.table
            m.fill_by_coord((table["RA"], table["DEC"]))

    path = Path(f"results/skyimage/{which}.fits.gz")
    path.parent.mkdir(exist_ok=True)
    log.info(f"Writing {path}")
    m.write(path, overwrite=True)


def sky_image_plot(which):
    import matplotlib.pyplot as plt

    # TODO: smooth and plot counts image
    # TODO: overplot on and off regions

    filename = f"results/skyimage/{which}.fits.gz"
    map = WcsNDMap.read(filename)
    fig, ax, im = map.plot(stretch="sqrt")

    filename = f"results/skyimage/{which}.png"
    log.info(f"Writing {filename}")
    fig.savefig(filename)
