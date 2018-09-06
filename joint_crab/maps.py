"""Make a sky image for each dataset."""
import logging
from pathlib import Path
from astropy.table import Table
from gammapy.data import DataStore
from gammapy.maps import WcsNDMap
from .conf import config

log = logging.getLogger(__name__)


def main():
    make_exclusion_mask()

    for name in config.all_datasets:
        sky_image_compute(name)
        sky_image_plot(name)


def make_map():
    return WcsNDMap.create(skydir=config.source_pos, proj="TAN", width=8, binsz=0.05)


# TODO: remove exclusion mask? Shouldn't be needed, no?
def make_exclusion_mask():
    m = make_map()
    m.data += 1  # zero is excluded, one is not excluded
    filename = "results/maps/exclusion_mask.fits.gz"
    log.info(f"Writing {filename}")
    m.write(filename, overwrite=True)


def sky_image_compute(which):
    m = make_map()

    if which == "fermi":
        table = Table.read("data/fermi/events.fits.gz", hdu="EVENTS")
        m.fill_by_coord((table["RA"], table["DEC"]))
    else:
        data_store = DataStore.from_dir(f"data/{which}")
        for obs_id in data_store.obs_table["OBS_ID"]:
            table = data_store.obs(obs_id).events.table
            m.fill_by_coord((table["RA"], table["DEC"]))

    path = Path(f"results/maps/{which}.fits.gz")
    log.info(f"Writing {path}")
    m.write(path, overwrite=True)


def sky_image_plot(which):
    import matplotlib.pyplot as plt

    # TODO: smooth and plot counts image
    # TODO: overplot on and off regions

    filename = f"results/maps/{which}.fits.gz"
    map = WcsNDMap.read(filename)
    fig, ax, im = map.plot(stretch="sqrt")

    filename = f"results/maps/{which}.png"
    log.info(f"Writing {filename}")
    fig.savefig(filename)
