"""Prepare H.E.S.S. data for this paper.

The four Crab runs here are the Crab data subset
of the first public H.E.S.S. text data release.

This script starts with that data release, and creates
index files for this subset, and copies the relevant
data files.
"""
import subprocess
from pathlib import Path
import shutil
import numpy as np
from astropy.table import Table
from gammapy.data import DataStore

subprocess.call("gammapy download datasets --src hess-dl3-dr1 --out .", shell=True)

path_in = Path("hess-dl3-dr1")
obs_ids = [23523, 23526, 23559, 23592]


def make_obs_index():
    table = Table.read(path_in / "obs-index.fits.gz")
    mask = np.array([_ in obs_ids for _ in table["OBS_ID"]])
    table = table[mask]

    filename = "obs-index.fits.gz"
    print(f"Writing {filename}")
    table.write(filename, overwrite=True)


def make_hdu_index():
    table = Table.read(path_in / "hdu-index.fits.gz")
    mask = np.array([_ in obs_ids for _ in table["OBS_ID"]])
    table = table[mask]

    filename = "hdu-index.fits.gz"
    print(f"Writing {filename}")
    table.write(filename, overwrite=True)


def copy_data():
    Path("data").mkdir(exist_ok=True)
    ds_in = DataStore.from_dir(path_in)
    ds_out = DataStore.from_dir(".")

    for obs_id in obs_ids:
        loc_in = ds_in.obs(obs_id).location(hdu_type="events")
        loc_out = ds_out.obs(obs_id).location(hdu_type="events")
        src = path_in / loc_in.file_dir / loc_in.file_name
        dst = Path(loc_out.file_dir) / loc_out.file_name
        print(f"cp {src} {dst}")
        shutil.copy(src, dst)


if __name__ == "__main__":
    make_obs_index()
    make_hdu_index()
    copy_data()
