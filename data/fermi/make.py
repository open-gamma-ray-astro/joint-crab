"""Prepare Fermi-LAT data for this paper."""
import subprocess
from astropy.io import fits
from astropy.coordinates import SkyCoord, Angle


def download():
    url = "https://github.com/gammapy/gammapy-fermi-lat-data/raw/master/3fhl"
    subprocess.call(f"wget {url}/fermi_3fhl_psf_gc.fits.gz -O psf.fits.gz", shell=True)
    subprocess.call(
        f"wget {url}/fermi_3fhl_exposure_cube_hpx.fits.gz -O exposure_cube.fits.gz",
        shell=True,
    )
    subprocess.call(
        f"wget {url}/fermi_3fhl_events_selected.fits.gz -O events_allsky.fits.gz",
        shell=True,
    )


def event_select():
    filename = "events_allsky.fits.gz"
    print(f"Reading {filename}")
    hdu_list = fits.open(filename)

    events = hdu_list["EVENTS"].copy()
    print("len(events):", len(events.data))
    event_coord = SkyCoord(events.data["RA"], events.data["DEC"], unit="deg")

    target_coord = SkyCoord(83.633_212, 22.01446, unit="deg")  # Crab position
    sep = target_coord.separation(event_coord)
    fov = Angle(10, "deg")
    mask = sep < fov
    data = events.data[mask]
    events2 = fits.BinTableHDU(data=data, header=events.header)
    events2.add_datasum()
    events2.add_checksum()
    print("len(events2):", len(events2.data))
    hdu_list2 = fits.HDUList([fits.PrimaryHDU(), events2, hdu_list["GTI"]])

    filename = "events.fits.gz"
    print("Writing", filename)
    hdu_list2.writeto(filename, overwrite=True)


if __name__ == "__main__":
    download()
    event_select()
