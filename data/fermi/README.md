# Fermi-LAT data

This folder contains the Fermi-LAT data used for this paper.

- `events.fits.gz` has events above 10 GeV, within 10 deg of the Crab position.
- `exposure_cube.fits.gz` is an all-sky exposure cube in HEALPix format,
  at 18 energies, equally spaced in `log(energy)` from 10 GeV to 2 TeV.
- `psf.fits.gz` contains the point spread function
  at 17 energies, equally spaced in `log(energy)` from 10 GeV to 2 TeV.

As you can see in `make.py`, these files were produced from the Fermi-LAT
dataset prepared at https://github.com/gammapy/gammapy-fermi-lat-data, matching
the 3FHL catalog analysis.

There, the Fermi Science Tools were used, specifically ``gtselect`` and
``gtmktime`` to prepare the events, and ``gtexpcube2`` to compute the exposure
cube.

The PSF was computed using ``gtpsf`` there at the Galactic center postion. Given
that the variation of the Fermi-LAT PSF with sky position is very small above 10
GeV, we just used that PSF for the analysis here at the Crab position.

A more detailed description of the Fermi-LAT data, and the analysis is in the
paper.
