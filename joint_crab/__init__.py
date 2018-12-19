"""Python package for the `joint-crab` paper.

This is just used to be able to organise the code
across multiple files, and to be able to import
helper functions and classes from `make.py` or
Jupyter notebooks.
"""
# Basics
from . import utils
from . import conf
from . import provenance
from . import maps

# Analysis
from . import extract_fermi
from . import extract_spectra
from . import fit_models
from . import fit_spectra
from . import fit_systematics
from . import fit_errorbands

# Text summary
from . import summary_data
from . import summary_results

# plots
from . import plot_counts
from . import plot_seds
from . import plot_errorbands
from . import plot_contours
