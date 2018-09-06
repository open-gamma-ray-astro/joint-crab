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
from . import models
from . import maps

# Analysis
from . import extract_lat
from . import extract_ogip_spectra
from . import fit_spec
from . import systematics
from . import errorbands

# Text summary
from . import summary_data
from . import summary_results

# plots
from . import plot_counts
from . import plot_seds
from . import plot_errorbands
from . import plot_contours
