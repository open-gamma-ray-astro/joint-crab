"""Gather provenance information.

The goal is to record information that helps to make
the analysis reproducible, e.g. on software versions.
"""
import platform
import sys
import getpass
import hashlib
from pathlib import Path
import numpy as np
import astropy
import gammapy


def get_provenance():
    """Compute provenance info about software and data used."""
    na = "not available"
    data = {}

    data["env"] = {}
    data["env"]["user"] = getpass.getuser()
    data["env"]["machine"] = platform.machine()
    data["env"]["system"] = platform.system()

    data["software"] = {}
    data["software"]["python_executable"] = sys.executable
    data["software"]["python_version"] = platform.python_version()
    data["software"]["numpy"] = np.__version__
    try:
        import scipy

        scipy_version = scipy.__version__
    except ImportError:
        scipy_version = na
    data["software"]["scipy"] = scipy_version
    data["software"]["astropy"] = astropy.__version__
    try:
        import sherpa

        sherpa_version = sherpa.__version__
    except ImportError:
        sherpa_version = na
    data["software"]["sherpa"] = sherpa_version
    data["software"]["gammapy"] = gammapy.__version__

    # Store MD5 checksum for all input files
    data["data"] = []
    for path in get_input_files():
        md5 = hashlib.md5(Path(path).read_bytes()).hexdigest()
        entry = dict(file=str(path), md5=md5)
        data["data"].append(entry)

    return data


def get_input_files():
    """Make a list of all input files from the `data` folder."""
    paths = Path("data").glob("**/*")
    paths = [_ for _ in paths if _.is_file()]

    return paths
