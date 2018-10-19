"""Common utility functions."""
import logging
from pathlib import Path
import ruamel.yaml

log = logging.getLogger(__name__)


def write_yaml(data, path):
    """Helper function to write data to a YAML file."""
    path = Path(path)
    log.info("Writing {}".format(path))
    with path.open("w") as fh:
        ruamel.yaml.round_trip_dump(data, fh)


def load_yaml(path):
    """Helper function to load data from a YAML file."""
    path = Path(path)
    log.debug("Reading {}".format(path))
    with path.open() as fh:
        data = ruamel.yaml.round_trip_load(fh)
    return data
