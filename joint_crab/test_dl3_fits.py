"""Check DL3 data files for spec compliance.

https://gamma-astro-data-formats.readthedocs.io/
"""
from astropy.io import fits
from astropy.units import Unit
from argparse import ArgumentParser
import logging
import yaml
import re
import json
from pathlib import Path

log = logging.getLogger("checkdl3")


def main():
    fmt = logging.Formatter("%(levelname)8s %(message)s")
    handler = logging.StreamHandler()
    handler.setFormatter(fmt)
    log.addHandler(handler)
    loglevels = {0: logging.WARNING, 1: logging.INFO, 2: logging.DEBUG}

    parser = ArgumentParser()
    parser.add_argument("inputfile")
    parser.add_argument("--verbose", "-v", action="count", default=0)
    parser.add_argument("--json", action="store_true")

    args = parser.parse_args()
    log.setLevel(loglevels.get(args.verbose, logging.DEBUG))
    errors = check_files([args.inputfile])

    if args.json:
        print(json.dumps(errors))


def check_files(filenames):
    path = Path(__file__).parent / "../config/gamma_dl3_specs.yaml"
    specs = yaml.safe_load(path.read_text())
    errors = {}
    for filename in filenames:
        errors.update(_check_file(filename, specs))
    return errors


def _check_file(filename, specs):
    errors = {}
    log.info(f"Opening {filename}")
    with fits.open(filename) as f:
        log.info("File contains the following hdus:")
        for hdu in f:
            log.info("* " + hdu.name)

        for hdu in f:
            if hdu.name == "PRIMARY":
                continue

            log.info(f'Checking HDU "{hdu.name}"')
            errors[hdu.name] = check_hdu(hdu, specs["extensions"])

    return errors


def check_hdu(hdu, specs):
    header = hdu.header
    name = hdu.name

    errors = {
        "columns": {
            "missing": [],
            "wrong_dtype": [],
            "wrong_dim": [],
            "wrong_unit": [],
            "additional": [],
        },
        "header": {"missing": []},
    }

    try:
        hduclas1 = header["HDUCLAS1"]
    except KeyError:
        log.error(f'HDU "{name}" is missing required keyword "HDUCLAS1"')
        errors["header"]["missing"].append("HDUCLAS1")

        log.info("Trying EXTNAME to identify hduclass")
        if name in specs:
            hduclas1 = name
        else:
            log.warning(f'Unknown extension "{name}"')
            return errors

    ext_specs = specs[hduclas1]

    # check if the HDUCLASS has subclasses and get the right specs
    # this is needed for e.g. RESPONSE, where subclasses are EFF_AREA, EDISP,
    # etc.
    i = 2
    while "subclasses" in ext_specs:
        try:
            hduclass = header[f"HDUCLAS{i}"]
        except KeyError:
            log.error(f'HDU "{name}" is missing required keyword "HDUCLAS{i}"')
            errors["header"]["missing"].append("HDUCLAS{i}")
            return
        try:
            ext_specs = ext_specs["subclasses"][hduclass]
        except KeyError:
            log.error(f'Unknown HDUCLAS{i}: "{hduclass}"')
            return errors
        i += 1

    # Check that all mandatory columns exist and have correct type and unit
    mandatory = ext_specs["columns"]["mandatory"]
    colnames = [col.name for col in hdu.columns]
    for name, spec in mandatory.items():
        if name not in colnames:
            errors["columns"]["missing"].append(name)
            log.error(f"Extension {hdu.name} is missing mandatory column {name}")
            continue

        fmt = hdu.columns[name].format
        dtype = re.match("(?:\d*)([LXBIJKAEDCMPQ])", fmt).groups()[0]
        if dtype not in spec["type"]:
            errors["columns"]["wrong_dtype"].append(name)
            log.error(f"Column {name} in Extension {hdu.name} has wrong dtype: {dtype}")

        dim = hdu.columns[name].dim
        ndim = len(dim.split(",")) if dim else None
        if spec.get("ndim") != ndim:
            errors["columns"]["wrong_dim"].append(name)
            log.error(f"Column {name} in Extension {hdu.name} has wrong dim: {ndim}")

        unit = hdu.columns[name].unit or ""
        unit_spec = spec["unit"] or ""
        if unit != Unit(unit_spec):
            errors["columns"]["wrong_unit"].append(name)
            log.error(f"Column {name} in Extension {hdu.name} has wrong unit: {unit!r}")

    optional = ext_specs["columns"].get("optional", {})
    additional = sorted(set(colnames) - set(mandatory) - set(optional))
    if additional:
        for name in additional:
            log.warning(f'HDU "{hdu.name}" contains unknown column: {name}')
        errors["columns"]["additional"].extend(additional)

    # check that all
    common_mandatory_header = ["HDUCLASS", "HDUDOC", "HDUVERS"]
    mandatory = ext_specs["header"]["mandatory"] + common_mandatory_header
    for kwd in mandatory:
        if kwd not in hdu.header:
            errors["header"]["missing"].append(kwd)
            log.error(f"Header-Keyword {kwd} in Extension {hdu.name} missing")

    return errors


if __name__ == "__main__":
    main()
