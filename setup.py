"""
This `setup.py` file let's us install the `joint_crab` package via:

    pip install .

 or in editable mode via

    pip install -e .

to make the code available for `import joint_crab`, e.g. from notebooks.

If you'd like to learn more about `setup.py` and `pip`:
- https://pip.pypa.io/
- https://python-packaging.readthedocs.io/
"""
from setuptools import setup

setup(
    name='join_crab',
    # We aren't using or incrementing this version number.
    # Only git version control is used for versioning.
    # This is a default or "dummy" value.
    version='1.0',
    description='Python package for the joint Crab paper',
    url='https://github.com/open-gamma-ray-astro/joint-crab',
    packages=['joint_crab'],
)
