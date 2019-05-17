# Towards open and reproducible multi-instrument analysis in gamma-ray astronomy
[![DOI](https://zenodo.org/badge/146204837.svg)](https://zenodo.org/badge/latestdoi/146204837)

[![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/open-gamma-ray-astro/joint-crab/master?urlpath=lab/tree/joint-crab)

This repository contains material related to the publication [Towards multi-instrument and reproducible gamma-ray analysis](https://www.aanda.org/articles/aa/full_html/2019/05/aa34938-18/aa34938-18.html) appeared in A&amp;A 625, A10 (2019)

- DOI: [10.1051/0004-6361/201834938](https://doi.org/10.1051/0004-6361/201834938)
- ADS: [2019A%26A...625A..10N](https://ui.adsabs.harvard.edu/abs/2019A%26A...625A..10N/abstract)
- Cite: [Bibtex entry](https://ui.adsabs.harvard.edu/abs/2019A%26A...625A..10N/exportcitation)

This `joint-crab` bundle allows to reproduce the published results by its installation and execution in the local desktop, using a Docker image or via the Binder cloud service.

##  Contents

####  Datasets

You may find the datasets used as inputs for the analysis scripts in the [folder data](data).
They are grouped in subfolders accounting for their origin. They consist of Crab observations from the Fermi-LAT gamma-ray space telescope, as well as from four ground-based gamma-ray telescopes (MAGIC, VERITAS, FACT, H.E.S.S.).

####  Results

The results produced in the analysis are placed in the [folder results](results).
The tables and figures published in the paper have been prepared using these results.
You may reproduce the results following the instructions given in the [analysis section](analysis.md).

* Table 1: [results/summary/data.md](results/summary/data.md) | [results/summary/data.tex](results/summary/data.tex)
* Table 2: [results/summary/results.md](results/summary/results.md) | [results/summary/results.tex](results/summary/results.tex)
* Fig. 1: [results/figures/counts_spectra.png](results/figures/counts_spectra.png) | [results/figures/counts_spectra.pdf](results/figures/counts_spectra.pdf)
* Fig. 2: [results/figures/crab_sed_fit.png](results/figures/crab_sed_fit.png) | [results/figures/crab_sed_fit.pdf](results/figures/crab_sed_fit.pdf)
* Fig. 3: [results/figures/contours.png](results/figures/contours.png) | [results/figures/contours.pdf](results/figures/contours.pdf)
* Fig. 4: [results/figures/errorband_sed_veritas.png](results/figures/errorband_sed_veritas.png) | [results/figures/errorband_sed_veritas.pdf](results/figures/errorband_sed_veritas.pdf)
* Fig. 5: [results/figures/contours_systematics.png](results/figures/contours_systematics.png) | [results/figures/contours_systematics.pdf](results/figures/contours_systematics.pdf)

####  Scripts

The scripts needed to perform the analysis and reproduce the results are placed in the [joint_crab](joint_crab) folder.

####  Jupyter notebooks

We also provide several Jupyter notebooks to explore the data or results, or contain extra analyses in addition to the Python scripts.

* [1_data.ipynb](1_data.ipynb) - overview of input data files, as well as sky maps and spectra, and how to access them from Python
* [2_results.ipynb](2_results.ipynb) - overview of results file and how to access them from Python
* [3_systematics.ipynb](3_systematics.ipynb) - a detailed description of the systematic error likelihood fit
* [4_naima.ipynb](4_naima.ipynb) - SED model fit using inverse Compton radiation from electrons using Naima
* [5_crab_pulsar_nebula_sed.ipynb](5_crab_pulsar_nebula_sed.ipynb) - A check that the Crab pulsar emission is faint above 30 GeV

You can read these notebooks on Github as static documents, or execute them on your local machine or the Binder cloud infrastrcuture (see instructions below).

##  Executing the `joint-crab` bundle

You can execute the provided scripts in a specific `joint-crab` *conda environment*. This will allow you to reproduce the results published in the paper on Linux or MacOS (not on Windows, because we do use ``healpy`` to access the Fermi-LAT exposure map in HEALPix format, and ``healpy`` doesn't support Windows).

####  Set-up the environment

As a first approach to assure the reproducibility of the results we provide a [conda configuration file](binder/environment.yml) to build a virtual environment with pinned versions for the software dependencies.

Hence, as a requisite you need [Anaconda](https://www.anaconda.com/download/) or
[miniconda](https://conda.io/miniconda.html) software installed in your desktop. Once you have installed this software, and downloaded the content of this `joint-crab` repository, you can type in the terminal the following commands at the top-level of the `joint-crab` folder:

    $ conda env create -f binder/environment.yml
    $ conda activate joint-crab

####  Reproduce the results
* [How to reproduce results and figures with `make.py`](analysis.md)

####  Check out the notebooks

If you have set-up the environment as it is described above, you can open and run the complemetary notebooks that we provide by typing in the terminal of your desktop:

```
jupyter lab
```

##  Binder space

Alternatively you can also open an run the notebooks in the Binder cloud infrastructure. You can also open a terminal tab in the Binder space and run the analysis as it is decribed in the [Reproduce the results](analysis.md) section.

[![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/open-gamma-ray-astro/joint-crab/master?urlpath=lab/tree/joint-crab)

## Docker
Since the availability of all the external software dependencies is not assured in the future, the *conda virtual environment* approach does not guarantee a mid-term preservation of the reproducibility. It is mainly because of this reason that we also provide a `joint-crab` docker container that may be accessed from the public [Gammapy DockerHub repository](https://hub.docker.com/u/gammapy/dashboard/).
* [How to reproduce results with Docker](docker.md)


## License

The code and content in this repository is shared under the [BSD-3-Clause license](LICENSE) (read more at [OSI](https://opensource.org/licenses/BSD-3-Clause)).
