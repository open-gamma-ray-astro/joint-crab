# Towards multi-instrument and reproducible gamma-ray analysis

[![DOI](https://zenodo.org/badge/146204837.svg)](https://zenodo.org/badge/latestdoi/146204837)

[![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/open-gamma-ray-astro/joint-crab/master?urlpath=lab/tree/joint-crab)

This repository contains material related to the paper *Towards multi-instrument and reproducible gamma-ray analysis*.

This `joint-crab` bundle allows to reproduce the published results by its installation and execution in the local desktop, using a Docker image or via the MyBinder cloud service.

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
* Fig. 3: [results/figures/errorband_sed_veritas.png](results/figures/errorband_sed_veritas.png) | [results/figures/errorband_sed_veritas.pdf](results/figures/errorband_sed_veritas.pdf)
* Fig. 4: [results/figures/contours.png](results/figures/contours.png) | [results/figures/contours.pdf](results/figures/contours.pdf)
* Fig. 5: [results/figures/contours_systematics.png](results/figures/contours_systematics.png) | [results/figures/contours_systematics.pdf](results/figures/contours_systematics.pdf)

####  Scripts

The scripts needed to perform the analysis and reproduce the results are placed in the [joint_crab](joint_crab) folder.

####  Jupyter notebooks

We also provide several complemetary jupyter notebooks at the root level. You can run these notebooks locally or in the MyBinder cloud infrastrcuture.

##  Executing the `joint-crab` bundle

You can execute the provided scripts in a specific `joint-crab` *conda environment*. This will allow you to reproduce the results published in the paper.

####  Set-up the environment

As a first approach to assure the reproducibility of the results we provide a [conda configuration file](binder/environment.yml) to build a virtual environment with pinned versions for the software dependencies.

Hence, as a requisite you need [Anaconda](https://www.anaconda.com/download/) or
[miniconda](https://conda.io/miniconda.html) - Python>=3.6 software installed in your desktop. Once you have installed this software, and downloaded the content of this `joint-crab` repository, you can type in the terminal the following commands at the top-level of the `joint-crab` folder:

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

Alternatively you can also open an run the notebooks in the MyBinder cloud infrastructure. You can also open a terminal tab in the MyBinder space and run the analysis as it is decribed in the [Reproduce the results](analysis.md) section.

[![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/open-gamma-ray-astro/joint-crab/master?urlpath=lab/tree/joint-crab)

## Docker
Since the availability of all the external software dependencies is not assured in the future, the *conda virtual environment* approach does not guarantee a mid-term preservation of the reproducibility. It is mainly because of this reason that we also provide a `joint-crab` docker container that may be accessed from the public [Gammapy DockerHub repository](https://hub.docker.com/u/gammapy/dashboard/).
* [How to reproduce results with Docker](docker.md)


## Licence

The code and content in this repository is shared under the [BSD-3-Clause license](LICENSE) (read more at [OSI](https://opensource.org/licenses/BSD-3-Clause)).
