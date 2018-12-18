# Towards multi-instrument and reproducible gamma-ray analysis

[![DOI](https://zenodo.org/badge/146204837.svg)](https://zenodo.org/badge/latestdoi/146204837)

[![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/open-gamma-ray-astro/joint-crab/master?urlpath=lab/tree/joint-crab)

This repository contains content related to the paper *Towards multi-instrument and reproducible gamma-ray analysis*.

This `joint-crab` bundle allows to reproduce the published results by installation and execution in the local desktop, using a Docker image or via the MyBinder cloud service.

<details>
<summary>Content</summary>
<pre>
.
├── binder
│   └── environment.yml
├── data
│   ├── README.md
│   ├── fact
│   ├── fermi
│   ├── hess
│   ├── magic
│   ├── other
│   └── veritas
├── joint_crab
│   ├── __init__.py
│   ├── conf.py
│   ├── errorbands.py
│   ├── extract_lat.py
│   ├── extract_ogip_spectra.py
│   ├── fit_spec.py
│   ├── maps.py
│   ├── models.py
│   ├── plot_contours.py
│   ├── plot_counts.py
│   ├── plot_errorbands.py
│   ├── plot_seds.py
│   ├── provenance.py
│   ├── summary_data.py
│   ├── summary_results.py
│   ├── systematics.py
│   └── utils.py
├── presentations
│   ├── README.md
│   └── 2018tevpa
│       ├── README.md
│       ├── poster.pdf
│       └── slides.pdf
├── results
│  ├── errorbands
│  ├── figures
│  ├── fits
│  ├── maps
│  ├── spectra
│  └── summary
├── 1_data.ipynb
├── 2_counts.ipynb
├── 3_systematics.ipynb
├── 4_naima.ipynb
├── 5_crab_pulsar_nebula_sed.ipynb
├── Dockerfile
├── LICENSE
├── README.md
├── analysis.md
├── docker.md
└── make.py
</pre>
</details>


##  The `joint-crab` bundle

You can execute the provided scripts in a specific `joint-crab` *conda environment*. This will allow you to reproduce the results published in the paper.

####  Set-up the environment

As a first approach to assure the reproducibility of the results we provide a [conda configuration file](./environment.yml) to build a virtual environment with pinned versions for the software dependencies.

Hence, as a requisite you need [Anaconda](https://www.anaconda.com/download/) or
[miniconda](https://conda.io/miniconda.html) - Python>=3.6 software installed in your desktop. Once you have installed this software, and downloaded the content of this `joint-crab` repository, you can type in the terminal the following commands at the top-level of the `joint-crab` folder:

    $ conda env create -f binder/environment.yml
    $ source activate joint-crab

####  Reproduce the results
* [How to reproduce results and figures with `make.py`](analysis.md)

####  Check out the  notebooks

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

The code anc content in this repository is shared under the [BSD-3-Clause license](LICENSE) (read more at [OSI](https://opensource.org/licenses/BSD-3-Clause)).
