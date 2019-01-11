## Analysis

The analysis results from this paper, including the plots, are reproducible.

If you have installed the `joint-crab` bundle as it is described in the
[README](README.md), you can type in the terminal the following commands to reproduce the analysis and figures. 
You may choose to perform the whole analysis with a single command or do it in a step-by-step basis.

#### All steps in a sigle command

    $ ./make.py all

#### List specific tasks or subcommands

    $ ./make.py --help

<details>
<summary>Subcommands</summary>
<pre>
  all              Run all steps.
  clean            Clean out results folder.
  extract-spectra  Extract 1d spectra
  fit-errorbands   Compute flux error bands.
  fit-spectra      Execute spectrum fit.
  fit-systematics  Fit that includes systematics.
  maps             Make and plot sky maps.
  plot-contours    Plot contours.
  plot-counts      Plot counts spectra.
  plot-errorbands  Plot SED error bands.
  plot-seds        Plot SEDs.
  provenance       Write `results/provenance.yaml`.
  summary-data     Write summary for data.
  summary-results  Write summary for results.
</pre>
</details>

#### Step-by-step

**1- Let's clean eventual previous results produced.**

    $ ./make.py clean
    
**2- We perform the analysis in a procedural way.**

    $ ./make.py maps
    $ ./make.py extract-spectra
    $ ./make.py fit-spectra
    $ ./make.py fit-systematics
    $ ./make.py fit-errorbands

**3- Then we gather the results in tables.**

    $ ./make.py summary-data
    $ ./make.py summary-results

* Table 1: [results/summary/data.md](results/summary/data.md) | [results/summary/data.tex](results/summary/data.tex)
* Table 2: [results/summary/results.md](results/summary/results.md) | [results/summary/results.tex](results/summary/results.tex)


**4- Finally, let's produce the plots.**

`$ ./make.py plot-counts`
    
* Fig. 1: [results/figures/counts_spectra.png](results/figures/counts_spectra.png) | [results/figures/counts_spectra.pdf](results/figures/counts_spectra.pdf)
    
`$ ./make.py plot-seds`

* Fig. 2: [results/figures/crab_sed_fit.png](results/figures/crab_sed_fit.png) | [results/figures/crab_sed_fit.pdf](results/figures/crab_sed_fit.pdf)
  
`$ ./make.py plot-errorbands`

* Fig. 4: [results/figures/errorband_sed_veritas.png](results/figures/errorband_sed_veritas.png) | [results/figures/errorband_sed_veritas.pdf](results/figures/errorband_sed_veritas.pdf)
  
`$ ./make.py plot-contours`
    
* Fig. 3: [results/figures/contours.png](results/figures/contours.png) | [results/figures/contours.pdf](results/figures/contours.pdf)
* Fig. 5: [results/figures/contours_systematics.png](results/figures/contours_systematics.png) | [results/figures/contours_systematics.pdf](results/figures/contours_systematics.pdf)
  

