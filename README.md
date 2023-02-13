The paper is now available online (https://doi.org/10.1017/pab.2023.1), please email me if you want a copy.

This repository contains all of the R code used in the paper, plus a couple of extra related scripts. Unfortunately the data file is too large for GitHub, but is available in the official Zenodo repository (https://doi.org/10.5281/zenodo.7355273), or email me and I can send it to you.


A brief summary of the code:

Evol rate simulation.R - the code that runs the first simulation (w/o extinction selectivity, with spatial bins)

Plot simulation results.R - code for evaluating and plotting the outputs of the first simulation

Statistical tests.R - the statistical tests used to evaluate simulation performance

Selectivity simulation.R - the code that runs the second simulation (with extinction selectivity)

Wrangle occurrences.R - code for cleaning the PBDB occurrence dataset

Rotating palaeo-occurrences.R - code for pre- and post-processing the occurrences for spatial rotation with GPlates

Sampling through space-time.R - code for generating summary statistics from the occurrence data

Latitudinal rates from fossils.R - applying the metrics to the occurrence data


Extras:

Raw, BC and GF global rates.R - code to calculate raw, boundary crosser and gap-filler rates on global datasets

Simulation demo.R - a much simpler version of the first simulation, developed for a seminar on building simulations in R for BES Palaeoecology
