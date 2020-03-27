# StigmergyDisease
[![DOI](https://zenodo.org/badge/185241410.svg)](https://zenodo.org/badge/latestdoi/185241410)

Code for: "A mechanistic, stigmergy model of territory formation in solitary animals: Territorial behavior can dampen disease prevalence but increase persistence"
* Simulation models to explore the consequences of stigmergy on territory formation and disease spread

## Results
* `stig_summary.csv`- compilation of simulation results with parameter sets and outcomes recorded
* `partyRF_logit500.csv`- random forest results for outbreak success
* `partyRF_logitdur500.csv`- random forest results for outbreak duration
* `partyRF_logitprev500.csv`- random forest results for maximum prevalence

## Simulation Code
* `StigmergyFunctions.R`- file containing all functions needed to run stigmergy based simulations
* `ToyModelStigmergy.R`- simulates a single simulation for a given parameter set; used in `PlottingMovement.R` to plot movement trajectories
* `StigmergyLoop.R`- Building up `ToyModelStigmergy.R` to run in a loop
* `StigmergyLoopFunction.R`- Setting up `StigmergyLoop.R` to run as a function `StigLoop` as a precursor to running scripts in parallel
* `StigmergyLoopRepFunction.R`- Allowing each parameter set to be divided into reps for particularly long running simulations
* `Rslurm_parallel.R`- run code in parallel on cluster using `rslurm` package

## Analaysis
* `CSVmerge.R`- read raw simulation output files into a single csv + exploratory plotting of results
* `RFParty.R`- random forest analysis using `party` package
* `RFStig.R`- random forest analysis using `randomforest` package
* `PlottingRFStigmergy.R`- plotting random forest results
* `StigmergyFigures.Rmd`- code to reproduce publication figures
