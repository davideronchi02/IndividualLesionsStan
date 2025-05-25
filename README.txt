Individual Lesion Model

This repository contains R scripts and required data for fitting tumor growth and survival models with individual lesion-level resolution. The models are implemented in Stan using the cmdstanr interface and apply to both continuous-time and finite-time modeling approaches.

Repository Structure
Copia
Modifica
.
├── RFiles/
│   ├── IndividualLesionModel_ContinuousTime.R
│   └── IndividualLesionModel_FiniteTime.R
├── Models/
│   ├── ModelMain_ContinuousTime.stan
│   └── ModelMain_FiniteTime.stan
├── datafiles/
│   ├── DataFileLesion.csv
│   ├── dataAUC.csv
│   └── datafileSurv.csv
├── TreatmentStart.csv
└── README.md

Description
RFiles/: Contains two main R scripts to run the model fitting.

IndividualLesionModel_ContinuousTime.R: Script for the continuous-time model variant.

IndividualLesionModel_FiniteTime.R: Script for the finite-time model variant.

Models/: Stan models used for MCMC sampling with cmdstanr.

datafiles/: Input data used for model fitting. These include:

Individual lesion dynamics,

AUC-derived treatment data,

Survival information.

TreatmentStart.csv: Treatment administration metadata for survival modeling.

Requirements
R with the following packages:

cmdstanr, posterior, bayesplot, ggplot2, plyr, gridExtra, dplyr, tidyverse, reshape2, huxtable, tidyr

A working installation of CmdStan (configured via cmdstanr::set_cmdstan_path()).

Running the Models
Make sure you have CmdStan installed and accessible via cmdstanr.

Prepare the folder structure as above.

Run either of the scripts in RFiles/ to fit the models:

R
Copia
Modifica
source("RFiles/IndividualLesionModel_ContinuousTime.R")
# or
source("RFiles/IndividualLesionModel_FiniteTime.R")
Each script:

Loads and preprocesses data,

Prepares data lists for Stan,

Specifies initial values,

Compiles the Stan model,

Runs MCMC sampling.

Contact
For questions, data requests, or troubleshooting, please contact:

Davide Ronchi
davide.ronchi@unipv.it