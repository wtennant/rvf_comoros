# Supporting data and code for *'Drivers of Rift Valley fever virus persistence and the impact of control measures in a spatially heterogeneous landscape: the case of the Comoros archipelago, 2004--2015'*

This repository contains the data and code which were used to produce the
manuscript *'Drivers of Rift Valley fever virus persistence and the impact of control measures in a spatially heterogeneous landscape: the case of the Comoros archipelago, 2004--2015'*. A preprint for this manuscript
will be available shortly.

In our manuscript, we fitted an epidemiological model to Rift Valley fever surveillance data conducted
in livestock across the Comoros archipelago between 2004 and 2015. Using a fitted
model, we then tested different transmission reduction scenarios (e.g inter-island movement restrictions)
on the epidemiology of Rift Valley fever in the archipelago.

Below we outline the data used to fit the epidemiological model and provide
a description for the code used to produce out manuscript.

## Data
### Serological data
All serological data used was secondary in our study. The Comoros archipelgo consists of the Union of Comoros and Mayotte.

Data for the Union of Comoros was is presented in:

- Roger, M., Girard, S., Faharoudine, A., Halifa, M., Bouloy, M., Cetre-Sossah, C., and Cardinale,
E. Rift Valley fever in ruminants, Republic of Comoros, 2009. Emerging infectious diseases, 17(7):
1319, 2011.

- Roger, M., Beral, M., Licciardi, S., Soule, M., Faharoudine, A., Foray, C., Olive, M.-M., Maquart,
M., Soulaimane, A., Kassim, A. M., et al. Evidence for circulation of the Rift Valley fever virus
among livestock in the union of Comoros. PLoS Negl Trop Dis, 8(7):e3045, 2014.

All data for Mayotte has been presented in:

- Metras, R., Cavalerie, L., Dommergues, L., Merot, P., Edmunds, W. J., Keeling, M. J., Cetre-Sossah,
C., and Cardinale, E. The epidemiology of Rift Valley fever in Mayotte: insights and perspectives
from 11 years of data. PLoS neglected tropical diseases, 10(6):e0004783, 2016.

- Metras, R., Fournie, G., Dommergues, L., Camacho, A., Cavalerie, L., Merot, P., Keeling, M. J.,
Cetre-Sossah, C., Cardinale, E., and Edmunds, W. J. Drivers for Rift Valley fever emergence in
Mayotte: a Bayesian modelling approach. PLoS neglected tropical diseases, 11(7):e0005767, 2017.


Summaries of this data are provided in `data/seroprev.csv`. These summaries were used to fit our epidemiological
model.  This data gives the number of animals which were tested and the number which were
positive for RVF-specific IgG antibodies per island, age group and time point. The columns of this data set are as follows:

- `ISLAND_ID`: Each island in the archipelago was given a unique identifier: Anjouan (0), Grande Comore (1), Mayotte (2), Mohéli (3).
- `N_POSITIVE`: Number of tested animals with RVF-specific IgG antibodies.
- `N_TESTED`: Number of tested animals.
- `AGE_GROUP`: The age group of tested animals. Values 1–9 correspond to animals ages 0–8. A value of 9 corresponds to animals greater than or equal to 9 years old. Animals only identified as either infant or adult are marked with a -3 and -2 respectively. Tests where animal age was not identified are under age group -1.
- `EPI_YEAR_START` and `EPI_YEAR_END`: The calendar year which the epidemiological year started and ended.
- `EPI_MONTH`: The month (from July) the data refers to.Epidemiological years began in July of each calendar. A value of 13 corresponds to the entire epidemiological year.

### Normalised Difference Vegetation Index (NDVI) data
Normalised Difference Vegetation Index (NDVI) was used to incur time-dependent transmission rates of RVF between livestock in our transmission model. The data provided in `data/ndvi.csv`. This file gives summarised NDVI data for each island in the Comoros archipelago and simulation week (note: there are 4 simulation weeks per month for ease of computation). The columns of this data are as follows:

- `ISLAND_ID`: Each island in the archipelago was given a unique identifier: Anjouan (0), Grande Comore (1), Mayotte (2), Mohéli (3).
- `DATE`: Not used in computation. Gives approximate corresponding date for NDVI data.
- `WEEK`: Simulation week of a month.
- `MONTH`: Month of the calendar year.
- `YEAR`: Calendar year
- `NDVI`: Value of NDVI at the given location and time point.

## Code
The custom code used in our manuscript was split between `C++` and `R`. The former was used for model fitting and simulation of control scenarios, whereas the latter is used for data visualisation. Please note that for these programs to work correctly, directories need to be appropriately set in the configuration files and script files.

### Model fitting
The summarised serological data was fitted in a fully Bayesian framework by executing custom `C++` code. This code consists of an epidemiological model as described in our manuscript, an observation model describing how empirical serological samples correspond to epidemiological model outputs and a fitting algorithm (namely MCMC) used to estimate the posterior distribution of model parameters.

This code can be found under the folder `cpp/mcmc/`. All files are fully commented and describe what they do. In summary, the *important* source files here are:

- `config.h`: Used to specify configurations of model fitting — including length of MCMC chain, number of independent chains, length of exploration period, adaption period, mixing period (see `mcmc.cpp`) — and epidemiological model specificications.
- `likelihood.cpp`: Functions to calculate the likelihood of empirical data given the simulated population-level seroprevalence from the epidemiological model.
- `mcmc.cpp`: Functions for running a multi-chain adaptive Metropolis-Hastings MCMC algorithm. This algorithm begins with a fixed number of independent chains. At each step in a chain, a new set of parameters (to be estimated) is proposed. The epidemiological model run with these parameters and likelihood calculated. The proposed parameter set is either rejected or accepted depending on the likelihood and prior of the proposed and previous parameter sets. Initially, proposals are taken as a random walk as to explore the parameter space (exploratory phase). After a set number of iterations, proposals of parameters are taken from an adaptive algorithm which considers correlations between parameters of interest. Finally, proposals are taken across all chains to ensure no chain became stuck in one (relatively lower likelihood) region of the parameter space. Full details of the algorithm are given in the manuscript.
- `parameters.cpp`: Parameters used in the epidemiological model.
- `simulation.cpp`: Function simulate the epidemiological model forward in time. Details of the epidemiological model are given in the manuscript.
- `theta.cpp`: Space containing the chains of parameters to be fitted to the empirical data.

### Preparation of posteriors for scenario testing
One `R` script was used to prepare posterior distributions for scenario testing. This script file can be found in `R/prepare_post.R`. After checking that all chains had converged, this file discards a set burn-in period of each chain generated by the model fitting algorithm and randomly (jointly) samples parameter values from the posterior distribution. This sampled posterior distribution is then used in scenario testing.

### Control scenario testing
Once a model had been fitted, we sampled from the posterior distributions of the parameters of interest and executed our epidemiological under different control scenarios. The scenarios that we tested were restrictions on movement between islands and fixed time-independent reductions in the transmision rate of the disease. This code can be found under the folders `cpp/movement_restrict` and `cpp/trasmission_reduce` respectively. There are two new *important* files in this folder:

- `posterior.cpp`: contains information on the posterior samples taken from a previously executed model fit **which has then been prepared in `R`**.
- `scenario.cpp`: parameters used during scenario testing.

### Data visualisation
`R` was used to visualise all empirical and simulation data. Each script file in the `R/` folder fully describes the function of each script. Some files have two versions, one with a suffix `_alt` and one without. Script files with the suffix `_alt` use simulated results generated from the control scenario testing `C++` programs. Script files without the suffix (yet there is one available) use simulated results from the model fitting program.
