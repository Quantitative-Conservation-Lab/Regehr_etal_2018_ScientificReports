# [Integrated Population Modeling Provides the First Empirical Estimates of Vital Rates and Abundance for Polar Bears in the Chukchi Sea](https://doi.org/10.1038/s41598-018-34824-7)

### Regehr, E. V., N. J. Hostetter, R. R. Wilson, K. D. Rode, M. S. Martin, and S. J. Converse

### Scientific Reports

### Code/Data DOI: [https://doi.org/10.5061/dryad.692jb15](https://doi.org/10.5061/dryad.692jb15)

### Please contact the first author for questions about the code or data: Eric Regehr (eregehr@uw.edu)
__________________________________________________________________________________________________________________________________________

## Abstract:
Large carnivores are imperiled globally, and characteristics making them vulnerable to extinction (e.g., low densities and expansive ranges) also make it difficult to estimate demographic parameters needed for management. Here we develop an integrated population model to analyze capture-recapture, radiotelemetry, and count data for the Chukchi Sea subpopulation of polar bears (Ursus maritimus), 2008–2016. Our model addressed several challenges in capture-recapture studies for polar bears by including a multievent structure reflecting location and life history states, while accommodating state uncertainty. Female breeding probability was 0.83 (95% credible interval [CRI] = 0.71–0.90), with litter sizes of 2.18 (95% CRI = 1.71–2.82) for age-zero and 1.61 (95% CRI = 1.46–1.80) for age-one cubs. Total adult survival was 0.90 (95% CRI = 0.86–0.92) for females and 0.89 (95% CRI = 0.83–0.93) for males. Spring on-ice densities west of Alaska were 0.0030 bears/km2 (95% CRI = 0.0016–0.0060), similar to 1980s-era density estimates although methodological differences complicate comparison. Abundance of the Chukchi Sea subpopulation, derived by extrapolating density from the study area using a spatially-explicit habitat metric, was 2,937 bears (95% CRI = 1,552–5,944). Our findings are consistent with other lines of evidence suggesting the Chukchi Sea subpopulation has been productive in recent years, although it is uncertain how long this will continue given sea-ice loss due to climate change.

## Code 
1. [IPM_analysis](./IPM_analysis/): This folder contains the code to load data and run the Integrated Population Model and goodness-of-fit for Chukchi Sea polar bears. It also contains the code to generate the JAGS file for the Bayesian implementation.

2. [recruitment_analysis](./recruitment_analysis/): This folder contains the code to load data and run the recruitment analysis (yearlings [c1s] per adult female). It also contains the code to generate the JAGS file for the Bayesian implementation.

3. [composition_analysis](./density_extrapolation_analysis/): This folder contains the code to run the density extrapolation.


## Data
Datasets used in this project are all found in the [data folder](./data):

1) Regehr2018_jags_data.txt       
Description: Formatted data to run the Integrated Population Model. See manuscript for detailed description of data.

2) Regehr2018_state_inits.txt         
Description: Initial values to run the Integrated Population Model. 

3) Regehr2018_C1perAF_data.txt      
Description: Formatted data to run the recruitment analysis.


3) Regehr2018_dextrp_data.txt      
Description: Formatted data to run the density extrapolation analysis.
