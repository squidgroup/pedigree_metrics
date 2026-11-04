# Data README

**Project title:** Pedigree metrics: What are they good for?

**Data License:** CC-BY 4.0


## Raw/

Directory containing raw data

### pedigrees/
Directory containing pedigrees. Information about the pedigrees (including their origin) is in *Data/Raw/dryad_ped_info.csv* and *Data/Raw/studies_info.csv*, where the filenames link to the IDs in these files. *Data/Raw/dryad_ped_info.csv* contains information for the files beginning with a 'D' and  *Data/Raw/studies_info.csv* contains information for the numeric filenames. Each file contains the variables:

id: individual ID

dam: mother ID

sire: father ID


### grms/
Directory containing GRMs.

#### dragon_GRM_2020.txt
#### koala_GRM_2021.txt
#### stickleback_GRM.csv
#### taz_devil_grm.csv
#### roe_deer_grm.Rdata



### dryad_ped_info.csv
Information about studies that used animal model to estiamte h2 from a search of Dryad repositories.

**Variables:**

id: pedigree ID

Link: url for the dryad repository

Paper: Article title

species_english: English species names

taxonomic class: Taxonomic group (bird/mammal/reptile/fish)

year: Year the article was published

Pedigree: Was a pedigree archived (Yes/No)

A: Was a A matrix archived (Yes/No)

InvA: Was a inverse A matrix archived (Yes/No)



### studies_info.csv
Information about subset of articles using animal models from Young and Postma 2023.

**Variables:**

id: article ID

reference: article reference

species_english: English species names

species_latin: Latin species names

population: Location of the population

country: Country the population is in

yp_study_ID: study ID in young and postma 2023

year: Year the article was published

taxonomic class: Taxonomic group (bird/mammal)


### studies_results.csv
Extracted data from subset of articles using animal models in Young and Postma 2023.


**Variables:**

Timestamp: Time and data of data extraction 

Study ID: article ID

paper link: article url 

Metrics reported?: Were any metrics reported

Metrics reported: List of metrics reported from a set list 

Other metrics: Any different metrics that were reported

Is pedigree error mentioned?: Did the article mention pedigree error (Yes/No)

Is pruning mentioned?: Did the article mention pedigree pruning (Yes/No)

Are pre-analysis simulations mentioned?: Did the article mention pre-analysis simulations (Yes/No)

Pedigree provided?: Has a pedigree been archived(Yes/No)

Pedigree link: url link to the pedigree

A or invA matrix provided?: Has a A or inv A matrix been archived(Yes/No)

Link to A/invA: url link to the pedigree

Notes: Any notes


### to_extract.csv
Subset of articles using animal models from young and Postma 2023

**Variables:**

study.ID: article ID

year: Year the article was published

reference: article reference

species.english: English species names

species.latin: Latin species names

taxonomic.class: Taxonomic group (bird/mammal)

population: Location of the population

country: Country the population is in

estimate.ID: ID of h2 estimate



### young_postma_2023_data.csv
Data from Young and Postma 2023
Original data repository: [https://doi.org/10.34894/RIVFHW](https://doi.org/10.34894/RIVFHW)

**Variables:**

estimate.ID: ID of the effect size

study.ID: ID of the study

method: method of estimating heritability

h2: Heritability estimate

SE.h2: Standard error of heritability estimate




## Processed/

Contains processed data files

### grm_dat_compare.Rdata
Data from comparison of analytical approximations of sampling variance h2 estimates from GRMs, with and without making the GRM positive definate

### grm_dragon_dat.Rdata
Sampling variance of h2 from water dragon GRM calculated via simulation and via analytical approximation

### grm_koala_dat.Rdata
Sampling variance of h2 from koala GRM calculated via simulation and via analytical approximation

### grm_roe_dat.Rdata
Sampling variance of h2 from roe deer GRM calculated via simulation and via analytical approximation

### grm_soays_dat.Rdata
Sampling variance of h2 from soay sheep GRM calculated via simulation and via analytical approximation

### grm_stickleback_dat.Rdata
Sampling variance of h2 from stickleback GRM calculated via simulation and via analytical approximation

### grm_td_dat.Rdata
Sampling variance of h2 from Tasmanian devil GRM calculated via simulation and via analytical approximation

### ped_h2_se_simple.Rdata
Analytical approximation of sampling variance of h2 estimates using method of Visscher et al. 2014, from 53 pedigrees

### ped_h2_se.Rdata
Analytical approximation of sampling variance of h2 estimates using method of Visscher & Goddard 2015, from 53 pedigrees

### ped_mbias.Rdata
Bias in h2 estimates from simulations, were maternal effects were simulated but not modelled, from 53 pedigrees

### ped_mde.Rdata
Minimum detectable effect size in h2 estimates from simulations, from 53 pedigrees

### ped_missing.Rdata
Precision in h2 estimates from simulations with and without missing phenotypic data, from 53 pedigrees

### ped_sims_high.Rdata
Precision in h2 estimates from simulations with a h2 of 0.5, from 53 pedigrees

### ped_sims_low.Rdata
Precision in h2 estimates from simulations with a h2 of 0.1, from 53 pedigrees

### ped_sims.Rdata
Precision in h2 estimates from simulations with a h2 of 0.3, from 53 pedigrees

### ped_stats.Rdata
Pedigree metrics from 53 pedigrees