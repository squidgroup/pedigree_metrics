# Code Readme

**Project title:** Pedigree metrics: What are they good for?

**Code Author:** Joel Pick

**Code License:** MIT (see R/LICENSE)



## R Code Files

### 00_functions.R
- Contains the functions that are commonly used in other code files. This does not need to be run - it is called within the other code files when needed.

### 01_young_postma.R
- Extracts the studies using animal models from the data of Young and Postma 2023
- Imports from *Data/Raw/young_postma_2023_data.csv*
- Creates
    - *Data/Raw/to_extract.csv*

### 02_pedigree_stats.R
- Generates pedigree statistics, analytical approximations of sampling variance, sampling variance from simulations, minimum detectable effect size and bias from unmodelled maternal effects for each of 53 pedigrees.
- Imports pedigree files from *Data/Raw/pedigrees*
- Creates
    - *Data/Processed/ped_stats.Rdata*
    - *Data/Processed/ped_h2_se.Rdata*
    - *Data/Processed/ped_h2_se_simple.Rdata*
    - *Data/Processed/ped_sims.Rdata*
    - *Data/Processed/ped_sims_low.Rdata*
    - *Data/Processed/ped_sims_high.Rdata*
    - *Data/Processed/ped_mde.Rdata*
    - *Data/Processed/ped_mbias.Rdata*


### 03_phenotyped.R
- Sampling variances from simulations and analytical approximations for pedigrees with missing phenotypic data
- Imports pedigree files from *Data/Raw/pedigrees*
- Creates
    - *Data/Processed/ped_missing.Rdata*

### 04_grms.R
- Sampling variances from simulations and analytical approximations for 6 GRMs
- Imports pedigree files from *Data/Raw/grms*
- Creates
    - */Data/Processed/grm_stickleback_dat.Rdata*
    - */Data/Processed/grm_td_dat.Rdata*
    - */Data/Processed/grm_dragon_dat.Rdata*
    - */Data/Processed/grm_koala_dat.Rdata*
    - */Data/Processed/grm_roe_dat.Rdata*
    - */Data/Processed/grm_soays_dat.Rdata*


### Figure_bias.R
- Creates Figures comparing pedigree metrics to biases in h2 estimation 
- Imports
    - *Data/Processed/ped_mbias.Rdata*
    - *Data/Processed/ped_stats.Rdata*
- Creates
    - Figure 8 - *Figures/Fig_bias.pdf*
    - Figure S8 - *Figures/FigSM_bias_all.pdf*

### Figure_correlations.R
- Creates figure of correlations between pedigree metrics
- imports 
    - *Data/Processed/ped_stats.Rdata*
- Creates 
    - Figure 2 - *Figures/Fig_cor.pdf*


### Figure_grms.R
- Creates figures comparing estimates of precision from GRMs and pedigrees
- Imports:
    - *Data/Processed/ped_stats.Rdata*
    - *Data/Processed/ped_h2_se.Rdata*
    - *Data/Processed/ped_sims.Rdata*
    - */Data/Processed/grm_stickleback_dat.Rdata*
    - */Data/Processed/grm_td_dat.Rdata*
    - */Data/Processed/grm_dragon_dat.Rdata*
    - */Data/Processed/grm_koala_dat.Rdata*
    - */Data/Processed/grm_roe_dat.Rdata*
    - */Data/Processed/grm_soays_dat.Rdata*
- Creates
    - Figure 7 - *Figures/Fig_grm.pdf*


### Figure_h2se.R
- Compares analytical approximations of sampling variance of h2 and expected sampling variance based on simulation
- Imports
    - *Data/Processed/ped_h2_se.Rdata*
    - *Data/Processed/ped_h2_se_simple.Rdata*
    - *Data/Processed/ped_sims.Rdata*
- Creates
    - Figure 4 - *Figures/Fig_h2_se.pdf*


### Figure_MDE.R
- Creates figures comparing pedigree metrics to minimum detectable effect size
- Imports 
    - *Data/Processed/ped_stats.Rdata*
    - *Data/Processed/ped_mde.Rdata*
- Creates
    - Figure 5 - *Figures/Fig_mde.pdf*
    - Figure S6 - *Figures/FigSM_mde_all.pdf*


### Figure_metrics.R
- Creates figures on reporting of pedigree metrics
- imports
    - *Data/Raw/studies_info.csv*
    - *Data/Raw/studies_results.csv*
- Creates 
    - Figure 1 - *Figures/Fig_metrics.pdf*
    - Figure S1 - *Fig_metrics_all.pdf*


### Figure_missing.R
- Create Figure showing the impact of missing phenotypic data
- Imports:
    - *Data/Processed/ped_missing.Rdata*
- Creates
    - Figure 6 - *Figures/Fig_missing.pdf*


### Figure_precision.R
- Creates Figures comparing pedigree metrics to precision in h2 estimates 
- Imports:
    - *Data/Processed/ped_sims.Rdata*
    - *Data/Processed/ped_stats.Rdata*
- Creates:
    - Figure 3 - *Figures/Fig_precision.pdf*
    - Figure S3 - *Figures/FigSM_precision_all.pdf*



### SM_grms.R
- Comparison of Analytical approximations of sampling variance from GRMs with and without making them positive definite
- Imports:
    - */Data/Processed/grm_stickleback_dat.Rdata*
    - */Data/Processed/grm_td_dat.Rdata*
    - */Data/Processed/grm_dragon_dat.Rdata*
    - */Data/Processed/grm_koala_dat.Rdata*
    - */Data/Processed/grm_roe_dat.Rdata*
    - */Data/Processed/grm_soays_dat.Rdata*
- Creates
    - */Data/Processed/grm_dat_compare.Rdata*
    - Figure 7 - */Data/Processed/grm_dat_compare.Rdata*
    - Figure S7 - *Figures/FigSM_grm_compare.pdf*

### SM_h2.R
- Compares relationship between pedigree metrics and sampling variance of h2 estimates at different simulated h2s
- Imports:
    - *Data/Processed/ped_stats.Rdata*
    - *Data/Processed/ped_sims_low.Rdata*
    - *Data/Processed/ped_sims_high.Rdata*
    - *Data/Processed/ped_sims.Rdata*
- Creates
    - Figure 4 - *Figures/FigSM_h2.pdf*

### SM_taxa.R
- Distribution of pedigree sizes across taxa
- Imports:
    - pedigree files from *Data/Raw/pedigrees*
    - *Data/Processed/ped_stats.Rdata*
    - *Data/Raw/studies_info.csv*
    - *Data/Raw/dryad_ped_info.csv*
- Creates
    - Figure S2 - *Figures/FigSM_taxa.pdf*


### SM_VG.R
- Compares analytical approximations of sampling variance of h2 estimates and expected sampling variance based on simulation cross different h2s
- Imports
    - *Data/Processed/ped_h2_se.Rdata*
    - *Data/Processed/ped_h2_se_simple.Rdata*
    - *Data/Processed/ped_sims.Rdata*
    - *Data/Processed/ped_sims_low.Rdata*
    - *Data/Processed/ped_sims_high.Rdata*
- Creates
    - Figure S5 - *Figures/FigSM_VG.pdf*




## Software
R 4.4.0

## R packages
MCMCglmm_2.36
nadiv_2.18.0
parallel_4.4.0
squidSim_0.2.6  
gremlin_1.1.0  
Matrix_1.7-0    
pedtricks_0.5.0
scales