# Data and R scripts for the paper 'TITLE'

This repository contains data and scripts used in the following paper:

TITLE

The repository is organised into four main directories: data, runscripts, results, and manuscript. 
The _data_ directory contains the rice radiocarbon database, pre-processing scripts, and simulated datasets; _runscript_  contains all scripts for running the analyses in the paper, _src_ contains custom R functions, _results_ contains R image files with the analyses output, and _manuscript_ contains scripts for reproducing figures and tables included in the manuscript. 

# Dataset

Dataset required for all analyses are contained in the directory _data_. The main dataset is a CSV file `R14CDB.csv` which contains site location data and radiocarbon ages of charred rice (_O. sativa_) from the Japanese islands collated from excavation reports, online databases, and journal articles. `prefecture_region_match.csv` is a look-up table defining the regions for the hierarchichal phase model analyses (see below). The file `prepare_data.R` contains an R script that pre-processes information contained in these two files and generates the R image file `c14rice.RData`, which contains R objects required for the Bayesian analyses. Finally the R image files `tactical_sim_gpqr.RData` and `tactical_sim_phase.RData` contains simulated dataset that were used to assess the robustness of the proposed methods. These were generated using the R scripts `tactical_gpqr_sim.R` and `tactical_phase_sim.R` contained in the directory `runscripts`.    


# Bayesian Analyses

## Quantile Regression and Gaussian Process Quantile Regression

## Hierarchichal Phase Model

# File Structure

### data
* `R14CDB.csv` ... Spreadsheet containing radiocarbon data on charred rice.
* `prefecture_region_match.csv` ... lookup table defining the regions for the hierarchical phase model analyses. 
* `key_sites_for_map.csv` ... contains coordinates of key site locations displayed on figure 1.
* `prepare_data.R` ... R scripts for pre-processing `R14CDB.csv` and `prefecture_region_match.csv`. Generates the R image file `c14rice.RData`.
* `c14rice.RData` ... R image file containing R objects required for analyses. Generated using `prepare_data.R`, `R14CDB.csv`, and `prefecture_region_match.csv`.
* `tactical_sim_gpqr.RData` ... R image file containing simulated dataset for testing purposes. Generated using the script `tactical_gpqr_sim.R` in the `runscripts` directory.
* `tactical_sim_phase.RData` ... R image file containing simulated dataset for testing purposes. Generated using the script `tactical_phase_sim.R` in the `runscripts` directory.

### src
* `gpqrRSim.R` ... Contains an R function (`gpqrSim()`) for simulating dispersal dates using Gaussian Process model.
* `diffplot.R` ... Contains an R function (`diffDens()`) for visualising the posterior density of a difference between two parameters.
* `orderPPlot.R` ... Contains an R function (`orderPPlot`) for cretating a matrix visualising the before/after probabilities of parameters defining particular chronological events.  

### runscripts
* `quantreg.R` ... Contains the R script for running Bayesian and non-Bayesian quantile regression. Results are stored in the R image file `quantreg_res.RData` in the _results_ directory. 
* `phasemodel_a.R` and `phasemodel_b.R` ... Contains the R script for running the unconstrained (model a) and constrained (model b) versions of the hierarhichal Bayesian phase models. Results are stored in the R image files `phase_model_a.RData` and `phase_model_b.RData` in the _results_ directory.
* `phasemodel_tactical.R` ... Contains the R script for running a hierarhichal Bayesian phase model on the simulated dataset `tactical_sim_phase.RData`. Results are stored in the R image file `phasemodel_tactsim.RData` in the _results_ directory. 
* `tactical_phase_sim.R` ... Contains the R script for generating a simulated dataset for testing hierarchichal phase models. The output is stored in the R image file `tactical_sim_phase.RData` in the _data_ directory. 
* `gpqr_tau90.R` ... Contains the R script for running a GPQR model with tau set at 0.90.  Results are stored in the R image file `gpqr_tau90.RData` in the _results_ directory. 
* `gpqr_tau99.R` ... Contains the R script for running a GPQR model with tau set at 0.99.  Results are stored in the R image file `gpqr_tau99.RData` in the _results_ directory. 
* `gpqr_tactical.R` ... Contains the R script for running a GPQR model on the simulated dataset `tactical_sim_gpqr.RData`. Results are stored in the R image file `gpqr_tactsim.RData` in the _results_ directory.   
* `tactical_gpqr_sim.R` ... Contains the R script for generating a simulated dataset for testing GPQR models. The output is stored in the R image file `tactical_sim_gpqr.RData` in the _data_ directory.


### results
* `quantreg_res.RData` ... Output of quantile regression analyses. Generated using the script contained in `quantreg.R` in the _runscripts_ directory. 
* `phase_model_a.RData` ... Output of hierarhichal Bayesian phase models with no constraints (model a). Generated using the script contained in `phasemodel_a.R`  in the _runscripts_ directory. 
* `phase_model_b.RData` ... Output of hierarhichal Bayesian phase models with constraints (model b). Generated using the script contained in `phasemodel_b.R`  in the _runscripts_ directory. 
* `phasemodel_tactsim.RData` ... Output of the hierarhichal Bayesian phase model on simulated dataset. Generated using the script contained in `phasemodel_tactical.R` in the _runscripts_ directory. 
* `gpqr_tau90.RData` ... Output of the GPQR model with tau=0.9. Generated using the script contained in `gpqr_tau90.R`  in the _runscripts_ directory. 
* `gpqr_tau99.RData` ... Output of the GPQR model with tau=0.99. Generated using the script contained in `gpqr_tau99.R`  in the _runscripts_ directory. 
* `gpqr_tactsim.RData` ... Output of the GPQR model on simulated dataset. Generated using the script contained in `gpqr_tactical.R` in the _runscripts_ directory. 

### manusript
* `main_figures_log.R` ... R scripts for generating figures for the main text. 
* `/main_figures/*.pdf` ... Sub-directory containing pdf files of main text figures. 
* `supplementary_figures_log.R` ... R scripts for generating supplementary figures.
* `/supplementary_figures/*.pdf` ... Sub-directory containing pdf files of supplementary figures. 
* `supplementary_tables_log.R` ... R scripts for generating supplementary tables.
* `/supplementary_tables/*.csv` ... Sub-directory containing CSV files of supplementary tables. 

# R Session Info
```
R version 4.1.2 (2021-11-01)
  
attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods  
[7] base     

other attached packages:
 [1] dplyr_1.0.7         coda_0.19-4         quantreg_5.86      
 [4] SparseM_1.81        diagram_1.6.5       shape_1.4.6        
 [7] gridExtra_2.3       latex2exp_0.5.0     viridis_0.6.1      
[10] viridisLite_0.4.0   rgeos_0.5-5         sf_1.0-5           
[13] maptools_1.1-1      sp_1.4-6            rcarbon_1.4.3      
[16] nimbleCarbon_0.2.1  nimble_0.12.1       rnaturalearth_0.1.0
[19] ggridges_0.5.3      ggplot2_3.3.5       corrplot_0.92      
[22] cascsim_0.4         truncnorm_1.0-8     here_1.0.1         
[25] nvimcom_0.9-115    

loaded via a namespace (and not attached):
  [1] colorspace_2.0-2      deldir_1.0-6         
  [3] ellipsis_0.3.2        class_7.3-19         
  [5] rprojroot_2.0.2       R2HTML_2.3.2         
  [7] proxy_0.4-26          spatstat.data_2.1-0  
  [9] listenv_0.8.0         MatrixModels_0.5-0   
 [11] gsl_2.1-7.1           lubridate_1.8.0      
 [13] prodlim_2019.11.13    fansi_0.5.0          
 [15] mvtnorm_1.1-2         codetools_0.2-18     
 [17] splines_4.1.2         knitr_1.36           
 [19] polyclip_1.10-0       pROC_1.18.0          
 [21] caret_6.0-90          spatstat.linnet_2.3-0
 [23] stabledist_0.7-1      copula_1.0-1         
 [25] spatstat.sparse_2.0-0 compiler_4.1.2       
 [27] assertthat_0.2.1      Matrix_1.3-4         
 [29] tools_4.1.2           igraph_1.2.7         
 [31] gtable_0.3.0          glue_1.5.0           
 [33] reshape2_1.4.4        Rcpp_1.0.7           
 [35] spatstat_2.2-0        vctrs_0.3.8          
 [37] nlme_3.1-153          conquer_1.2.1        
 [39] iterators_1.0.13      timeDate_3043.102    
 [41] xfun_0.28             gower_0.2.2          
 [43] stringr_1.4.0         globals_0.14.0       
 [45] lifecycle_1.0.1       goftest_1.2-3        
 [47] future_1.23.0         MASS_7.3-54          
 [49] scales_1.1.1          ipred_0.9-12         
 [51] spatstat.core_2.3-2   doSNOW_1.0.19        
 [53] spatstat.utils_2.2-0  parallel_4.1.2       
 [55] rpart_4.1-15          stringi_1.7.6        
 [57] pcaPP_1.9-74          foreach_1.5.1        
 [59] e1071_1.7-9           lava_1.6.10          
 [61] matrixStats_0.61.0    rlang_0.4.12         
 [63] pkgconfig_2.0.3       moments_0.14         
 [65] lattice_0.20-45       purrr_0.3.4          
 [67] tensor_1.5            recipes_0.1.17       
 [69] tidyselect_1.1.1      parallelly_1.28.1    
 [71] plyr_1.8.6            magrittr_2.0.1       
 [73] R6_2.5.1              snow_0.4-4           
 [75] generics_0.1.1        ADGofTest_0.3        
 [77] DBI_1.1.1             pillar_1.6.4         
 [79] foreign_0.8-81        withr_2.4.2          
 [81] mgcv_1.8-36           fitdistrplus_1.1-6   
 [83] units_0.7-2           survival_3.2-13      
 [85] scatterplot3d_0.3-41  abind_1.4-5          
 [87] nnet_7.3-16           tibble_3.1.5         
 [89] future.apply_1.8.1    pspline_1.0-18       
 [91] crayon_1.4.2          KernSmooth_2.23-20   
 [93] utf8_1.2.2            spatstat.geom_2.3-0  
 [95] grid_4.1.2            data.table_1.14.2    
 [97] ModelMetrics_1.2.2.2  digest_0.6.28        
 [99] classInt_0.4-3        numDeriv_2016.8-1.1  
[101] stats4_4.1.2          munsell_0.5.0        
```

# Funding
This research was funded by the ERC grant _Demography, Cultural Change, and the Diffusion of Rice and Millets during the Jomon-Yayoi transition in prehistoric Japan (ENCOUNTER)_ (Project N. 801953, PI: Enrico Crema).

# Licence
CC-BY 3.0


