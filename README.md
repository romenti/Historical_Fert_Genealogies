## Bayesian Indirect Estimation of Historical Fertility in Europe and US Using Online Genealogical Data

[![OSF](https://img.shields.io/badge/OSF-project-blue)](https://osf.io/bpnyz/overview)
[![Generic badge](https://img.shields.io/badge/R-4.3.1-orange.svg)](https://cran.r-project.org/bin/macosx/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository contains code and materials to replicate the paper "Bayesian Indirect Estimation of Historical Fertility in Europe and US Using Online Genealogical Data"

**Abstract**: A growing number of social scientists use online genealogical data as an alternative digital census of historical populations to study past demographic dynamics. However, the non-representativeness of this data source requires the development of bias-adjusting methods to obtain accurate demographic estimates. We address this challenge by proposing an indirect estimation framework to investigate fertility trends in seven European countries and the United States of America for the historical period 1751-1910, integrating data from the big genealogical database FamiLinx with more conventional data sources. The proposed methods allow for the indirect estimation of the total fertility rate using the number of women aged 15-49 and children under age 5, while accounting for child mortality, age-specific fertility patterns, and biases. Our methodological approaches demonstrate that, when combined with reliable demographic data, online genealogical data can be fruitfully used to examine fertility patterns in countries and periods lacking well-functioning national civil registration systems.


### Authors

- [Riccardo Omenti](https://romenti.github.io/)
- [Monica Alexander](https://monicaalexander.com/)
- [Nicola Barban](https://nicolabarban.com/)

### Structure of the folder
This repository contains the code required to replicate all figures and tables presented in the main paper. To reproduce the analyses, follow the steps below:

- Download all data sources into the `data` folder. Most .RData files are available directly in the GitHub repository, while larger datasets are stored in the compressed `data.zip` folder available through the Open Science Framework repository.

- The `code` folder contains all R scripts required to replicate the analyses and results presented in the paper.

- The `results` folder includes .RData files containing the main outputs of the paper, including posterior samples and indirect TFR estimates.

- The `table` folder contains the two main tables of the paper saved in .tex format.
  
- The `figure` folder contains the four main figures of the paper saved in .pdf format.

### Code
After downloading all the data sources, researchers can replicate the figures and tables of the main paper by running the R scripts within the `code` folder.

1. `data_preparetion`: script to transform data from FamiLinx in a suitable format for the main data analysis
2. `upload_mortality_bounds`: script to generate lower and upper bounds for infant mortality
3. `upload_packages`: script to upload the packages needed for the data analysis
4. `helpers`: script to upload R functions needed for the data analysis
5. `data_input_model`: script to create the inputs needed to fit the Bayesian models
6. `model_code`: script with the main Bayesian model of the paper written in Nimble
7. `fit_model_main_paper`: script to run the main Bayesian model
9. `fit_model_no_gen`: script to run the Bayesian model without FamiLinx data
10. `indirect_estimation`: script to run the indirect estimation method
11. `model_tfr_accuracy`: script to compute the accuracy of the TFR estimation methods
12. `model_validation`: script to conduct out-of-sample validation of the main Bayesian model
13. `figure_smr`: script to create the figures of the main paper
    

### Replication

All analyses and computations were carried out on 2023 MacBook Pro with an Apple M5 Pro chip, 24 GB memory. All analyses were originally conducted using R version 4.5.3 and the package versions recorded in the attached session info at the bottom of the README. 
Re-running the pipeline with updated R or package versions, or a different seed, may yield minor numerical differences. These do not affect the paper’s results or conclusions.


## Session Info

```

R version 4.5.3 (2026-03-11)
Platform: aarch64-apple-darwin20
Running under: macOS Tahoe 26.3.1

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: Europe/Rome
tzcode source: internal

attached base packages:
[1] parallel  splines   stats     graphics  grDevices
[6] utils     datasets  methods   base     

other attached packages:
 [1] gt_1.3.0            knitr_1.51         
 [3] ggthemes_5.2.0      patchwork_1.3.2    
 [5] mgcv_1.9-4          nlme_3.1-168       
 [7] MCMCvis_0.16.5      nimbleHMC_0.2.4    
 [9] nimble_1.4.2        rstan_2.32.7       
[11] StanHeaders_2.32.10 corrr_0.4.5        
[13] zoo_1.8-15          viridis_0.6.5      
[15] viridisLite_0.4.3   babynames_1.0.1    
[17] assertthat_0.2.1    LearnBayes_2.15.2  
[19] plotly_4.12.0       hrbrthemes_0.9.2   
[21] kableExtra_1.4.0    cowplot_1.2.0      
[23] readxl_1.4.5        h2o_3.44.0.3       
[25] mapdata_2.3.1       maps_3.4.3         
[27] reshape2_1.4.5      lubridate_1.9.5    
[29] forcats_1.0.1       stringr_1.6.0      
[31] dplyr_1.2.1         purrr_1.2.2        
[33] readr_2.2.0         tidyr_1.3.2        
[35] tibble_3.3.1        ggplot2_4.0.3      
[37] tidyverse_2.0.0     data.table_1.18.2.1
[39] here_1.0.2         

loaded via a namespace (and not attached):
 [1] bitops_1.0-9            gridExtra_2.3          
 [3] inline_0.3.21           rlang_1.2.0            
 [5] magrittr_2.0.5          otel_0.2.0             
 [7] matrixStats_1.5.0       compiler_4.5.3         
 [9] loo_2.9.0               systemfonts_1.3.2      
[11] vctrs_0.7.3             pkgconfig_2.0.3        
[13] fastmap_1.2.0           backports_1.5.1        
[15] labeling_0.4.3          utf8_1.2.6             
[17] rmarkdown_2.31          tzdb_0.5.0             
[19] pracma_2.4.6            xfun_0.57              
[21] jsonlite_2.0.0          splines2_0.5.4         
[23] terra_1.9-11            broom_1.0.12           
[25] R6_2.6.1                stringi_1.8.7          
[27] RColorBrewer_1.1-3      extrafontdb_1.1        
[29] car_3.1-5               cellranger_1.1.0       
[31] numDeriv_2016.8-1.1     Rcpp_1.1.1-1.1         
[33] extrafont_0.20          Matrix_1.7-4           
[35] igraph_2.2.3            timechange_0.4.0       
[37] tidyselect_1.2.1        rstudioapi_0.18.0      
[39] abind_1.4-8             yaml_2.3.12            
[41] codetools_0.2-20        pkgbuild_1.4.8         
[43] lattice_0.22-9          plyr_1.8.9             
[45] withr_3.0.2             S7_0.2.2               
[47] coda_0.19-4.1           evaluate_1.0.5         
[49] survival_3.8-6          RcppParallel_5.1.11-2  
[51] xml2_1.5.2              pillar_1.11.1          
[53] ggpubr_0.6.3            carData_3.0-6          
[55] stats4_4.5.3            generics_0.1.4         
[57] sp_2.2-1                rprojroot_2.1.1        
[59] RCurl_1.98-1.18         hms_1.1.4              
[61] scales_1.4.0            glue_1.8.1             
[63] gdtools_0.5.0           lazyeval_0.2.3         
[65] tools_4.5.3             ggsignif_0.6.4         
[67] fs_2.1.0                grid_4.5.3             
[69] Rttf2pt1_1.3.14         QuickJSR_1.9.2         
[71] colorspace_2.1-2        raster_3.6-32          
[73] Formula_1.2-5           cli_3.6.6              
[75] textshaping_1.0.5       fontBitstreamVera_0.1.1
[77] svglite_2.2.2           gtable_0.3.6           
[79] rstatix_0.7.3           fontquiver_0.2.1       
[81] digest_0.6.39           htmlwidgets_1.6.4      
[83] farver_2.1.2            htmltools_0.5.9        
[85] lifecycle_1.0.5         httr_1.4.8             
[87] fontLiberation_0.1.0   
