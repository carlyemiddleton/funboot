# funboot

[![DOI](https://zenodo.org/badge/1154788071.svg)](https://doi.org/10.5281/zenodo.21403975)

<img src="man/figures/logo.png" align="right" height="139" alt="" /> 

## Installation

`funboot` uses version `0.1-37` of `refund`.  Install that using the following command:

```
devtools::install_version(
  "refund",
  version = "0.1-37",
  repos = "https://cloud.r-project.org",
  upgrade = "never"
)
```

Then, `funboot` can be installed as follows:

```
devtools::install_github("carlyemiddleton/funboot", upgrade = "never")
```

Load the package and its dependencies:

```
library(ggplot2)
library(spatstat.geom)
library(spatstat.explore)
library(refund)
library(mgcv)
library(funboot)
```

It is recommended to load only `spatstat.geom` and `spatstat.explore` instead of loading `spatstat` all together, as other packages in the `spatstat` family may interfere with `refund::pffr().`



## Required Format for Input Data

The package requires single-cell input data in the form of an R `data.frame` containing variables named exactly as: `patient_id`, `image_number`, `cell_id`, `cell_x`, `cell_y`, `cell_type`, and any covariates.  Two example datasets following this format are included in the package:  

```
data("breastcancer_data", package = "funboot")
data("melanoma_data", package = "funboot")
```

## Vignette

The package contains a vignette, which can be built during the installation:

```
detach("package:funboot")
devtools::install_github("carlyemiddleton/funboot", upgrade = "never", build_vignettes = TRUE)
browseVignettes("funboot")
```
Alternatively, a pre-rendered version of it is accessible [here](https://htmlpreview.github.io/?https://github.com/carlyemiddleton/funboot/blob/main/vignettes/funboot.html).


## Reproducibility

Code for reproducing all study results is located in the [`case-studies/`](case-studies/) and [`simulations/`](simulations/) folders.


## Session Info 

```
R version 4.4.0 (2024-04-24 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 11 x64 (build 26200)

Matrix products: default


locale:
[1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8   
[3] LC_MONETARY=English_United States.utf8 LC_NUMERIC=C                          
[5] LC_TIME=English_United States.utf8    

time zone: America/New_York
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] mgcv_1.9-3             refund_0.1-37          spatstat.explore_3.4-3 nlme_3.1-164          
 [5] spatstat.random_3.4-1  spatstat.geom_3.4-1    spatstat.univar_3.1-3  spatstat.data_3.1-6   
 [9] ggplot2_4.0.1          funboot_0.1.0         

loaded via a namespace (and not attached):
 [1] tidyselect_1.2.1      hdrcde_3.4            RLRsim_3.1-8          dplyr_1.1.4           farver_2.1.2         
 [6] S7_0.2.1              bitops_1.0-9          fastmap_1.2.0         RCurl_1.98-1.17       pracma_2.4.4         
[11] promises_1.3.2        digest_0.6.37         mime_0.13             lifecycle_1.0.4       cluster_2.1.6        
[16] ellipsis_0.3.2        processx_3.8.6        magrittr_2.0.3        compiler_4.4.0        rlang_1.1.6          
[21] tools_4.4.0           grpreg_3.5.0          htmlwidgets_1.6.4     pkgbuild_1.4.7        mclust_6.1.1         
[26] curl_6.2.3            RColorBrewer_1.1-3    rainbow_3.8           pkgload_1.4.0         abind_1.4-8          
[31] KernSmooth_2.23-22    gamm4_0.2-7           fda_6.3.0             miniUI_0.1.2          withr_3.0.2          
[36] purrr_1.0.4           desc_1.4.3            polyclip_1.10-7       grid_4.4.0            pcaPP_2.0-5          
[41] urlchecker_1.0.1      profvis_0.4.0         xtable_1.8-4          colorspace_2.1-1      spatstat.utils_3.1-4 
[46] scales_1.4.0          MASS_7.3-60.2         cli_3.6.5             mvtnorm_1.3-3         reformulas_0.4.1     
[51] generics_0.1.4        remotes_2.5.0         rstudioapi_0.17.1     magic_1.6-1           sessioninfo_1.2.3    
[56] fds_1.8               minqa_1.2.8           cachem_1.1.0          splines_4.4.0         parallel_4.4.0       
[61] vctrs_0.6.5           devtools_2.4.5        boot_1.3-30           Matrix_1.7-0          callr_3.7.6          
[66] tensor_1.5            goftest_1.2-3         glue_1.8.0            nloptr_2.2.1          ps_1.9.1             
[71] gtable_0.3.6          deldir_2.0-4          later_1.4.2           lme4_1.1-37           tibble_3.2.1         
[76] pillar_1.10.2         htmltools_0.5.8.1     deSolve_1.40          R6_2.6.1              Rdpack_2.6.4         
[81] ks_1.15.1             shiny_1.10.0          lattice_0.22-6        rbibutils_2.3         memoise_2.0.1        
[86] httpuv_1.6.16         pbs_1.1               Rcpp_1.0.14           spatstat.sparse_3.1-0 fs_1.6.6             
[91] usethis_3.1.0         pkgconfig_2.0.3        
```

