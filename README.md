# funboot

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

Then, `funboot` can be installed and loaded as follows:

```
devtools::install_github('carlyemiddleton/funboot', upgrade = 'never')
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
data("breastcancer_data", package="funboot")
data("melanoma_data", package="funboot")
```

## Vignette

The package contains a vignette, which can be built during the installation:

```
devtools::install_github('carlyemiddleton/funboot', upgrade='never', build_vignettes = TRUE)
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

loaded via a namespace (and not attached):
 [1] magic_1.6-1        ks_1.15.1          generics_0.1.4     bitops_1.0-9       KernSmooth_2.23-22
 [6] lattice_0.22-6     hdrcde_3.4         pracma_2.4.4       lme4_1.1-37        magrittr_2.0.3    
[11] gamm4_0.2-7        grid_4.4.0         RColorBrewer_1.1-3 mvtnorm_1.3-3      Matrix_1.7-0      
[16] deSolve_1.40       mclust_6.1.1       mgcv_1.9-3         scales_1.4.0       abind_1.4-8       
[21] reformulas_0.4.1   Rdpack_2.6.4       cli_3.6.5          rlang_1.1.6        fda_6.3.0         
[26] rbibutils_2.3      refund_0.1-37      splines_4.4.0      RLRsim_3.1-8       parallel_4.4.0    
[31] tools_4.4.0        nloptr_2.2.1       funboot_0.1.0      minqa_1.2.8        dplyr_1.1.4       
[36] colorspace_2.1-1   ggplot2_4.0.1      rainbow_3.8        boot_1.3-30        vctrs_0.6.5       
[41] R6_2.6.1           lifecycle_1.0.4    fds_1.8            MASS_7.3-60.2      pcaPP_2.0-5       
[46] cluster_2.1.6      pkgconfig_2.0.3    grpreg_3.5.0       pillar_1.10.2      gtable_0.3.6      
[51] glue_1.8.0         Rcpp_1.0.14        pbs_1.1            tibble_3.2.1       tidyselect_1.2.1  
[56] rstudioapi_0.17.1  farver_2.1.2       nlme_3.1-164       compiler_4.4.0     S7_0.2.1          
[61] RCurl_1.98-1.17   
```

