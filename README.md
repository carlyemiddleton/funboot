# FunBoot

<img src="man/figures/logo.png" align="right" height="139" alt="" /> 

## Installation

First, install and load the dependencies.

```
#Dependencies for the analysis functions
library(spatstat.data)
library(spatstat.univar)
library(spatstat.geom)
library(spatstat.random)
library(spatstat.explore)
library(refund)
library(dplyr)
library(parallel)

#Dependencies for the plotting functions
library(ggplot2)
library(patchwork)

```

For the spatstat package, it is recommended to load only the above spatstat modules instead of loading library(spatstat) all together.  In our experience, some of the spatstat modules besides the ones above have interfered with refund::pffr().

Once all dependencies are present, install the package:

```
devtools::install_github('carlyemiddleton/funboot')
library(funboot)
```

## Example Data

Our package requires input data in the form of an R data.frame containing variables named exactly as: **patient_id**, **image_number**, **cell_id**, **cell_x**, **cell_y**, **cell_type**, and any covariates to be adjusted for.  An example dataset following this format is built into the package:  

```
data('breastcancer_data_subset_100images.rda', package='funboot')
```

## Package Vignette

This package contains a vignette, which can be built during the install:

```
devtools::install_github('carlyemiddleton/funboot', build_vignettes = TRUE, force=TRUE)
browseVignettes("funboot")
```

The vignette takes about 30 minutes to build on a laptop computer with 8 cores and a 2.3 GHz processor. 

## Two Intercepts

Note: pffr() estimates two different intercepts:  one constant intercept, and one radius-varying intercept.  The intercept we refer to in our manuscript, $\beta_0(r)$, manuscript is the sum of these two intercepts, and its confidence band can be estimated via funboot::EY_CB() with argument covar.list=list().  funboot::wildBS_CB() produces separate confidence bands for each intercept. 
