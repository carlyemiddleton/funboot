# FunBoot

Package implementing the proposed method from the manuscript "FunBoot: Pairwise cell type colocalization analysis for multiplex imaging data via functional regression and
wild bootstrap."  

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

For the spatstat package, it is recommended to load only the above spatstat modules instead of loading library(spatstat) all together.  In our experience, some of the spatstat modules besides the ones above can interfere with with the function pffr() from the package *refund*, which we use to fit our method's functional regression model.

Once all dependencies are present, install the package:

```
devtools::install_github('carlyemiddleton/funboot')
library(funboot)
```

## Example Data

Our package requires input data in the form of an R data.frame containing variables named exactly as: **patient_id**, **image_number**, **cell_id**, **cell_x**, **cell_y**, **cell_type**, and any covariates to be adjusted for.  An example dataset following this format is built into the package:  

```
data('breastcancer_data_subset_100images.rda')
```

## Package Vignette

This package contains a vignette, which can be accessed via:

```
devtools::build_vignettes()
```

The vignette takes about 30 minutes to build on a laptop computer with 8 cores and a 2.3 GHz processor. 

## Two Intercepts

Note: pffr() estimates two different intercepts:  one constant intercept, and one radius-varying intercept.  The $\beta_0(r)$ in our manuscript is the sum of these two intercepts, and its confidence band can be estimated via funboot::lin_comb_CB() with argument lin.comb=list().  funboot::wildBS_CB() produces separate confidence bands for each intercept. 
