# funboot

<img src="man/figures/logo.png" align="right" height="139" alt="" /> 

## Installation

First, load all package dependencies:

```
library(ggplot2)
library(spatstat.geom)
library(spatstat.explore)
library(refund)
library(mgcv)
```

It is recommended to load only `spatstat.geom` and `spatstat.explore` instead of loading `spatstat` all together, as other packages in the `spatstat` family may interfere with `refund::pffr().`

Once all dependencies are present, `funboot` can be installed and loaded as follows:

```
devtools::install_github('carlyemiddleton/funboot')
library(funboot)
```

## Required Format for Input Data

The package requires single-cell input data in the form of an R `data.frame` containing variables named exactly as: `patient_id`, `image_number`, `cell_id`, `cell_x`, `cell_y`, `cell_type`, and any covariates.  Two example datasets following this format are included in the package:  

```
data("breastcancer_data", package="funboot")
data("melanoma_data", package="funboot")
```

## Vignette

The package contains a vignette, which can be built during the install:

```
devtools::install_github('carlyemiddleton/funboot', build_vignettes = TRUE, force=TRUE)
browseVignettes("funboot")
```
Alternatively, a pre-rendered verison of it is accessible [here](https://htmlpreview.github.io/?https://github.com/carlyemiddleton/funboot/blob/main/vignettes/funboot.html).

The vignette takes about 4 minutes to build on a standard laptop computer with 8 cores and a 2.3 GHz processor.  

