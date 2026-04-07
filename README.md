# funboot

<img src="man/figures/logo.png" align="right" height="139" alt="" /> 

## Installation

The package's dependencies are below.

```
library(ggplot2)
library(spatstat.geom)
library(spatstat.explore)
library(refund)
library(mgcv)

```

For the spatstat package, it is recommended to load only the above spatstat modules instead of loading library(spatstat) all together.  In our experience, some of the spatstat modules besides the ones above have interfered with `refund::pffr().`

Once all dependencies are present, `funboot` can be installed and loaded as follows:

```
devtools::install_github('carlyemiddleton/funboot', build_vignettes = TRUE, force=TRUE)
```

```
library(funboot)
```

```
devtools::install_github('carlyemiddleton/funboot')
library(funboot)
```

## Example Data

Our package requires input data in the form of an R data.frame containing variables named exactly as: `patient_id`, `image_number`, `cell_id`, `cell_x`, `cell_y`, `cell_type`, and any covariates to be adjusted for.  Two example datasets following this format are included in the package:  

```
data("breastcancer_data", package="funboot")
data("melanoma_data", package="funboot")
```

## Vignette

This package contains a vignette, which can be built during the install:

```
devtools::install_github('carlyemiddleton/funboot', build_vignettes = TRUE, force=TRUE)
browseVignettes("funboot")
```
Alternatively, a pre-rendered verison of it is accessible [here](https://htmlpreview.github.io/?https://github.com/carlyemiddleton/funboot/blob/main/vignettes/funboot.html).

The vignette takes about 4 minutes to build on a laptop computer with 8 cores, a 2.3 GHz processor, and 16 GB of RAM.  

