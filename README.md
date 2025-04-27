# FunBoot

Carly Middleton

## Installation

First, install and load the dependencies:

```
library(spatstat.data)
library(spatstat.univar)
library(spatstat.geom)
library(spatstat.random)
library(spatstat.explore)
library(refund)
library(dplyr)
library(parallel)
```

For the spatstat package, it is recommended to load only the above spatstat modules instead of loading library(spatstat) all together.  In our experience, some of the spatstat modules besides the ones above have interfered with with the pffr() function.

## How to load the example data 

```
library(devtools)
devtools::install_github('carlyemiddleton/funboot')
library(funboot)
data(data_example)
```

Vignette things to cover:  
1. preprocess_data(), and demonstration of the permuted outcome being corrected
2. make a spaghetti plot function for the summary functions?
3. make a plot.images() function for color coded cell image plots?
4. fit_model()
5. multiplicity adjustment
6. wildBS_CB() and plot.wildBS_CB()
7. Ftest() 


Should make a note about the required variable naming format for funboot

Also make a not about pffr() having 2 intercepts and how to calculate them together

make a note that the vignette takes about 20-30minutes to build on a laptop computer with 8 cores and 2.3 GHz processor