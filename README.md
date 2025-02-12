# Loca-Cola

Carly Middleton

## Installation

Loca-Cola has several dependencies.  First, install and load the dependencies:

```
library(spatstat.data)
library(spatstat.univar)
library(spatstat.geom)
library(spatstat.random)
library(spatstat.explore)
library(refund)
```

For the spatstat package, it is recommended to load only the above spatstat modules instead of loading library(spatstat) all together.  In our experience, some of the spatstat modules besides the ones above have interfered with with the pffr() function.

## How to load the example data 

```
library(devtools)
devtools::install_github('carlyemiddleton/phantem')
library(phantem)
data(data_example)
```
