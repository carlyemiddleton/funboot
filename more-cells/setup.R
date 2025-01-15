# Load packages
library(spatstat)
library(BiocParallel)
library(parallel)
#library(tidyverse)
#library(plotROC)
#library(spicyR)
library(stringr)
library(tvem)

ifelse(dir.exists('environments'),
       unlink('environments', recursive = T),
       "Directory doesn't exist")

res <- NULL
save(res, file='res.RData')

## INITIALISE
s1 <- Sys.time()

seed = 456
set.seed(seed)
window <- owin(xrange = c(0, 1000),
               yrange = c(0, 1000))
nPatients <- 100
nIm <- 1
nSim <- 2
nCores <- 50
counts <- seq(from = 100, to = 400, by = 10)
Rs <- seq(from = 0, to = 200, by = 1)

lambda = 100 


## SIGNAL
s1 <- Sys.time()

#sim  <- function(i, counts, nPatients, nIm, window, lambda){
for(i in seq_len(1000)+seed){
  set.seed(i)
  
  g1 <- rpois(nPatients/2, lambda)
  g2 <- rpois(nPatients/2, lambda + 0 )
  adjustSigma = c(g1,g2)+1
  
  x <- c()
  y <- c()
  cellType <- c()
  imageID <- c()
  
  for (p in 1:nPatients) {
    for (j in 1:nIm) {
      sCount1 <- sample(counts,1)
      sCount2 <- sample(counts,1)
      a <- rpoispp(sCount1/1000^2, win = window)
      aDens <- density(a, sigma = adjustSigma[p], kernel = "disc")
      aDens$v <- pmax(aDens$v,0)*sCount2/sCount1
      b <- rpoispp(aDens)
      
      x <- c(x, a$x, b$x)
      y <- c(y, a$y, b$y)
      
      cellType <- c(cellType, rep("A", a$n), rep("B", b$n))
      imageID <- c(imageID, rep(paste(p,j,sep = "_"), a$n+b$n))
    }
  }
  
  imageID <- factor(imageID)
  
  cellExp <- data.frame(
    x = x,
    y = y,
    cellType = factor(cellType),
    imageID = imageID
  )
  #cellExp2 <- spicyR::SegmentedCells(cellExp, verbose = FALSE)
  
  
  phenoData <- data.frame(imageID = unique(imageID),
                          condition = rep(c("Group1", "Group2"), each = nIm*nPatients/2),
                          subject = rep(1:nPatients, each = nIm))
  #spicyR::imagePheno(cellExp2) <- phenoData
  cellExp2 <- merge(cellExp, phenoData, by='imageID')
  
  numcores <- detectCores() - 1
  parallel_backend <- MulticoreParam(workers = numcores) 
  
  ifelse(!dir.exists('environments'),
         dir.create('environments'),
         "Directory Exists")
  save.image(file=paste0('environments/myEnvironment',i,'.RData'))
  
}
  