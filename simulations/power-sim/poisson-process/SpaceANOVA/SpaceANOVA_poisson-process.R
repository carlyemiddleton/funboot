library(spatstat)
library(tidyr)
library(dplyr)
library(fda.usc)
library(ggplot2)
library(cowplot)
library(gridExtra)
library(SpaceANOVA)
library(spatstat.geom)
library(spatstat.random)
library(spatstat.explore)
library(stringr)
library(parallel)

#Setup
seed <- 12345

counts <- seq(from = 100, to = 400, by = 10)
nPatients <- 100 #number of patients per dataset
nIm <- 1 #number of images per patient

nSim <- 500 #number of datasets

n_cores <- parallel::detectCores(logical=F)
cl <- makeCluster(n_cores-1)

sim  <- function(i, counts, nPatients, nIm, sigma, delta, cl){

  set.seed(i)

  #Generate data from Poisson process
  #Same setup as in https://github.com/nickcee/spicyRPaper/blob/main/baseSim.R:
  window <- spatstat.geom::owin(xrange = c(0, 1000), yrange = c(0, 1000))
  g1 <- rpois(nPatients/2, sigma) #g1 is the control group
  g2 <- rpois(nPatients/2, sigma + delta) #g2 is the treatment group
  adjustSigma = c(g1,g2)+1

  x <- y <- cellType <- imageID <- NULL
  for (p in 1:nPatients) {
    for (j in 1:nIm) {
      sCount1 <- sample(counts,1)
      sCount2 <- sample(counts,1)
      a <- spatstat.random::rpoispp(sCount1/1000^2, win = window)
      aDens <- spatstat.explore::density.ppp(a, sigma = adjustSigma[p], kernel = "disc")
      aDens$v <- pmax(aDens$v,0)*sCount2/sCount1
      b <- spatstat.random::rpoispp(aDens)
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

  phenoData <- data.frame(
    imageID = unique(imageID),
    condition = rep(c("Group1", "Group2"), each = nIm*nPatients/2),
    subject = rep(1:nPatients, each = nIm)
    )

  table <- cellExp
  names(table) <- c('Xcoord','Ycoord','cell.type','image')
  table$patient_id <- sapply(table$image, function(x){str_split(x, fixed("_"))[[1]][1]})
  table$roi <- sapply(table$image, function(x){str_split(x, fixed("_"))[[1]][2]})
  names(phenoData) <- c('image','group','patient_id')
  table <- merge(table, phenoData, by=c('image','patient_id'), all=T)
  colnames(table) <- c('imageID','ID','x', 'y', "cellType",  'roi', 'Group')

  #Apply SpaceANOVA
  Final_result <- SpaceANOVA::All_in_one(
                              data = table,
                              fixed_r = seq(0, 200, by = 1),
                              Summary_function = "L",
                              Hard_ths = 0,
                              homogeneous = T,
                              interaction_adjustment = T,
                              perm = F,
                              cores = n_cores)

  p_res <- p_extract(Final_result)
  Univ_p <- p_res[[1]][1,2] #select the p-value for B cells around A cells
  Mult_p <- p_res[[2]][1,2]

  c(Univ_p = Univ_p, Mult_p = Mult_p)

}

#Simulate nSim datasets for each sigma and delta
sigma_vals <- c(20, 40, 60, 80, 100)
for (sigma_val in sigma_vals) {
  delta_vals <- c(0, sigma_val/10, sigma_val/5, sigma_val/2)
  for (delta_val in delta_vals) {
    Sys.time()
    res <- lapply(
      as.list(seq_len(nSim) + seed),
      sim,
      counts = counts,
      nPatients = nPatients,
      nIm = nIm,
      sigma = sigma_val,
      delta = delta_val,
      cl=cl
    )
    Sys.time()
    res <- do.call("rbind", res)
    save(res, file = paste0("res_sigma", sigma_val, "_delta", delta_val, ".RData"))
  }
}




