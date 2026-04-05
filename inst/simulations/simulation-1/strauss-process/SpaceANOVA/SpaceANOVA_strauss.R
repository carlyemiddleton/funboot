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
nPatients <- 50
nIm <- 3

nSim <- 500

n_cores <- parallel::detectCores(logical=F)
cl <- makeCluster(n_cores-1)

sim  <- function(i, counts, nPatients, nIm, tauB, tauB_diff, cl){
  
  set.seed(i)
  
  ##Generate data from Strauss process
  x <- y <- cellType <- imageID <- NULL
  tauB_G1 <- tauB
  tauB_G2 <- tauB + tauB_diff
  for (p in 1:nPatients){
    group <- ifelse(p <= nPatients/2, "Group1", "Group2")
    if (group == "Group1") {
      tauB  <- abs(rnorm(1, mean = tauB_G1, sd = sqrt(tauB_G1 / 10)))
    } else {
      tauB  <- abs(rnorm(1, mean = tauB_G2, sd = sqrt(tauB_G2 / 10)))
    }
    nIm <- nIm
    for (j in 1:nIm){
      Rhard  <- sample(c(0.02, 0.03,  0.04,  0.05, 0.06), 1)
      beta_p <- sample(seq(50, 100, by=10), 1)
      muB    <- sample(seq(10, 30, by=1), 1)
      A_unit <- spatstat.random::rHardcore(
        beta = beta_p,
        R    = Rhard,
        W    = owin(c(0,1),c(0,1))
      )
      A_x <- A_unit$x * 1000
      A_y <- A_unit$y * 1000
      nA  <- A_unit$n
      B_xy <- matrix(numeric(0), ncol = 2)
      if (nA > 0){
        for (i in 1:nA){
          kB <- rpois(1, muB)
          if (kB > 0) {
            ptsB <- cbind(
              rnorm(kB, A_x[i], tauB),
              rnorm(kB, A_y[i], tauB)
            )
            B_xy <- rbind(B_xy, ptsB)
          }
        }
      }
      #remove B cells generated outside of the window
      if (nrow(B_xy) > 0) {
        inside <- B_xy[,1] >= 0 & B_xy[,1] <= 1000 &
          B_xy[,2] >= 0 & B_xy[,2] <= 1000
        B_xy <- B_xy[inside, , drop = FALSE]
      }
      ID <- paste(p, j, sep = "_")
      if (nA > 0) {
        x <- c(x, A_x)
        y <- c(y, A_y)
        cellType <- c(cellType, rep("A", nA))
        imageID  <- c(imageID, rep(ID, nA))
      }
      if (nrow(B_xy) > 0) {
        x <- c(x, B_xy[,1])
        y <- c(y, B_xy[,2])
        cellType <- c(cellType, rep("B", nrow(B_xy)))
        imageID  <- c(imageID, rep(ID, nrow(B_xy)))
      }
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
    imageID  = unique(imageID),
    condition = ifelse(as.numeric(sapply(unique(imageID), function(x) str_split(x, "_")[[1]][1]))
                       <= nPatients/2, "Group1", "Group2"),
    subject   = as.numeric(sapply(unique(imageID), function(x) str_split(x, "_")[[1]][1]))
  )
  table <- cellExp
  names(table) <- c('Xcoord','Ycoord','cell.type','image') 
  table$patient.id <- sapply(table$image, function(x) str_split(x, "_")[[1]][1])
  table$roi        <- sapply(table$image, function(x) str_split(x, "_")[[1]][2])
  names(phenoData) <- c('image','group','patient.id')
  table <- merge(table, phenoData, by=c('image','patient.id'), all=T)
  colnames(table) <- c('imageID','ID','x', 'y', "cellType",  'roi', 'Group')  #variable names as in SpaceANOVA's readme
  
  Final_result = SpaceANOVA::All_in_one(data = table, 
                                        fixed_r = seq(0, 200, by = 1), 
                                        Summary_function = "L", 
                                        Hard_ths = 0, 
                                        homogeneous = FALSE, 
                                        interaction_adjustment = T, 
                                        perm = F, 
                                        cores = n_cores)
  p_res = p_extract(Final_result)
  Mult_p = p_res[[2]][1,2]
  
  c(Mult_p = Mult_p)
  
}

#Simulate for each tauB and tauB_diff
tauB_vals <- c(2, 3, 4, 5, 6)
for (tauB_val in tauB_vals) {
  tauB_diff_vals <- c(0, tauB_val*1/3, tauB_val*2/3, tauB_val)
  for (tauB_diff_val in tauB_diff_vals) {
    Sys.time()
    res <- lapply(
      as.list(seq_len(nSim) + seed),
      sim,
      counts = counts,
      nPatients = nPatients,
      nIm = nIm,
      tauB = tauB_val,
      tauB_diff = tauB_diff_val,
      cl=cl
    )
    Sys.time()
    res <- do.call("rbind", res)
    save(res, file = paste0("res_tauB", tauB_val, "_tauB_diff", tauB_diff_val, ".RData"))
  }
}




