library(spatstat.geom)
library(spatstat.random)
library(stringr)

seed <- 12345

nPatients <- 50
nIm <- 3

sigmaB <- 2        #Change to 3,4,5,6
sigmaB_diff <- 0   #Change to 1/3*sigmaB, 1/3*sigmaB, sigmaB

if (dir.exists("environments")) {
  unlink("environments", recursive = TRUE)
}

iters <- 1:175 #Change to 176:350 and 351:500 for full results. Running in chunks due to cluster timeout

for(iter in iters){
  set.seed(iter+seed)

  #Generate data
  x <- y <- cellType <- imageID <- NULL
  sigmaB_G1 <- sigmaB
  sigmaB_G2 <- sigmaB + sigmaB_diff
  for (p in 1:nPatients){
    group <- ifelse(p <= nPatients/2, "Group1", "Group2")
    if (group == "Group1") {
      sigmaB_patient  <- abs(rnorm(1, mean = sigmaB_G1, sd = sqrt(sigmaB_G1 / 10)))
    } else {
      sigmaB_patient  <- abs(rnorm(1, mean = sigmaB_G2, sd = sqrt(sigmaB_G2 / 10)))
    }
    for (j in 1:nIm){
      Rhard  <- sample(c(0.02, 0.03,  0.04,  0.05, 0.06), 1)
      beta_p <- sample(seq(50, 100, by=10), 1)
      muB    <- sample(seq(10, 30, by=1), 1)
      A_unit <- spatstat.random::rHardcore(
        beta = beta_p,
        R    = Rhard,
        W    = spatstat.geom::owin(c(0,1),c(0,1))
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
              rnorm(kB, A_x[i], sigmaB_patient),
              rnorm(kB, A_y[i], sigmaB_patient)
            )
            B_xy <- rbind(B_xy, ptsB)
          }
        }
      }
      #Delete type B cells generated outside of the image window
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
    condition = ifelse(
      as.numeric(sapply(unique(imageID),function(x){
        stringr::str_split(x, "_")[[1]][1]})) <= nPatients/2,
      "Group1",
      "Group2"
      ),
    subject = as.numeric(sapply(unique(imageID),function(x){
      stringr::str_split(x, "_")[[1]][1]
      }))
    )

  cellExp2 <- merge(cellExp, phenoData, by='imageID')

  if(!dir.exists('environments')){
    dir.create('environments')
  }
  save(cellExp2, file = paste0("environments/cellExp2_", iter, ".RData"))
}
