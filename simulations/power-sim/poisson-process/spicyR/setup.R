library(spatstat.geom)
library(spatstat.random)
library(spatstat.explore)
library(stringr)

seed <- 12345

nPatients <- 100 #number of patients per dataset
nIm <- 1 #number of images per patient

sigma <- 20   #Also change to 40,60,80,100
delta <- 0    #Also change to sigma/10, sigma/5, sigma/2

counts <- seq(100, 400, by = 10)

if (dir.exists("environments")) {
  unlink("environments", recursive = TRUE)
}

iters <- 1:500 #number of datasets

for(iter in iters){
  set.seed(iter+seed)

  #Generate data from Poisson process
  #Same setup as in https://github.com/nickcee/spicyRPaper/blob/main/baseSim.R:
  window <- spatstat.geom::owin(xrange = c(0, 1000), yrange = c(0, 1000))
  g1 <- rpois(nPatients / 2, sigma) #g1 is the control group
  g2 <- rpois(nPatients / 2, sigma + delta) #g2 is the treatment group
  adjustSigma <- c(g1, g2) + 1

  x <- y <- cellType <- imageID <- NULL
  for (p in 1:nPatients){
    for (j in 1:nIm){
      sCount1 <- sample(counts, 1)
      sCount2 <- sample(counts, 1)
      a <- spatstat.random::rpoispp(sCount1/1000^2, win = window)
      aDens <- spatstat.explore::density.ppp(a, sigma = adjustSigma[p], kernel = "disc")
      aDens$v <- pmax(aDens$v, 0) * sCount2 / sCount1
      b <- spatstat.random::rpoispp(aDens)
      x <- c(x, a$x, b$x)
      y <- c(y, a$y, b$y)
      cellType <- c(cellType, rep("A", a$n), rep("B", b$n))
      imageID  <- c(imageID, rep(paste(p, j, sep = "_"), a$n + b$n))
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
