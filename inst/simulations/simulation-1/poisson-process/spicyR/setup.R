library(spatstat.geom)
library(spatstat.random)
library(spatstat.explore)
library(stringr)

seed <- 12345

nPatients <- 100
nIm <- 1

sigma <- 20   #Change to 40,60,80,100
delta <- 0   #Change to sigma/10, sigma/5, sigma/2

counts <- seq(100, 400, by = 10)

if (dir.exists("environments")) {
  unlink("environments", recursive = TRUE)
}

iters <- 1:500

for(iter in iters){
  set.seed(iter+seed)

  #Generate data from Poisson process
  window <- spatstat.geom::owin(xrange = c(0, 1000), yrange = c(0, 1000))
  g1 <- rpois(nPatients / 2, sigma)
  g2 <- rpois(nPatients / 2, sigma + delta)
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
