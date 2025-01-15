# Load packages
library(spicyR)

NOWEIGHTS <- NULL
save(NOWEIGHTS, file='NOWEIGHTS.RData')

## INITIALISE
s1 <- Sys.time()

seed = 456
set.seed(seed)
#window <- owin(xrange = c(0, 1000),
#               yrange = c(0, 1000))
nPatients <- 100
nIm <- 1
nSim <- 2
nCores <- 50
counts <- seq(from = 100, to = 400, by = 10)
Rs <- seq(from = 0, to = 200, by = 1)

#lambda = 60 


## SIGNAL
print('starting spicy()')
Sys.time()

sim  <- function(i, counts, nPatients, nIm, window, lambda){
  
  set.seed(i)
  
  load(file=paste0(getwd(),'/environments/myEnvironment',i,'.RData'))
  
  test.NoWeights <- spicyR::spicy(cellExp2,
                                  condition = "condition",
                                  subject = "subject",
                                  from = "B",
                                  to = "A",
                                  fast = TRUE, Rs = Rs,
                                  edgeCorrect = TRUE,
                                  weights = FALSE,
                                  verbose = FALSE,
                                  BPPARAM = parallel_backend)
  load(file='NOWEIGHTS.RData')
  NOWEIGHTS <- c(NOWEIGHTS, test.NoWeights$p.value[1,'conditionGroup2'])
  save(NOWEIGHTS, file='NOWEIGHTS.RData')
  #save(test.NoWeights, file=paste0('test.NoWeights',i,'.RData'))
}

lapply(as.list(seq_len(1000)+seed),sim)

print('finished spicy()')
Sys.time()
  