library(BiocParallel)
library(spicyR)

seed <- 12345

Rs <- seq(from = 0, to = 200, by = 1) #evaluation radii
parallel_backend <- BiocParallel::MulticoreParam(workers = 24)

iters <- 1:500

print('starting spicy()')
Sys.time()

sim  <- function(iter, Rs=Rs){

  set.seed(iter+seed)

  load(file=paste0(getwd(),'/environments/cellExp2_',iter,'.RData'))

  test_no_weights <- spicyR::spicy(cellExp2,
                                   condition = "condition",
                                   subject = "subject",
                                   from = "A",
                                   to = "B",
                                   fast = TRUE,
                                   Rs = Rs,
                                   edgeCorrect = TRUE,
                                   weights = FALSE,
                                   verbose = FALSE,
                                   BPPARAM = parallel_backend)

  test_no_weights$p.value[1, "conditionGroup2"]
}

pvals <- vapply(iters, sim, FUN.VALUE = numeric(1), Rs = Rs)
names(pvals) <- iters

saveRDS(list(iters = iters, Rs = Rs, pvals = pvals), file = "NOWEIGHTS.rds")

print('finished spicy()')
Sys.time()
