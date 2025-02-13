#' @title Preprocesses data
#'
#' @param sumfun.data is
#' @return It returns res.wildBS
#'
#' @export

test_wildBS <- function(sumfun.data, image.dims, seed=456, parallel=T){
  set.seed(seed)
  image.xmax <- image.dims[2]
  image.ymax <- image.dims[4]
  image.xmin <- image.dims[1]
  image.ymin <- image.dims[3]
  W <- spatstat.geom::owin(c(image.xmin,image.xmax), c(image.ymin,image.ymax))
  #nPatients <- length(unique(sumfun.data$patient_id))

  counts <- seq(from = 100, to = 400, by = 10)
  Rs <- seq(from = 0, to = 200, by = 1)

  lambda = 100


  #put in the documentation:  this package uses one constant and one functional intercept, and does not allow suppression of either
  # beta0_constant <- rep(summary(pffrmodel)$p.coeff, length(grid))
  # beta0_hat <- coef(pffrmodel,n1=length(grid))$smterms$`Intercept(grid)`$coef$value
  # for(i in 1:length(formula)){
  #   list <- coef(pffrmodel,n1=length(grid))$smterms[
  #     which(gsub('.{6}$', '', names(coef(pffrmodel,n1=length(grid))$smterms))==all.vars(formula)[i+1])]
  #   assign(paste0('beta_hat_',all.vars(formula)[i+1]),
  #          list[[1]]$coef$value)
  # }




























}
