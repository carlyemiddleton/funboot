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






























}
