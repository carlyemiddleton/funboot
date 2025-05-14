#' @title Preprocesses raw data
#'
#' @description This function preprocess raw data by: 1.) filtering out sparse images using **qc.cellcount.cutoff**, 2.) calculating the desired spatial summary function, 3.) calculating the outcome variable to be used in *pffr()*
#'
#' @param data Data frame containing the variables **patient_id**, **image_number**, **cell_id**, **cell_x**, **cell_y**, **cell_type**, and any covariates to be adjusted for, as in the example data.
#' @param from.cell The "from" cell type to be used in the calculation of the spatial summary function
#' @param to.cell The "to" cell type to be used in the calculation of the spatial summary function
#' @param qc.cellcount.cutoff If an image has **qc.cellcount.cutoff** or less "to" cells, or **qc.cellcount.cutoff** or less "from" cells, the image is excluded from the calculation of the spatial summary function.  (Analogous to SpaceANOVA's **Hard_ths** parameter.)
#' @param P The number of permutations to use in the simulation envelope procedure, as long as **perm.yn**=T.
#' @param perm.yn If TRUE, the permuted mean of the spatial summary function is calculated using the permutation envelope procedure described in the [Wilson paper](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009900).  The outcome variable to be used in *pffr()* is then calculated by subtracting the value of the observed spatial summary function from its permuted mean.  If *perm.yn*=F, the outcome variable is calculated by subtracting the observed spatial summary function from its expected value assuming complete spatial randomness.
#' @param R The maximum radius at which the observed spatial summary function should be evaluated.
#' @param inc The increment of units by which the spatial summary function should be calculated.  For example, if **R**=200 and **inc**=1, the function is evaluated on the grid 0:200.
#' @param image.dims Vector containing the dimensions of the images on which cell locations are recorded, structured as c(x.min, x.max, y.min, y.max).  All images must be of the same dimension.
#' @param summary.function The spatial summary function to be calculated.  Must be one of 'K', 'L', or 'g'.
#' @param seed Random seed to be used during the permutation envelope procedure (optional).
#'
#' @return A data frame containing values of the observed and expected spatial summary function, as well as values of the covariates from the input dataset. L.pmean = permuted mean from the permutation envelope, outcome = L.expect - L.obs if **perm.yn**=F, or L.expect - L.pmean if **perm.yn**=T.
#'
#' @export

preprocess_data <- function(data, from.cell, to.cell, qc.cellcount.cutoff=20, P=50, perm.yn=F,
                            R=200, inc=1, image.dims, summary.function='L',seed=NULL){
  if(!is.null(seed)){set.seed(seed)}
  image.xmax <- image.dims[2]
  image.ymax <- image.dims[4]
  image.xmin <- image.dims[1]
  image.ymin <- image.dims[3]
  W <- spatstat.geom::owin(c(image.xmin,image.xmax), c(image.ymin,image.ymax))

  #########################
  ##Calculate sumfun.data #
  #########################

  if(!(summary.function %in% c('K', 'L', 'g'))){stop("summary.function must be one of 'K', 'L', or g")}

  Kdata <- K.pmean.vec <- NULL
  if(summary.function=='L'){Ldata<- L.pmean.vec <- NULL}
  if(summary.function=='g'){gdata<- g.pmean.vec <- NULL}
  for(i in unique(data[['image_number']])){
    if(sum(data[['image_number']] == i & data[['cell_type']] == from.cell)>qc.cellcount.cutoff & #quality control criteria
       sum(data[['image_number']] == i & data[['cell_type']] == to.cell)>qc.cellcount.cutoff){
      qc.data <- data[data[['image_number']] == i,]
      if(perm.yn==T){
        permuted.K <- NULL
        if(summary.function=='L'){permuted.L <- NULL}
        if(summary.function=='g'){permuted.g <- NULL}
        for(p in 1:P){
          ppp <- spatstat.geom::ppp(x = qc.data$cell_x, y = qc.data$cell_y,
                     marks = factor(sample(qc.data[['cell_type']],
                                           size=length(qc.data[['cell_type']]), replace = F)), window = W)
          Kdata.temp <- data.frame(spatstat.explore::Kcross(ppp, i = from.cell, j = to.cell, correction='Ripley', r=seq(0,R,by=inc) ), image=i)
          permuted.K <- cbind(permuted.K, Kdata.temp$iso)
          if(summary.function=='L'){permuted.L <- sqrt(permuted.K/pi)}
          if(summary.function=='g'){
            gdata.temp <- data.frame(spatstat.explore::pcf(spatstat.explore::Kcross(ppp, i = from.cell, j = to.cell, correction='Ripley', r=seq(0,R,by=inc) ),
                                         method='c', divisor ="d"), image=i)
            permuted.g <- cbind(permuted.g, gdata.temp$pcf)
          }
        }
        K.pmean <- apply(permuted.K, 1, mean);  K.pmean.vec <- c(K.pmean.vec, K.pmean)
        if(summary.function=='L'){
          L.pmean <- apply(permuted.L, 1, mean); L.pmean.vec <- c(L.pmean.vec, L.pmean)
        }
        if(summary.function=='g'){
          g.pmean <- apply(permuted.g, 1, mean); g.pmean.vec <- c(g.pmean.vec, g.pmean)
        }
        print(paste0('permuted outcome for image ',i,' calculated'))
      }
      ppp <- spatstat.geom::ppp(x = qc.data$cell_x, y = qc.data$cell_y,
                 marks = factor(qc.data[['cell_type']]), window = W)
      Kdata.temp <- data.frame(spatstat.explore::Kcross(ppp, i = from.cell, j = to.cell, correction='Ripley', r=seq(0,R,by=inc) ), image=i)
      names(Kdata.temp) <- c('r','K.expect','K.obs','image_number')
      Kdata <- rbind(Kdata, Kdata.temp)
      if(summary.function=='L'){
        Ldata <- Kdata
        Ldata$K.expect <- sqrt(Ldata$K.expect/pi)
        Ldata$K.obs <- sqrt(Ldata$K.obs/pi)
        names(Ldata) <- c('r','L.expect','L.obs','image_number')
      }
      if(summary.function=='g'){
        gdata.temp <- data.frame(spatstat.explore::pcf(spatstat.explore::Kcross(ppp, i = from.cell, j = to.cell, correction='Ripley', r=seq(0,R,by=inc) ),
                                     method='c', divisor ="d"), image=i)
        names(gdata.temp) <- c('r','g.expect','g.obs','image_number')
        gdata <- rbind(gdata, gdata.temp)
      }
    }
  }
  #Calculate the outcome variable
  if(summary.function=='L'){
    sumfun.data <- Ldata
    if(perm.yn==T){sumfun.data$L.pmean <- L.pmean.vec;
    sumfun.data$outcome <- sumfun.data$L.obs - sumfun.data$L.pmean
    }else{
      sumfun.data$L.pmean <- NA;
      sumfun.data$outcome <- sumfun.data$L.obs - sumfun.data$L.expect}
    names(sumfun.data) <- c('r','L.expect','L.obs','image_number','L.pmean','outcome')
  }else if(summary.function=='g'){
    sumfun.data <- gdata
    if(perm.yn==T){sumfun.data$g.pmean <- g.pmean.vec;
    sumfun.data$outcome <- sumfun.data$g.obs - sumfun.data$g.pmean
    }else{
      sumfun.data$g.pmean <- NA;
      sumfun.data$outcome <- sumfun.data$g.obs - sumfun.data$g.expect}
    names(sumfun.data) <- c('r','g.expect','g.obs','image_number','g.pmean','outcome')
  }else{
    sumfun.data <- Kdata
    if(perm.yn==T){sumfun.data$K.pmean <- K.pmean.vec;
    sumfun.data$outcome <- sumfun.data$K.obs - sumfun.data$K.pmean
    }else{
      sumfun.data$K.pmean <- NA;
      sumfun.data$outcome <- sumfun.data$K.obs - sumfun.data$K.expect}
    names(sumfun.data) <- c('r','K.expect','K.obs','image_number','K.pmean','outcome')
  }
  covariate.df <- data
  covariate.df$cell_id <- covariate.df$cell_x <- covariate.df$cell_y <- covariate.df[['cell_type']] <- NULL
  covariate.df <- unique(covariate.df)
  sumfun.data <- merge(sumfun.data, covariate.df, by='image_number'
                       , all=F) #all=F:  if an image doesn't meet the qc.cutoff, don't include its covariates
  return(sumfun.data)
}
