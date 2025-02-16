# This function preprocesses the data by:
#   1.) filtering out sparse images using qc.cellcount.cutoff
#   2.) calculating the observed spatial summary function to be modeled
#   3.) calculating either the expected or the permuted spatial summary function
#   4.) compiling the outcome variable and covariates into one dataset (sumfun.data)

#' @title Preprocesses data
#'
#' @param data is
#' @param from.cell is
#' @param to.cell is
#' @param qc.cellcount.cutoff is
#' @param P is
#' @param perm.yn is
#' @param R is
#' @param inc is
#' @param image.dims is
#' @param summary.function is
#' @param seed only makes a difference if perm.yn is switched on
#' @return It returns sumfun.data
#'
#' @export

preprocess_data <- function(data=melanoma_data_subset, from.cell="Macrophage", to.cell="CD8- T cell",
                            qc.cellcount.cutoff=20, P=50, perm.yn=T, cell.label = 'cell_type', image.id='image_number',
                            cell.xcoord = 'cell_x', cell.ycoord = 'cell_y', R=200, inc=1, image.dims=c(0,1200,0,1200), summary.function='g',seed=456){
  set.seed(seed)
  image.xmax <- image.dims[2]
  image.ymax <- image.dims[4]
  image.xmin <- image.dims[1]
  image.ymin <- image.dims[3]
  W <- spatstat.geom::owin(c(image.xmin,image.xmax), c(image.ymin,image.ymax))

  #########################
  ##Calculate sumfun.data ##
  #########################

  if(!(summary.function %in% c('K', 'L', 'g'))){stop("summary.function must be one of 'K', 'L', or g")}

  Kdata <- K.pmean.vec <- NULL
  if(summary.function=='L'){Ldata<- L.pmean.vec <- NULL}
  if(summary.function=='g'){gdata<- g.pmean.vec <- NULL}
  for(i in unique(data[[paste0(image.id)]])){
    if(sum(data[[paste0(image.id)]] == i & data[[paste0(cell.label)]] == from.cell)>qc.cellcount.cutoff & #quality control criteria
       sum(data[[paste0(image.id)]] == i & data[[paste0(cell.label)]] == to.cell)>qc.cellcount.cutoff){
      qc.data <- data[data[[paste0(image.id)]] == i,]
      if(perm.yn==T){
        permuted.K <- NULL
        if(summary.function=='L'){permuted.L <- NULL}
        if(summary.function=='g'){permuted.g <- NULL}
        for(p in 1:P){
          ppp <- ppp(x = qc.data[[paste0(cell.xcoord)]], y = qc.data[[paste0(cell.ycoord)]],
                     marks = factor(sample(qc.data[[paste0(cell.label)]],
                                           size=length(qc.data[[paste0(cell.label)]]), replace = F)), window = W)
          Kdata.temp <- data.frame(Kcross(ppp, i = from.cell, j = to.cell, correction='Ripley', r=seq(0,R,by=inc) ), image=i)
          permuted.K <- cbind(permuted.K, Kdata.temp$iso)
          if(summary.function=='L'){permuted.L <- sqrt(permuted.K/pi)}
          if(summary.function=='g'){
            gdata.temp <- data.frame(pcf(Kcross(ppp, i = from.cell, j = to.cell, correction='Ripley', r=seq(0,R,by=inc) ),
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
      ppp <- ppp(x = qc.data[[paste0(cell.xcoord)]], y = qc.data[[paste0(cell.ycoord)]],
                 marks = factor(qc.data[[paste0(cell.label)]]), window = W)
      Kdata.temp <- data.frame(Kcross(ppp, i = from.cell, j = to.cell, correction='Ripley', r=seq(0,R,by=inc) ), image=i)
      names(Kdata.temp) <- c('r','K.expect','K.obs','image_number')
      Kdata <- rbind(Kdata, Kdata.temp)
      if(summary.function=='L'){
        Ldata <- Kdata
        Ldata$K.expect <- sqrt(Ldata$K.expect/pi)
        Ldata$K.obs <- sqrt(Ldata$K.obs/pi)
        names(Ldata) <- c('r','L.expect','L.obs','image_number')
      }
      if(summary.function=='g'){
        gdata.temp <- data.frame(pcf(Kcross(ppp, i = from.cell, j = to.cell, correction='Ripley', r=seq(0,R,by=inc) ),
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
  covariate.df$cell_id <- covariate.df[[paste0(cell.xcoord)]] <- covariate.df[[paste0(cell.ycoord)]] <- covariate.df[[paste0(cell.label)]] <- NULL
  covariate.df <- unique(covariate.df)
  sumfun.data <- merge(sumfun.data, covariate.df, by='image_number'
                      , all=F) #all=F:  if an image doesn't meet the qc.cutoff, don't include its covariates
  return(sumfun.data)
}
