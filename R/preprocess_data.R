# This function preprocesses the data by:
#   1.) filtering out sparse images using qc.cellcount.cutoff
#   2.) calculating the observed spatial summary function to be modeled
#   3.) calculating either the expected or the permuted spatial summary function
#   4.) compiling the outcome variable and covariates into one dataset (model.data)

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
#' @return It returns model.data
#'
#' @export

preprocess_data <- function(data, from.cell, to.cell, qc.cellcount.cutoff=20, P=50, perm.yn=F,
                            R=200, inc=1, image.dims, summary.function='L'){

  image.xmax <- image.dims[2]
  image.ymax <- image.dims[4]
  image.xmin <- image.dims[1]
  image.ymin <- image.dims[3]
  W <- owin(c(image.xmin,image.xmax), c(image.ymin,image.ymax))

  #########################
  ##Calculate model.data ##
  #########################

  if(!(summary.function %in% c('K', 'L', 'g'))){stop("summary.function must be one of 'K', 'L', or g")}

  Kdata <- K.pmean.vec <- NULL
  if(summary.function=='L'){Ldata<- L.pmean.vec <- NULL}
  if(summary.function=='g'){gdata<- g.pmean.vec <- NULL}
  for(i in unique(data$image_number)){
    if(sum(data$image_number == i & data$cell_metacluster == from.cell)>qc.cellcount.cutoff & #quality control criteria
       sum(data$image_number == i & data$cell_metacluster == to.cell)>qc.cellcount.cutoff){
      qc.data <- data[data$image_number == i,]
      if(perm.yn==T){
        permuted.K <- NULL
        if(summary.function=='L'){permuted.L <- NULL}
        if(summary.function=='g'){permuted.g <- NULL}
        for(p in 1:P){
          ppp <- ppp(x = qc.data$cell_x, y = qc.data$cell_y,
                     marks = factor(sample(qc.data$cell_metacluster,
                                           size=length(qc.data$cell_metacluster), replace = F)), window = W)
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
      ppp <- ppp(x = qc.data$cell_x, y = qc.data$cell_y,
                 marks = factor(qc.data$cell_metacluster), window = W)
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
    model.data <- Ldata
    if(perm.yn==T){model.data$L.pmean <- L.pmean.vec;
                   model.data$outcome <- model.data$L.obs - model.data$L.pmean
                  }else{
                   model.data$L.pmean <- NA;
                   model.data$outcome <- model.data$L.obs - model.data$L.expect}
    names(model.data) <- c('r','L.expect','L.obs','image_number','L.pmean','outcome')
  }else if(summary.function=='g'){
    model.data <- gdata
    if(perm.yn==T){model.data$g.pmean <- g.pmean.vec;
                   model.data$outcome <- model.data$g.obs - model.data$g.pmean
                  }else{
                   model.data$g.pmean <- NA;
                   model.data$outcome <- model.data$g.obs - model.data$g.expect}
    names(model.data) <- c('r','g.expect','g.obs','image_number','g.pmean','outcome')
  }else{
    model.data <- Kdata
    if(perm.yn==T){model.data$K.pmean <- K.pmean.vec;
                   model.data$outcome <- model.data$K.obs - model.data$K.pmean
                  }else{
                   model.data$K.pmean <- NA;
                   model.data$outcome <- model.data$K.obs - model.data$K.expect}
    names(model.data) <- c('r','K.expect','K.obs','image_number','K.pmean','outcome')
  }
  covariate.df <- data
  covariate.df$cell_id <- covariate.df$cell_x <- covariate.df$cell_y <- covariate.df$cell_metacluster <- NULL
  covariate.df <- unique(covariate.df)
  model.data <- merge(model.data, covariate.df, by='image_number'
                      , all=F) #all=F:  if an image doesn't meet the qc.cutoff, don't include its covariates
  return(model.data)
}
