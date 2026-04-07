#' @title Preprocesses single-cell spatial data for colocalization analysis
#'
#' @description Discards sparse images, calculates the selected spatial summary function, and calculates the colocalization curve to be used as the outcome in downstream functional regression.
#'
#' @param data Data frame containing the variables `patient_id`, `image_number`, `cell_id`, `cell_x`, `cell_y`, `cell_type`, and any covariates to be adjusted for, as in the example data.
#' @param from_cell The "from" cell type to be used in the calculation of the spatial summary function
#' @param to_cell The "to" cell type to be used in the calculation of the spatial summary function
#' @param qc_cellcount_cutoff If an image has `qc_cellcount_cutoff` or less "to" cells, or `qc_cellcount_cutoff` or less "from" cells, the image is excluded from the calculation of the spatial summary function.  (Similar to SpaceANOVA's `Hard_ths` parameter.)
#' @param n_perm The number of permutations to use in the permutation adjustment to the colocalization curve, if `perm_yn = TRUE`.
#' @param perm_yn If `TRUE`, the permuted mean of the spatial summary function is calculated using the procedure described in the [Wilson paper](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009900).  The colocalization curve is then calculated by subtracting the value of the permuted mean from the value of the observed spatial summary function.  If `perm_yn = FALSE`, the curve is calculated by subtracting the summary function's expected value assuming complete spatial randomness from the observed spatial summary function.
#' @param r_max The maximum radius at which the observed spatial summary function should be evaluated.
#' @param inc The increment of units by which the spatial summary function should be evaluated.  For example, if `r_max = 200` and `inc=1`, the spatial summary function is evaluated over the grid `0:200`.
#' @param image_dims Vector containing the dimensions of the images in the observed data. For example, `c(x_min, x_max, y_min, y_max)`.  All images must be of the same dimension.
#' @param summary_function The spatial summary function to be calculated.  Must be one of `K`, `L`, or `g`.
#' @param verbose If `TRUE`, prints progress
#'
#' @return A data frame containing values of the observed and expected spatial summary function for each image at each radius, as well as values of the covariates carried forward from the input data. Depending on `summary_function`, the output includes `K_*`, `L_*`, or `g_*` columns. When `perm_yn = TRUE`, the output also includes the corresponding permuted mean column (`K_pmean`, `L_pmean`, or `g_pmean`). The `outcome` column is defined as the observed summary function minus the expected summary function when `perm_yn = FALSE`, or the observed summary function minus the permuted mean when `perm_yn = TRUE`.
#'
#' @export

preprocess_data <- function(data,
                            from_cell,
                            to_cell,
                            qc_cellcount_cutoff=0,
                            n_perm=50,
                            perm_yn=FALSE,
                            r_max=200,
                            inc=1,
                            image_dims,
                            summary_function='L',
                            verbose=TRUE){
  image_xmax <- image_dims[2]
  image_ymax <- image_dims[4]
  image_xmin <- image_dims[1]
  image_ymin <- image_dims[3]
  W <- spatstat.geom::owin(c(image_xmin,image_xmax),
                           c(image_ymin,image_ymax))

  #########################
  ##Calculate sumfun_data #
  #########################

  if(!(summary_function %in% c('K', 'L', 'g'))){
    stop("summary_function must be one of 'K', 'L', or 'g' ")
    }

  Kdata <- K_pmean_vec <- NULL
  if(summary_function=='L'){
    Ldata<- L_pmean_vec <- NULL
    }
  if(summary_function=='g'){
    gdata<- g_pmean_vec <- NULL
    }
  for(i in unique(data$image_number)){
    if(sum(data$image_number == i & data$cell_type == from_cell)>qc_cellcount_cutoff & #quality control criteria
       sum(data$image_number == i & data$cell_type == to_cell)>qc_cellcount_cutoff){
      qc_data <- data[data$image_number == i,]
      if(perm_yn==TRUE){
        permuted_K <- NULL
        if(summary_function=='L'){
          permuted_L <- NULL
          }
        if(summary_function=='g'){
          permuted_g <- NULL
          }
        for(p in 1:n_perm){
          ppp <- spatstat.geom::ppp(x = qc_data$cell_x,
                                    y = qc_data$cell_y,
                                    marks = factor(sample(qc_data$cell_type,
                                                          size=length(qc_data$cell_type),
                                                          replace = FALSE)
                                                   ),
                                    window = W)
          Kdata_temp <- data.frame(spatstat.explore::Kcross(ppp, i = from_cell, j = to_cell, correction='Ripley', r=seq(0,r_max,by=inc) ), image=i)
          permuted_K <- cbind(permuted_K, Kdata_temp$iso)
          if(summary_function=='L'){permuted_L <- sqrt(permuted_K/pi)}
          if(summary_function=='g'){
            gdata_temp <- data.frame(spatstat.explore::pcf(spatstat.explore::Kcross(ppp, i = from_cell, j = to_cell, correction='Ripley', r=seq(0,r_max,by=inc) ),
                                         method='c', divisor ="d"), image=i)
            permuted_g <- cbind(permuted_g, gdata_temp$pcf)
          }
        }
        K_pmean <- apply(permuted_K, 1, mean);  K_pmean_vec <- c(K_pmean_vec, K_pmean)
        if(summary_function=='L'){
          L_pmean <- apply(permuted_L, 1, mean); L_pmean_vec <- c(L_pmean_vec, L_pmean)
        }
        if(summary_function=='g'){
          g_pmean <- apply(permuted_g, 1, mean); g_pmean_vec <- c(g_pmean_vec, g_pmean)
        }
        if (verbose) message(paste0('permuted outcome for image ',i,' calculated'))
      }
      ppp <- spatstat.geom::ppp(x = qc_data$cell_x, y = qc_data$cell_y,
                 marks = factor(qc_data$cell_type), window = W)
      Kdata_temp <- data.frame(spatstat.explore::Kcross(ppp, i = from_cell, j = to_cell, correction='Ripley', r=seq(0,r_max,by=inc) ), image=i)
      names(Kdata_temp) <- c('r','K_expect','K_obs','image_number')
      Kdata <- rbind(Kdata, Kdata_temp)
      if(summary_function=='L'){
        Ldata <- Kdata
        Ldata$K_expect <- sqrt(Ldata$K_expect/pi)
        Ldata$K_obs <- sqrt(Ldata$K_obs/pi)
        names(Ldata) <- c('r','L_expect','L_obs','image_number')
      }
      if(summary_function=='g'){
        gdata_temp <- data.frame(spatstat.explore::pcf(spatstat.explore::Kcross(ppp, i = from_cell, j = to_cell, correction='Ripley', r=seq(0,r_max,by=inc) ),
                                     method='c', divisor ="d"), image=i)
        names(gdata_temp) <- c('r','g_expect','g_obs','image_number')
        gdata <- rbind(gdata, gdata_temp)
      }
    }
  }
  #Calculate the outcome variable
  if(summary_function=='L'){
    sumfun_data <- Ldata
    if(perm_yn==TRUE){sumfun_data$L_pmean <- L_pmean_vec;
    sumfun_data$outcome <- sumfun_data$L_obs - sumfun_data$L_pmean
    }else{
      sumfun_data$L_pmean <- NA;
      sumfun_data$outcome <- sumfun_data$L_obs - sumfun_data$L_expect}
    names(sumfun_data) <- c('r','L_expect','L_obs','image_number','L_pmean','outcome')
  }else if(summary_function=='g'){
    sumfun_data <- gdata
    if(perm_yn==TRUE){sumfun_data$g_pmean <- g_pmean_vec;
    sumfun_data$outcome <- sumfun_data$g_obs - sumfun_data$g_pmean
    }else{
      sumfun_data$g_pmean <- NA;
      sumfun_data$outcome <- sumfun_data$g_obs - sumfun_data$g_expect}
    names(sumfun_data) <- c('r','g_expect','g_obs','image_number','g_pmean','outcome')
  }else{
    sumfun_data <- Kdata
    if(perm_yn==TRUE){sumfun_data$K_pmean <- K_pmean_vec;
    sumfun_data$outcome <- sumfun_data$K_obs - sumfun_data$K_pmean
    }else{
      sumfun_data$K_pmean <- NA;
      sumfun_data$outcome <- sumfun_data$K_obs - sumfun_data$K_expect}
    names(sumfun_data) <- c('r','K_expect','K_obs','image_number','K_pmean','outcome')
  }
  covariate_df <- data
  covariate_df$cell_id <- covariate_df$cell_x <- covariate_df$cell_y <- covariate_df$cell_type <- NULL
  covariate_df <- unique(covariate_df)
  sumfun_data <- merge(sumfun_data, covariate_df, by='image_number'
                       , all=FALSE) #all=FALSE:  if an image doesn't meet the qc_cutoff, don't include it
  return(sumfun_data)
}
