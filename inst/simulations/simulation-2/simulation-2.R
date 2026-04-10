library(spatstat.data)
library(spatstat.univar)
library(spatstat.geom)
library(spatstat.random)
library(spatstat.explore)
library(stringr)
library(refund)
library(parallel)
library(ggplot2)
library(patchwork)

#------------ Setup ----------------------------------------------------------------------------------------------
nCores <- 24
nthreads <- 1
B1   <- 50
B2   <- 50
alpha <- 0.05
seed <- 12345

nPatients <- 100 #number of patients
nIm <- 1         #number of images per patient
counts <- seq(400, 700, by = 10)
sigmab  <- 10
deltab <- 7

#------------ Define utility functions ---------------------------------------------------------------------------
format_spatial_variable <- function(data, spatial_variable, grid, id){
  column <- function(i){data[data[id] == i,][spatial_variable][1] } #calling the outcome function Y
  mat_transpose <- data.frame(grid = grid)
  for(i in sort(unique(data[id])[[1]])){
    mat_transpose <- cbind(mat_transpose, column(i))
  }
  mat <- t(mat_transpose[,-c(1)])
  return(mat)
}

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
    L_pmean_vec <- NULL
  }
  if(summary_function=='g'){
    gdata<- g_pmean_vec <- NULL
  }
  for(i in unique(data$image_number)){
    from_count <- sum(data$image_number == i & data$cell_type == from_cell)
    to_count <- sum(data$image_number == i & data$cell_type == to_cell)

    if(from_count > qc_cellcount_cutoff &
       to_count > qc_cellcount_cutoff){ #quality control criteria
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
      if(summary_function=='g'){
        gdata_temp <- data.frame(spatstat.explore::pcf(spatstat.explore::Kcross(ppp, i = from_cell, j = to_cell, correction='Ripley', r=seq(0,r_max,by=inc) ),
                                                       method='c', divisor ="d"), image=i)
        names(gdata_temp) <- c('r','g_expect','g_obs','image_number')
        gdata <- rbind(gdata, gdata_temp)
      }
    }else if(from_count > qc_cellcount_cutoff &
             to_count == 0){
      r_grid <- seq(0, r_max, by = inc)

      Kdata_temp <- data.frame(
        r = r_grid,
        K_expect = NA_real_,
        K_obs = 0,
        image_number = i
      )
      Kdata <- rbind(Kdata, Kdata_temp)

      if(summary_function=='g'){
        gdata_temp <- data.frame(
          r = r_grid,
          g_expect = NA_real_,
          g_obs = 0,
          image_number = i
        )
        gdata <- rbind(gdata, gdata_temp)
      }

      if(perm_yn==TRUE){
        K_pmean_vec <- c(K_pmean_vec, rep(NA_real_, nrow(Kdata_temp)))
        if(summary_function=='L'){
          L_pmean_vec <- c(L_pmean_vec, rep(NA_real_, nrow(Kdata_temp)))
        }
        if(summary_function=='g'){
          g_pmean_vec <- c(g_pmean_vec, rep(NA_real_, nrow(Kdata_temp)))
        }
      }
    }
  }
  #Calculate the outcome variable
  if(summary_function=='L'){
    sumfun_data <- Kdata
    sumfun_data$K_expect <- sqrt(sumfun_data$K_expect/pi)
    sumfun_data$K_obs <- sqrt(sumfun_data$K_obs/pi)
    names(sumfun_data) <- c('r','L_expect','L_obs','image_number')
    if(perm_yn==TRUE){
      sumfun_data$L_pmean <- L_pmean_vec
      sumfun_data$outcome <- ifelse(is.na(sumfun_data$L_pmean) & sumfun_data$L_obs == 0,
                                    0,
                                    sumfun_data$L_obs - sumfun_data$L_pmean)
    }else{
      sumfun_data$L_pmean <- NA
      sumfun_data$outcome <- ifelse(is.na(sumfun_data$L_expect) & sumfun_data$L_obs == 0,
                                    0,
                                    sumfun_data$L_obs - sumfun_data$L_expect)
    }
    names(sumfun_data) <- c('r','L_expect','L_obs','image_number','L_pmean','outcome')
  }else if(summary_function=='g'){
    sumfun_data <- gdata
    if(perm_yn==TRUE){
      sumfun_data$g_pmean <- g_pmean_vec
      sumfun_data$outcome <- ifelse(is.na(sumfun_data$g_pmean) & sumfun_data$g_obs == 0,
                                    0,
                                    sumfun_data$g_obs - sumfun_data$g_pmean)
    }else{
      sumfun_data$g_pmean <- NA
      sumfun_data$outcome <- ifelse(is.na(sumfun_data$g_expect) & sumfun_data$g_obs == 0,
                                    0,
                                    sumfun_data$g_obs - sumfun_data$g_expect)
    }
    names(sumfun_data) <- c('r','g_expect','g_obs','image_number','g_pmean','outcome')
  }else{
    sumfun_data <- Kdata
    if(perm_yn==TRUE){
      sumfun_data$K_pmean <- K_pmean_vec
      sumfun_data$outcome <- ifelse(is.na(sumfun_data$K_pmean) & sumfun_data$K_obs == 0,
                                    0,
                                    sumfun_data$K_obs - sumfun_data$K_pmean)
    }else{
      sumfun_data$K_pmean <- NA
      sumfun_data$outcome <- ifelse(is.na(sumfun_data$K_expect) & sumfun_data$K_obs == 0,
                                    0,
                                    sumfun_data$K_obs - sumfun_data$K_expect)
    }
    names(sumfun_data) <- c('r','K_expect','K_obs','image_number','K_pmean','outcome')
  }
  covariate_df <- data
  covariate_df$cell_id <- covariate_df$cell_x <- covariate_df$cell_y <- covariate_df$cell_type <- NULL
  covariate_df <- unique(covariate_df)
  sumfun_data <- merge(sumfun_data, covariate_df, by='image_number'
                       , all=FALSE) #all=FALSE:  if an image doesn't meet the qc_cutoff, don't include it
  return(sumfun_data)
}

extract_target_curve <- function(model,
                                 grid,
                                 target,
                                 form,
                                 covar_list = list(),
                                 func) {

  #Possible errors
  if (length(target) == 0) stop("No symbols found in target.")
  if (length(form) != length(target)) {
    stop("form must have length equal to the number of variables in target (",
         length(target), "). The length given was ", length(form), ".")
  }
  if (any(!form %in% c("coef", "term"))) {
    stop('Each entry of form must be either "coef" or "term".')
  }
  if (!is.function(func)) stop("func must be a function.")
  if (any(names(formals(func)) != target)) stop("Arguments of func must match target.")
  if (any(names(covar_list) != target)) stop("Names of covar_list must match target.")

  #Extract components
  sm <- stats::coef(model, n1 = length(grid))$smterms
  get_beta <- function(name) {
    name <- paste0(name, "(grid)")
    if (name %in% names(sm)) result <- sm[[name]]$coef$value
    if (length(result) > length(grid)) stop("Multiple smooth terms match for ", name)
    if (length(result) == 0)stop("No smooth term found matching ", name)
    if (name == 'Intercept(grid)'){
      result <- result + rep(summary(model)$p.coeff, length(grid)) #include the constant component of beta0
    }
    return(result)
  }

  #Prepare the list to perform func on
  components <- vector("list", length(target))
  names(components) <- target
  for (i in seq_along(target)) {
    v <- target[[i]]
    beta_v <- get_beta(v)
    components[[i]] <- if (form[[i]] == "coef") beta_v else beta_v * covar_list[[i]]
  }
  return(unname(do.call(func, components)))
}

get_fixed_and_re <- function(model, n_Im_total, grid, id) {
  Tm <- mgcv::predict.gam(model, type = "terms")
  nm <- colnames(Tm)
  re_col <- grep(id, nm)
  if (length(re_col) == 0) re_col <- NULL #if no RE, delete re_col and return NA for b0
  fx_cols <- grep("^(?!ti\\()", nm, perl = TRUE) #Get all columns that are NOT random effects
  if (length(fx_cols) == 0) stop("Could not find fixed effects.")

  b0 <- matrix(Tm[, re_col], nrow = n_Im_total, ncol = length(grid), byrow = TRUE)
  fixed <- matrix(rowSums(Tm[, fx_cols, drop = FALSE]) + as.numeric(summary(model)$p.coeff),
                  nrow = n_Im_total, ncol = length(grid), byrow = TRUE)
  list(fixed = fixed, b0 = b0)
}


wildBS_CB <- function(formula,
                      data,
                      spatial_covars = NULL,
                      target,
                      form,
                      covar_list=list(),
                      func,
                      B1,
                      B2,
                      alpha=.05,
                      re=NULL,
                      id,
                      nthreads = 1){

  ####################
  ## Format dataset ##
  ####################

  n_patients <- dim(unique(data['patient_id']))[1]
  n_Im_total <- dim(unique(data['image_number']))[1]
  grid <- sort(unique(data$r))

  #format image-level covariates
  image_level_covars <- unique(
    data.frame(
      data[["image_number"]],
      data[all.vars(formula)[c(-1, -which(all.vars(formula) %in% spatial_covars))]]
    )
  )
  names(image_level_covars)[1] <- 'image_number'
  pffr_data <- image_level_covars
  if(dim(pffr_data)[1]!=n_Im_total){stop('spatial_covars argument may be invalid')}

  #format outcome function
  outcome <- format_spatial_variable(
    data = data,
    spatial_variable = "outcome",
    grid = grid,
    id = "image_number"
  )
  pffr_data$outcome <- outcome

  #format spatial covariates
  if (!is.null(spatial_covars)) {
    for (i in seq_along(spatial_covars)) {
      pffr_data$temp <- format_spatial_variable(
        data = data,
        spatial_variable = spatial_covars[i],
        grid = grid,
        id = "image_number"
      )
      names(pffr_data)[names(pffr_data) == "temp"] <- paste0(spatial_covars[i])
    }
  }

  pffr_data <- as.data.frame(pffr_data)

  ###############
  ## Fit model ##
  ###############

  pffrmodel <- refund::pffr(
    formula   = formula,
    data      = pffr_data,
    yind      = grid,
    algorithm = "bam",
    discrete  = TRUE,
    nthreads  = nthreads
  )

  #Get estimates and extract target curve
  pred_mat <- stats::predict(pffrmodel)
  epsilon <- outcome - pred_mat
  fixed_and_b0 <- get_fixed_and_re(model=pffrmodel, n_Im_total=n_Im_total, grid=grid, id=id)
  fixed <- fixed_and_b0$fixed
  if(!is.null(re)){b0 <- fixed_and_b0$b0}
  if(is.null(re)){b0 <- matrix(0, nrow=n_Im_total, ncol=length(grid))}
  target_curve <- extract_target_curve(model=pffrmodel,
                                       grid=grid,
                                       target=target,
                                       form=form,
                                       covar_list = covar_list,
                                       func=func)


  ###########################
  ## Outer bootstrap layer ##
  ###########################

  #Image-level multipliers
  c_img <- 2 * stats::rbinom(n_Im_total * B1, size = 1, prob = .5) - 1 #1 or -1
  c_img <- matrix(c_img, nrow = n_Im_total, ncol = B1)
  #Patient-level multipliers
  if(!is.null(re)){
    c_patient <- 2 * stats::rbinom(n_patients * B1, size = 1, prob = .5) - 1 #1 or -1
    c_patient <- matrix(c_patient, nrow = n_patients, ncol = B1)
    c_patient_df <- merge(pffr_data[, c("image_number", "patient_id")],
                          data.frame(patient_id = unique(pffr_data$patient_id),
                                     c_patient),
                          by = "patient_id",
                          all.x = TRUE)
    ord <- match(pffr_data$image_number, c_patient_df$image_number)
    c_patient <- as.matrix(c_patient_df[ord, -(1:2), drop = FALSE])
  }else{
    c_patient <- matrix(0, nrow = n_Im_total, ncol = B1)
  }
  #Storage
  M <- numeric(B1)
  target_curve_outer <- matrix(NA, nrow = B1, ncol = length(grid))

  #Do the outer layer bootstrap
  for (b1 in 1:B1) {
    #resample data
    Y_mat_bs <- fixed + c_patient[, b1] * b0 + c_img[, b1] * epsilon
    #fit model and get estimates
    pffr_data_bs <- pffr_data
    pffr_data_bs$Y_mat_bs <- Y_mat_bs
    formula_bs <- stats::update.formula(formula, stats::as.formula('Y_mat_bs ~ .'))
    pffrmodel_bs <- refund::pffr(
      formula_bs,
      data = pffr_data_bs,
      yind = grid,
      algorithm = "bam",
      discrete = TRUE,
      nthreads = nthreads
    )
    pred_mat_bs <- stats::predict(pffrmodel_bs)
    epsilon_bs <- Y_mat_bs - pred_mat_bs
    fixed_and_b0_bs <- get_fixed_and_re(model=pffrmodel_bs, n_Im_total=n_Im_total, grid=grid, id=id)
    fixed_bs <- fixed_and_b0_bs$fixed
    if(!is.null(re)){b0_bs <- fixed_and_b0_bs$b0}
    if(is.null(re)){b0_bs <- matrix(0, nrow=n_Im_total, ncol=length(grid))}
    target_curve_bs <- extract_target_curve(model=pffrmodel_bs,
                                            grid=grid,
                                            target=target,
                                            form=form,
                                            covar_list = covar_list,
                                            func=func)

    ###########################
    ## Inner bootstrap layer ##
    ###########################

    #Image-level multipliers
    c_img_bs <- 2 * stats::rbinom(n_Im_total * B2, size = 1, prob = .5) - 1 #1 or -1
    c_img_bs <- matrix(c_img_bs, nrow = n_Im_total, ncol = B2)
    #Patient-level multipliers
    if(!is.null(re)){
      c_patient_bs <- 2 * stats::rbinom(n_patients * B2, size = 1, prob = .5) - 1 #1 or -1
      c_patient_bs <- matrix(c_patient_bs, nrow = n_patients, ncol = B2)
      c_patient_df_bs <- merge(pffr_data_bs[, c("image_number", "patient_id")],
                               data.frame(patient_id = unique(pffr_data_bs$patient_id),
                                          c_patient_bs),
                               by = "patient_id",
                               all.x = TRUE)
      ord_bs <- match(pffr_data_bs$image_number, c_patient_df_bs$image_number)
      c_patient_bs <- as.matrix(c_patient_df_bs[ord_bs, -(1:2), drop = FALSE])
    }else{
      c_patient_bs <- matrix(0, nrow = n_Im_total, ncol = B2)
    }

    M_inner <- matrix(NA, nrow = B2, ncol = length(grid))
    for (b2 in 1:B2) {
      #Resample data
      Y_mat_bs_bs <- fixed_bs + c_patient_bs[, b2] * b0_bs + c_img_bs[, b2] * epsilon_bs
      #fit model and get estimates
      pffr_data_bs_bs <- pffr_data_bs
      pffr_data_bs_bs$Y_mat_bs_bs <- Y_mat_bs_bs
      formula_bs_bs <- stats::update.formula(formula_bs, stats::as.formula('Y_mat_bs_bs ~ .'))
      pffrmodel_bs_bs <- refund::pffr(
        formula_bs_bs,
        data = pffr_data_bs_bs,
        yind = grid,
        algorithm = "bam",
        discrete = TRUE,
        nthreads = nthreads
      )
      target_curve_bs_bs <- extract_target_curve(model=pffrmodel_bs_bs,
                                                 grid=grid,
                                                 target=target,
                                                 form=form,
                                                 covar_list = covar_list,
                                                 func=func)
      M_inner[b2,] <- target_curve_bs_bs
    }
    SEs_inner <- apply(M_inner, 2, stats::sd)
    SEs_inner <- pmax(SEs_inner, 1e-8) #avoid potential division by 0
    M[b1] <- max(abs((target_curve_bs-target_curve)/SEs_inner))
    target_curve_outer[b1,] <- target_curve_bs
  }

  #########################
  ## Get confidence band ##
  #########################

  q <- stats::quantile(M, 1-alpha)

  target_curve_SEs <- apply(target_curve_outer, 2, stats::sd)
  MoE <- q * target_curve_SEs
  CB_lower <- target_curve - MoE
  CB_upper <- target_curve + MoE
  CBs <- list(CB_lower, CB_upper, grid, target_curve, target_curve_SEs, q, M)
  names(CBs) <- c('CB_lower','CB_upper', 'grid', 'target_curve_estimates', 'target_curve_SEs', 'q', 'M')
  return(CBs)
}

plot_wildBS_CB <- function(CB_object){
  p <- ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::geom_line(ggplot2::aes(x = CB_object$grid,
                                    y = CB_object$target_curve_estimates,
                                    col = "Estimated Target Curve")) +
    ggplot2::geom_line(ggplot2::aes(x = CB_object$grid,
                                    y = CB_object$CB_lower,
                                    col = "Confidence Band")) +
    ggplot2::geom_line(ggplot2::aes(x = CB_object$grid,
                                    y = CB_object$CB_upper,
                                    col = "Confidence Band")) +
    ggplot2::labs(x = "Radius",
                  y = "Target Curve",
                  col = " ") +
    ggplot2::theme(axis.title.y = ggplot2::element_text(size = 20),
                   axis.title.x = ggplot2::element_text(size = 20)) +
    ggplot2::geom_hline(yintercept = 0, lty = 2) +
    ggplot2::scale_color_manual(values = c("Estimated Target Curve" = "#F8766D",
                                           "Confidence Band" = "#00BA38"))
  return(p)
}

#------------ Simulate 1 draw of data -------------------------------------------------------------------------
set.seed(seed)

tumor_prob <- .9 #probability that a given image contains a tumor
tumor_heights <- seq(from = 500, to = 700, by = 100)
tumor_widths <- seq(from = 500, to = 700, by = 100)
sigmat = 10
window <- spatstat.geom::owin(xrange = c(0, 1000),
                              yrange = c(0, 1000))

g1 <- rpois(nPatients/2, sigmab)
g2 <- rpois(nPatients/2, sigmab + deltab )
adjustSigma = c(g1,g2)+1

x <- c()
y <- c()
cellType <- c()
imageID <- c()

for (p in 1:nPatients) {
  for (j in 1:nIm) {
    if(rbinom(1,size=1,prob=tumor_prob)==1){ #if the image contains a tumor
      #Determine tumor region
      tumor_height <- sample(tumor_heights, 1)
      tumor_width <- sample(tumor_widths, 1)
      tumor_center_x <- runif(1, min=window$xrange[1], max=window$xrange[2])
      tumor_center_y <- runif(1, min=window$yrange[1], max=window$yrange[2])
      tumor_xrange <- pmin(pmax(c(tumor_center_x-tumor_width/2, tumor_center_x+tumor_width/2),0),1000)
      tumor_yrange <- pmin(pmax(c(tumor_center_y-tumor_height/2, tumor_center_y+tumor_height/2),0),1000)
      tumor_counts <- seq(from = 100, to = 400, by = 10)*(tumor_height*tumor_width/1000^2)
      tumor_window <- spatstat.geom::owin(xrange = tumor_xrange,
                                          yrange = tumor_yrange)
      #Generate Type A cells in stroma using homogeneous Poisson point process
      sCount1 <- sample(counts,1)
      a <- rpoispp(sCount1/1000^2, win = window)
      #Generate Type A cells in tumor using homogeneous Poisson point process
      a_tumor <- rpoispp(sCount1/1000^2, win = tumor_window)
      aDens_tumor <- density(a_tumor, sigma = adjustSigma[p], kernel = "disc")
      sCount2 <- sCount1
      sCount_tumor <- sCount1
      aDens_tumor$v <- pmax(aDens_tumor$v,0)*sCount2/sCount1
      #Generate B in tumor as colocalized with A
      if(a_tumor$n > 0){b_tumor <- rpoispp(aDens_tumor)}else{b_tumor <- a_tumor}
      #Generate B in stroma using homogeneous Poisson point process
      b <- rpoispp(sCount2/1000^2, win = window)  #B is random in stroma
      #Generate T in tumor as colocalized with A
      aDens_tumor <- density(a_tumor, sigma = sigmat, kernel = "disc")
      aDens_tumor$v <- pmax(aDens_tumor$v,0)*sCount_tumor/sCount1
      if(a_tumor$n > 0){t <- rpoispp(aDens_tumor)} else {t <- a_tumor}

      #Patch tumor and stroma regions together
      num_a_stroma <- length(a$x[a$x < tumor_xrange[1] | a$x > tumor_xrange[2] |
                                   a$y < tumor_yrange[1] | a$y > tumor_yrange[2] ])
      num_b_stroma <- length(b$x[b$x < tumor_xrange[1] | b$x > tumor_xrange[2] |
                                   b$y < tumor_yrange[1] | b$y > tumor_yrange[2] ])

      stroma_x <- c(a$x[a$x < tumor_xrange[1] | a$x > tumor_xrange[2] |
                          a$y < tumor_yrange[1] | a$y > tumor_yrange[2] ],
                    b$x[b$x < tumor_xrange[1] | b$x > tumor_xrange[2] |
                          b$y < tumor_yrange[1] | b$y > tumor_yrange[2] ] )
      tumor_x <- c(a_tumor$x, b_tumor$x, t$x)
      stroma_y <- c(a$y[a$x < tumor_xrange[1] | a$x > tumor_xrange[2] |
                          a$y < tumor_yrange[1] | a$y > tumor_yrange[2] ],
                    b$y[b$x < tumor_xrange[1] | b$x > tumor_xrange[2] |
                          b$y < tumor_yrange[1] | b$y > tumor_yrange[2] ] )
      tumor_y <- c(a_tumor$y, b_tumor$y, t$y)

      x <- c(x, stroma_x, tumor_x)
      y <- c(y, stroma_y, tumor_y)

      cellType <- c(cellType, rep("A", num_a_stroma), rep("B", num_b_stroma),
                    rep("A", a_tumor$n), rep("B", b_tumor$n), rep("T", t$n))
      imageID <- c(imageID, rep(paste(p,j,sep = "_"), num_a_stroma+num_b_stroma+a_tumor$n+b_tumor$n+t$n))
    } else {
      #Generate A and B in entire image region using a homogeneous Poisson point process
      sCount1 <- sample(counts,1)
      sCount2 <- sample(counts,1)
      a <- rpoispp(sCount1/1000^2, win = window)
      b <- rpoispp(sCount2/1000^2, win = window)

      x <- c(x, a$x, b$x)
      y <- c(y, a$y, b$y)

      cellType <- c(cellType, rep("A", a$n), rep("B", b$n))
      imageID <- c(imageID, rep(paste(p,j,sep = "_"), a$n+b$n))
    }
  }
}

imageID <- factor(imageID)

cellExp <- data.frame(
  x = x,
  y = y,
  cellType = factor(cellType),
  imageID = imageID
)


phenoData <- data.frame(imageID = unique(imageID),
                        condition = rep(c("Group1", "Group2"), each = nIm*nPatients/2),
                        subject = rep(1:nPatients, each = nIm))

table <- cellExp
names(table) <- c('Xcoord','Ycoord','cell_type','image')
table$patient_id <- sapply(table$image, function(x){stringr::str_split(x, stringr::fixed("_"))[[1]][1]})
table$roi <- sapply(table$image, function(x){stringr::str_split(x, stringr::fixed("_"))[[1]][2]})
names(phenoData) <- c('image','group','patient_id')
table <- merge(table, phenoData, by=c('image','patient_id'), all=T)
colnames(table) <- c('image_number','patient_id','cell_x', 'cell_y', "cell_type",  'roi', 'group')

#Calculate spatial summary functions
sumfun_data <- preprocess_data(table,
                               from_cell='A',
                               to_cell='B',
                               qc_cellcount_cutoff=0,
                               perm_yn=F,
                               r_max=200,
                               inc=1,
                               image_dims=c(0,1000,0,1000),
                               summary_function='L')
sumfun_data$group <- ifelse(sumfun_data$group=='Group1', 1, 0)
#var_mat represents the values of the spatial covariate \hat{L}_{AT}(r)
var_mat_data <- preprocess_data(table,
                                from_cell='A',
                                to_cell='T',
                                qc_cellcount_cutoff=0,
                                perm_yn=F,
                                r_max=200,
                                inc=1,
                                image_dims=c(0,1000,0,1000),
                                summary_function='L')
var_mat_data$group <- ifelse(var_mat_data$group=='Group1', 1, 0)
names(var_mat_data)[names(var_mat_data)=='L_obs'] <- 'var_mat'
var_mat_data$L_expect <- var_mat_data$L_pmean <- var_mat_data$outcome <- NULL
sumfun_data_full <- merge(sumfun_data, var_mat_data, by=c('image_number','r','patient_id','roi','group'),all=F)
sumfun_data_full <- sumfun_data_full[
  order(sumfun_data_full$image_number,
        sumfun_data_full$r),
]

sumfun_data_full$var_group_mat <- sumfun_data_full$var_mat * sumfun_data_full$group

#For the prediction, use \hat{L}_{AT}(r) from patients 10 and 51 to represent groups 1 and 2, respectively.
sumfun_data_trtpatient <- sumfun_data_full[sumfun_data_full$image_number=='5_1',]
sumfun_data_ctrlpatient <- sumfun_data_full[sumfun_data_full$image_number=='59_1',]

#-------------- Get confidence band for E[Y|X] for a trt patient (overall image) ---------------------------------
covar_list_trt_overall <- list(Intercept = rep(1, length(sumfun_data_trtpatient$group)),
                                  group = sumfun_data_trtpatient$group,
                                  var_mat = sumfun_data_trtpatient$var_mat,
                                  var_group_mat = sumfun_data_trtpatient$var_group_mat)

CB1 <- wildBS_CB(formula=outcome ~ group + var_mat + var_group_mat,
                data=sumfun_data_full,
                spatial_covars = c('var_mat','var_group_mat'),
                target = c("Intercept","group", "var_mat", "var_group_mat"),
                form   = c("coef","term","term","term"),
                covar_list = covar_list_trt_overall,
                func = function(Intercept, group, var_mat, var_group_mat){
                  Intercept + group + var_mat + var_group_mat
                },
                B1=B1,
                B2=B2,
                alpha=.05,
                re=NULL,
                id='patient_id',
                nthreads = 1)
save(CB1, file='CB1.RData')
plot_wildBS_CB(CB1)

#-------------- Get confidence band for E[Y|X] for a ctrl patient (overall image) ---------------------------------

covar_list_ctrl_overall <- list(Intercept = rep(1, length(sumfun_data_ctrlpatient$group)),
                                  group = sumfun_data_ctrlpatient$group,
                                  var_mat = sumfun_data_ctrlpatient$var_mat,
                                  var_group_mat = sumfun_data_ctrlpatient$var_group_mat)

CB2 <- wildBS_CB(formula=outcome ~ group + var_mat + var_group_mat,
                data=sumfun_data_full,
                spatial_covars = c('var_mat','var_group_mat'),
                target = c("Intercept","group", "var_mat", "var_group_mat"),
                form   = c("coef","term","term","term"),
                covar_list = covar_list_ctrl_overall,
                func = function(Intercept, group, var_mat, var_group_mat){
                  Intercept + group + var_mat + var_group_mat
                },
                B1=B1,
                B2=B2,
                alpha=.05,
                re=NULL,
                id='patient_id',
                nthreads = 1)
save(CB2, file='CB2.RData')
plot_wildBS_CB(CB2)

#-------------- Get confidence band for E[Y|X] for a group 1 patient (stroma region only) ---------------------------------
covar_list_trt_stroma <- list(Intercept = rep(1, length(sumfun_data_trtpatient$group)),
                                  group = sumfun_data_trtpatient$group,
                                  var_mat = rep(0, length(sumfun_data_trtpatient$group)),
                                  var_group_mat = rep(0, length(sumfun_data_trtpatient$group)))

CB3 <- wildBS_CB(formula=outcome ~ group + var_mat + var_group_mat,
                data=sumfun_data_full,
                spatial_covars = c('var_mat','var_group_mat'),
                target = c("Intercept","group", "var_mat", "var_group_mat"),
                form   = c("coef","term","term","term"),
                covar_list = covar_list_trt_stroma,
                func = function(Intercept, group, var_mat, var_group_mat){
                  Intercept + group + var_mat + var_group_mat
                },
                B1=B1,
                B2=B2,
                alpha=.05,
                re=NULL,
                id='patient_id',
                nthreads = 1)
save(CB3, file='CB3.RData')
plot_wildBS_CB(CB3)

#-------------- Get confidence band for E[Y|X] for a ctrl patient (stroma region only) ---------------------------------

covar_list_ctrl_stroma <- list(Intercept = rep(1, length(sumfun_data_ctrlpatient$group)),
                                 group = sumfun_data_ctrlpatient$group,
                                 var_mat = rep(0, length(sumfun_data_ctrlpatient$group)),
                                 var_group_mat = rep(0, length(sumfun_data_ctrlpatient$group)))

CB4 <- wildBS_CB(formula=outcome ~ group + var_mat + var_group_mat,
                data=sumfun_data_full,
                spatial_covars = c('var_mat','var_group_mat'),
                target = c("Intercept","group", "var_mat", "var_group_mat"),
                form   = c("coef","term","term","term"),
                covar_list = covar_list_ctrl_stroma,
                func = function(Intercept, group, var_mat, var_group_mat){
                  Intercept + group + var_mat + var_group_mat
                },
                B1=B1,
                B2=B2,
                alpha=.05,
                re=NULL,
                id='patient_id',
                nthreads = 1)
save(CB4, file='CB4.RData')
plot_wildBS_CB(CB4)

#-------------- Plot results---------------------------------------------------------------------------------------------
lw <- 1

pimage <- ggplot() +
  theme_bw() +
  geom_line(aes(x = CB1$grid, y = CB1$target_curve_estimates, col = "Estimated E[Y(r) | X]", linetype = "Treatment Group"), linewidth = lw) +
  geom_line(aes(x = CB2$grid, y = CB2$target_curve_estimates, col = "Estimated E[Y(r) | X]", linetype = "Control Group"), linewidth = lw) +
  geom_line(aes(x = CB1$grid, y = CB1$CB_lower, col = "95% Confidence Band", linetype = "Treatment Group"), linewidth = lw) +
  geom_line(aes(x = CB2$grid, y = CB2$CB_lower, col = "95% Confidence Band", linetype = "Control Group"), linewidth = lw) +
  geom_line(aes(x = CB1$grid, y = CB1$CB_upper, col = "95% Confidence Band", linetype = "Treatment Group"), linewidth = lw) +
  geom_line(aes(x = CB2$grid, y = CB2$CB_upper, col = "95% Confidence Band", linetype = "Control Group"), linewidth = lw) +
  labs(
    x = "r",
    y = expression(E * "[" * Y(r) * "|" * X * "]"),
    col = ""
  ) +
  theme(
    axis.title.y = element_text(size = 20),
    axis.title.x = element_text(size = 20)
  ) +
  scale_color_manual(
    values = c(
      "Estimated E[Y(r) | X]" = "#F8766D",
      "95% Confidence Band" = "#00BA38"
    )
  ) +
  scale_linetype_manual(name = "", values = c("Treatment Group" = 1, "Control Group" = 2)) +
  geom_hline(yintercept = 0, lty = 3) +
  ggtitle("Overall Image Colocalization Curve") +
  ylim(c(-1, 4)) +
  theme(
    title = element_text(size = 15),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 15)
  ) +
  theme(
    legend.key.width = unit(1, "cm")
  )

pstroma <- ggplot() +
  theme_bw() +
  geom_line(aes(x = CB3$grid, y = CB3$target_curve_estimates, col = "Estimated E[Y(r) | X]", linetype = "Treatment Group"), linewidth = lw) +
  geom_line(aes(x = CB4$grid, y = CB4$target_curve_estimates, col = "Estimated E[Y(r) | X]", linetype = "Control Group"), linewidth = lw) +
  geom_line(aes(x = CB3$grid, y = CB3$CB_lower, col = "95% Confidence Band", linetype = "Treatment Group"), linewidth = lw) +
  geom_line(aes(x = CB4$grid, y = CB4$CB_lower, col = "95% Confidence Band", linetype = "Control Group"), linewidth = lw) +
  geom_line(aes(x = CB3$grid, y = CB3$CB_upper, col = "95% Confidence Band", linetype = "Treatment Group"), linewidth = lw) +
  geom_line(aes(x = CB4$grid, y = CB4$CB_upper, col = "95% Confidence Band", linetype = "Control Group"), linewidth = lw) +
  labs(
    x = "r",
    y = expression(E * "[" * Y(r) * "|" * X * "]"),
    col = ""
  ) +
  theme(
    legend.position = "none",
    axis.title.y = element_text(size = 20),
    axis.title.x = element_text(size = 20)
  ) +
  scale_color_manual(
    values = c(
      "Estimated E[Y(r) | X]" = "#F8766D",
      "95% Confidence Band" = "#00BA38"
    )
  ) +
  scale_linetype_manual(name = "", values = c("Treatment Group" = 1, "Control Group" = 2)) +
  geom_hline(yintercept = 0, lty = 3) +
  ggtitle("Colocalization Curve in Stroma Regions") +
  ylim(c(-2, 3)) +
  theme(
    title = element_text(size = 15),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 15)
  ) +
  theme(
    legend.key.width = unit(1, "cm")
  )

png('sim2-bands.png', width=4000, height=1000, res=300)
(pimage + pstroma) + plot_annotation(
  theme = theme(plot.title = element_text(size = 22),
                plot.subtitle = element_text(size = 22))
)
dev.off()


#---Plot of Data Generation-------------------------------------------------------------------------------------

png(paste0('simplot5_1.png'), width=1700, height=1500, res=300)
ggplot(data=cellExp[cellExp$imageID=='5_1',]) + theme_bw() +
  geom_point(aes(x=x, y=y, col=cellType)) + labs(title='Representative Treatment Group Image',
                                                 x=expression('X Coordinate (' ~ mu ~ 'm)'),y=expression('Y Coordinate (' ~ mu ~ 'm)'), col='Cell Type') +
  theme(plot.title = element_text(size = 14.5, face = "bold"),
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=15), #change legend text font size
        axis.title=element_text(size=14)) +
  guides(colour = guide_legend(override.aes = list(size=3)))
dev.off()


png(paste0('simplot59_1.png'), width=1700, height=1500, res=300)
ggplot(data=cellExp[cellExp$imageID=='59_1',]) + theme_bw() +
  geom_point(aes(x=x, y=y, col=cellType)) + labs(title='Representative Control Group Image',
                                                 x=expression('X Coordinate (' ~ mu ~ 'm)'),y=expression('Y Coordinate (' ~ mu ~ 'm)'), col='Cell Type')  +
  theme(plot.title = element_text(size = 14.5, face = "bold"),
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=15), #change legend text font size
        axis.title=element_text(size=14)) +
  guides(colour = guide_legend(override.aes = list(size=3)))
dev.off()


## Then annotate in powerpoint
## Then convert to 300 dpi .png:

library(pdftools)
library(magick)

img <- image_read_pdf("annotate.pdf", density = 300)
image_write(img, "annotated.png")
