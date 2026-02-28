library(parallel)

#------------ Setup ----------------------------------------------------------------------------------------------
nCores <- 2
nthreads <- 1
B1   <- 100
B2   <- 100
alpha <- 0.05
seed <- 12345

nSim <- 500
nPatients <- 50  #number of patients
nIm <- 3         #number of images per patient
sigmaB <- 4
sigmaB_diff <- 4*2/3

progress_file <- "sim_progress.RData"

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

extract_package_SEs <- function(name, model, grid) {
  sm <- stats::coef(model, n1 = length(grid))$smterms
  if (name %in% names(sm)) result <- sm[[name]]$coef$se
  ind <- grep(name, names(sm), value = TRUE)
  if (length(ind) == 1) result <- sm[[ind]]$coef$se
  if (length(ind) > 1) stop("Multiple smooth terms match for ", name)
  if (length(ind) == 0)stop("No smooth term found matching ", name)
  return(result)
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
                      B2=NULL,
                      alpha=.05,
                      re=NULL,
                      id,
                      nthreads = 1,
                      inner_SEs_from_package=FALSE){

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

  ## Get package SE's
  SE_pkg_orig <- NULL
  if (inner_SEs_from_package==TRUE) {
    if (length(target) != 1) {
      stop("inner_SEs_from_package cannot be used when target has length greater than 1.")
    }
    if (any(form != "coef")) {
      stop("inner_SEs_from_package can only be used if form = 'coef'. ")
    }
    if (n_patients != n_Im_total) {
      stop("inner_SEs_from_package cannot be used with multiple images per subject.")
    }
    if (!is.null(re)) {
      stop("inner_SEs_from_package cannot be used with random intercept models.")
    }

    SE_pkg_orig <- extract_package_SEs(name = target, model = pffrmodel, grid = grid)
    SE_pkg_orig <- pmax(SE_pkg_orig, 1e-8) # avoid division by 0
  }


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

    ##If inner_SEs_from_package == TRUE, skip the inner layer and use pffr()'s package SEs instead:
    if (inner_SEs_from_package==TRUE) {
      SEs_inner <- extract_package_SEs(name = target, model = pffrmodel_bs, grid = grid)
      SEs_inner <- pmax(SEs_inner, 1e-8) #avoid potential division by 0
      M[b1] <- max(abs((target_curve_bs - target_curve) / SEs_inner))
      target_curve_outer[b1, ] <- target_curve_bs
    }else{
      ##Do the inner layer:
      if (!inner_SEs_from_package && is.null(B2)) {
        stop("B2 must be provided when inner_SEs_from_package = FALSE.")
      }
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
  }

  #########################
  ## Get confidence band ##
  #########################

  q <- stats::quantile(M, 1-alpha)

  if(inner_SEs_from_package){
    target_curve_SEs <- SE_pkg_orig
    MoE <- q * target_curve_SEs
  } else {
    target_curve_SEs <- apply(target_curve_outer, 2, stats::sd)
    MoE <- q * target_curve_SEs
  }
  CB_lower <- target_curve - MoE
  CB_upper <- target_curve + MoE
  CBs <- list(CB_lower, CB_upper, grid, target_curve, target_curve_SEs, q, M)
  names(CBs) <- c('CB_lower','CB_upper', 'grid', 'target_curve_estimates', 'target_curve_SEs', 'q', 'M')
  return(CBs)
}

#------------ Define simulation functions -------------------------------------------------------------------------
#Function to simulate 1 draw of data and calculate CB. To be repeated nSim times:
sim <- function(iter,
                nPatients,
                nIm,
                sigmaB,
                sigmaB_diff,
                grid,
                alpha,
                B1,
                B2,
                nthreads,
                seed){

  set.seed(iter+seed)

  ##Settings for data generation
  x <- y <- cellType <- imageID <- NULL
  #Rhard_G1  <- .05
  #beta_G1   <- 40
  #muB_G1    <- 8
  sigmaB_G1 <- sigmaB

  #Rhard_G2  <- .05
  #beta_G2   <- 40
  #muB_G2    <- 8
  sigmaB_G2 <- sigmaB + sigmaB_diff

  ##Generate data
  for (p in 1:nPatients){
    group <- ifelse(p <= nPatients/2, "Group1", "Group2")
    if (group == "Group1") {
      #Rhard   <- Rhard_G1
      #beta_p  <- beta_G1
      #muB     <- muB_G1
      sigmaB  <- abs(rnorm(1, mean = sigmaB_G1, sd = sqrt(sigmaB_G1 / 10)))
    } else {
      #Rhard   <- Rhard_G2
      #beta_p  <- beta_G2
      #muB     <- muB_G2
      sigmaB  <- abs(rnorm(1, mean = sigmaB_G2, sd = sqrt(sigmaB_G2 / 10)))
    }
    nIm <- nIm
    for (j in 1:nIm){
      Rhard  <- sample(c(0.02, 0.03,  0.04,  0.05, 0.06), 1)
      beta_p <- sample(seq(50, 100, by=10), 1)
      muB    <- sample(seq(10, 30, by=1), 1)
      A_unit <- spatstat.random::rHardcore(
        beta = beta_p,
        R    = Rhard,
        W    = owin(c(0,1),c(0,1))
      )
      A_x <- A_unit$x * 1000
      A_y <- A_unit$y * 1000
      nA  <- A_unit$n
      B_xy <- matrix(numeric(0), ncol = 2)
      if (nA > 0){
        for (i in 1:nA){
          kB <- rpois(1, muB)
          if (kB > 0) {
            ptsB <- cbind(
              rnorm(kB, A_x[i], sigmaB),
              rnorm(kB, A_y[i], sigmaB)
            )
            B_xy <- rbind(B_xy, ptsB)
          }
        }
      }
      #remove B cells generated outside of the window
      if (nrow(B_xy) > 0) {
        inside <- B_xy[,1] >= 0 & B_xy[,1] <= 1000 &
          B_xy[,2] >= 0 & B_xy[,2] <= 1000
        B_xy <- B_xy[inside, , drop = FALSE]
      }
      ID <- paste(p, j, sep = "_")
      if (nA > 0) {
        x <- c(x, A_x)
        y <- c(y, A_y)
        cellType <- c(cellType, rep("A", nA))
        imageID  <- c(imageID, rep(ID, nA))
      }
      if (nrow(B_xy) > 0) {
        x <- c(x, B_xy[,1])
        y <- c(y, B_xy[,2])
        cellType <- c(cellType, rep("B", nrow(B_xy)))
        imageID  <- c(imageID, rep(ID, nrow(B_xy)))
      }
    }
  }
  imageID <- factor(imageID)

  #build table
  cellExp <- data.frame(
    x = x,
    y = y,
    cellType = factor(cellType),
    imageID = imageID
  )
  phenoData <- data.frame(
    imageID  = unique(imageID),
    condition = ifelse(as.numeric(sapply(unique(imageID),function(x) str_split(x, "_")[[1]][1])) <= nPatients/2,
                        "Group1",
                        "Group2"),
    subject   = as.numeric(sapply(unique(imageID),function(x) str_split(x, "_")[[1]][1]))
  )
  table <- cellExp
  names(table) <- c("Xcoord", "Ycoord", "cell.type", "image")
  table$patient_id <- sapply(table$image, function(x) str_split(x, "_")[[1]][1])
  table$roi        <- sapply(table$image, function(x) str_split(x, "_")[[1]][2])
  names(phenoData) <- c("image", "group", "patient_id")
  table <- merge(table, phenoData,
                 by = c("image", "patient_id"), all = TRUE)
  colnames(table) <- c("image_number", "patient_id", "cell_x", "cell_y",
                       "cell_type", "roi", "group")

  ##Calculate spatial summary functions
  sumfun_data <- preprocess_data(
    table,
    from_cell = "A",
    to_cell = "B",
    qc_cellcount_cutoff = 0,
    perm_yn = FALSE, #no holes in simulated data
    r_max = 200,
    inc = 1,
    image_dims = c(0, 1000, 0, 1000),
    summary_function = "L"
  )
  sumfun_data$group <- ifelse(sumfun_data$group == "Group1", 1, 0)
  sumfun_data$patient_id <- factor(sumfun_data$patient_id)

  ##Calculate confidence band
  CBs <- wildBS_CB(formula=outcome ~ group + s(patient_id, bs = "re"),
                   data=sumfun_data,
                   target = c("group"),
                   form   = c("coef"),
                   func = function(group){group},
                   B1=B1,
                   B2=B2,
                   alpha=.05,
                   re='patient_id',
                   nthreads = nthreads,
                   id='patient_id',
                   inner_SEs_from_package=FALSE)
  excludes0 <- any(CBs$CB_lower > 0 | CBs$CB_upper < 0)
  list(
    excludes0 = excludes0,
    target_curve_estimates = CBs$target_curve_estimates,
    CB_lower = CBs$CB_lower,
    CB_upper = CBs$CB_upper
  )
}

sim_with_save <- function(iter,
                              nPatients,
                              nIm,
                              sigmaB,
                              sigmaB_diff,
                              grid,
                              alpha,
                              B1,
                              B2,
                              nthreads,
                              seed,
                              progress_file,
                              nSim){
  res <- sim(iter, nPatients, nIm, sigmaB, sigmaB_diff, grid, alpha, B1, B2, nthreads, seed)

  #append progress
  attempt <- 1
  repeat {
    ok <- try({
      if (file.exists(progress_file)) {
        load(progress_file)
        if (length(excludes0_vec) != nSim) {
          excludes0_vec <- rep(NA, nSim)
        }
      } else {
        excludes0_vec <- rep(NA, nSim)
      }
      excludes0_vec[iter] <- res$excludes0
      save(excludes0_vec, file = progress_file)
      TRUE
    }, silent = TRUE)
    if (isTRUE(ok)) break
    if (attempt > 10) {
      warning("Failed to update progress after 10 attempts.")
      break
    }
    attempt <- attempt + 1
    Sys.sleep(runif(1, 0.05, 0.2))
  }

  res
}

#------------ Simulate--------------------------------------------------------------------------------------------
set.seed(seed)
iters <- as.list(1:nSim)
if (file.exists(progress_file)) file.remove(progress_file)
cl <- makeCluster(nCores, type = "PSOCK")
clusterEvalQ(cl, { #load the necessary packages on the parallel workers
  library(spatstat.geom)
  library(spatstat.random)
  library(spatstat.explore)
  library(stringr)
  library(refund)
  library(mgcv)
})
clusterExport(cl, #export the necessary variables and functions to the parallel workers
              varlist = c(
                          #Variables
                          "nPatients", "nIm", "sigmaB", "sigmaB_diff",
                          "grid", "alpha", "B1", "B2", "nthreads", "seed",
                          "nSim", "progress_file",

                          #Functions
                          "format_spatial_variable", "preprocess_data",
                          "extract_target_curve", "extract_package_SEs", "get_fixed_and_re",
                          "wildBS_CB",
                          "sim", "sim_with_save"
                          ),
              envir = environment()
)

res_list <- parLapplyLB(cl, iters, sim_with_save,
                        nPatients = nPatients, nIm=nIm, sigmaB = sigmaB, sigmaB_diff=sigmaB_diff,
                        grid = grid, alpha = alpha, B1 = B1, B2 = B2,
                        nthreads = nthreads, seed=seed,
                        progress_file = progress_file, nSim = nSim)

target_curve <- do.call(rbind, lapply(res_list, `[[`, "target_curve_estimates"))
CB_lower     <- do.call(rbind, lapply(res_list, `[[`, "CB_lower"))
CB_upper     <- do.call(rbind, lapply(res_list, `[[`, "CB_upper"))
excludes0    <- sapply(res_list, `[[`, "excludes0") #does the CB exclude 0?
grid         <- 0:200

save(target_curve, CB_lower, CB_upper, excludes0, grid, file = "sim_results.RData")



















