library(ggplot2)
library(patchwork)

load('melanoma_data.rda')

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
    ggplot2::scale_color_manual(values = c("Estimate" = "#F8766D",
                                           "Wild Bootstrap CB" = "#00BA38"))
  return(p)
}

#------------ CB for overall -------------------------------------------------------------------------------------
set.seed(12345)

sumfun_data <- preprocess_data(data=melanoma_data,
                               from_cell="CD8- T cell",
                               to_cell="CD8+ T cell",
                               qc_cellcount_cutoff=20,
                               n_perm=100,
                               perm_yn=TRUE,
                               r_max=500,
                               inc=1,
                               image_dims=c(0,1200,0,1200),
                               summary_function='L',
                               verbose=TRUE)
sumfun_data$patient_id <- as.factor(sumfun_data$patient_id)
save(sumfun_data, file='sumfun_data.RData')

CB_overall <- wildBS_CB(formula=outcome ~ s(patient_id, bs = "re"),
                       data=sumfun_data,
                       spatial_covars = NULL,
                       target = c("Intercept"),
                       form   = c("coef"),
                       covar_list = list(),
                       func = function(Intercept){
                         Intercept
                       },
                       B1=50,
                       B2=50,
                       alpha=.05,
                       re='patient_id',
                       nthreads = 1,
                       id='patient_id')
save(CB_overall, file='CB_overall.RData')
plot_wildBS_CB(CB_overall)

#------------ CB for avascular regions -------------------------------------------------------------------------------------
set.seed(12345)

#var_mat represents the values of the spatial covariate \hat{L}_{CD8-,Vas}(r)
var_mat_data <- preprocess_data(melanoma_data,
                                from_cell='CD8- T cell',
                                to_cell='Vasculature',
                                qc_cellcount_cutoff=20,
                                n_perm=100,
                                perm_yn=TRUE,
                                r_max=500,
                                inc=1,
                                image_dims=c(0,1200,0,1200),
                                summary_function='L',
                                verbose=TRUE)
var_mat_data$patient_id <- as.factor(var_mat_data$patient_id)
names(var_mat_data)[names(var_mat_data)=='L_obs'] <- 'var_mat'
var_mat_data$L_expect <- var_mat_data$L_pmean <- var_mat_data$outcome <- NULL
sumfun_data_full <- merge(sumfun_data, var_mat_data, by=c('image_number','r','patient_id','patient_age'),all=T)
sumfun_data_full$var_mat <- ifelse(is.na(sumfun_data_full$var_mat), 0, sumfun_data_full$var_mat)

covar_list <- list(Intercept = rep(1, 501),
                   var_mat = rep(0, 501))

CB_avascular <- wildBS_CB(formula=outcome ~ var_mat + s(patient_id, bs = "re"),
                         data=sumfun_data_full,
                         spatial_covars = c('var_mat'),
                         target = c("Intercept","var_mat"),
                         form   = c("coef","term"),
                         covar_list = covar_list,
                         func = function(Intercept, var_mat){
                           Intercept + var_mat
                         },
                         B1=50,
                         B2=50,
                         alpha=.05,
                         re='patient_id',
                         nthreads = 1,
                         id='patient_id')
save(CB_avascular, file='CB_avascular.RData')
plot_wildBS_CB(CB_avascular)

#-------------- Plot results---------------------------------------------------------------------------------------------

pimage <- ggplot() + theme_bw() +
  geom_line(aes(x=CB_overall$grid,y=CB_overall$target_curve_estimates, col='Estimated E[Y(r)|X]'), size=1.5) +
  geom_line(aes(x=CB_overall$grid,y=CB_overall$CB_lower, col='95% Confidence Band'), size=1.5) +
  geom_line(aes(x=CB_overall$grid,y=CB_overall$CB_upper, col='95% Confidence Band'), size=1.5) +
  labs(x='r', y='E[Y(r)|X]' #expression(
    #hat(L)^{CD8 * "-," * CD8 * "+,"} * "(" * r * ")" -
    #  hat(L)^{CD8 * "-," * CD8 * "+,"} * "(" * r * ")"
  #),
  #col = " "
  ) +
  theme(legend.position = 'none',
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20))  +
  scale_color_manual(values = c("Estimated E[Y(r)|X]" = '#F8766D',
                                "95% Confidence Band"='#00BA38')) +
  geom_hline(yintercept=0, lty=3, size=1) + ylim(c(-100,200)) +
  theme(title = element_text(size = 18),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        plot.subtitle = element_text(size=18)) +
  ggtitle('Overall Image')

pavascular <- ggplot() + theme_bw() +
  geom_line(aes(x=CB_avascular$grid,y=CB_avascular$target_curve_estimates, col='Estimated E[Y(r)|X]'), size=1.5) +
  geom_line(aes(x=CB_avascular$grid,y=CB_avascular$CB_lower, col='95% Confidence Band'), size=1.5) +
  geom_line(aes(x=CB_avascular$grid,y=CB_avascular$CB_upper, col='95% Confidence Band'), size=1.5) +
  labs(x='r', y='E[Y(r)|X]' #expression(
   # hat(L)^{CD8 * "-," * CD8 * "+,"} * "(" * r * ")" -
  #    hat(L)^{CD8 * "-," * CD8 * "+,"} * "(" * r * ")"
  #),
  #col = " "
  ) +
  theme(
    axis.title.y = element_text(size = 20),
    axis.title.x = element_text(size = 20))  +
  scale_color_manual(name = " ",
                     values = c("Estimated E[Y(r)|X]" = '#F8766D',
                                "95% Confidence Band"='#00BA38')) +
  geom_hline(yintercept=0, lty=3, size=1) + ylim(c(-100,200)) +
  theme(title = element_text(size = 18),
        legend.title = element_text(size=12),
        legend.text = element_text(size=18),
        plot.subtitle = element_text(size=18)) +
  ggtitle('Avascular Regions')

png('melanoma-bands.png', width=4000, height=1050, res=300)
(pimage + pavascular) + plot_annotation( #title = 'Colocalization Curve',
  theme = theme(plot.title = element_text(size = 22),
                plot.subtitle = element_text(size = 22))
)
dev.off()
#
# data <- melanoma_data
# cols <- c('Other' = "gray",
#           'CD8- T cell' = 'purple',
#           'CD8+ T cell' = "orange",
#           'Vasculature'='cyan4')
# data$cells_interested <- ifelse(data$cell_type %in% c('CD8- T cell',
#                                                       'CD8+ T cell',
#                                                       'Vasculature'), data$cell_type, 'Other')
# for(i in c(45, 54, 65, 137)){
#   plot.data <- data[data$image_number == i,]
#   plot.data$cells_interested <- factor(plot.data$cells_interested)
#   assign(paste0('p',i), ggplot(data=plot.data, aes(x = cell_x, y = cell_y, col=cells_interested)) +
#            geom_point() + theme_bw() + ggtitle(paste0('Image ',i)) +
#            labs(x=expression('X Coordinate (' ~ mu ~ 'm)'), y = expression('Y Coordinate (' ~ mu ~ 'm)')) +
#            scale_color_manual('Cell Type', values=cols, labels = c('CD8- T Cells', 'CD8+ T Cells', 'Other', 'Vasculature Cells') ) +
#            theme(title = element_text(size = 18),
#                  axis.title.x = element_text(size = 14),
#                  axis.title.y = element_text(size = 14),
#                  legend.key.height = unit(1, 'cm'), #change legend key height
#                  legend.key.width = unit(1, 'cm'), #change legend key width
#                  legend.title = element_text(size=16), #change legend title font size
#                  legend.text = element_text(size=16)) +
#            guides(colour = guide_legend(override.aes = list(size=5))) )
# }
# p45 <- p45 + theme(legend.position = 'none')
# p54 <- p54 + theme(legend.position = 'none')
# p65 <- p65 + theme(legend.position = 'none')
#
# png('melanoma-images.png', width=4800, height=1325, res=300)
# (p45 + p54 + p65 + p137) + plot_layout(ncol = 4) + plot_annotation(
#   title = "Interaction Between CD8+ T Cells and CD8- T Cells",
#   theme = theme(plot.title = element_text(size = 22, face = "bold"),
#                 plot.subtitle = element_text(size = 22))
# )
# dev.off()
