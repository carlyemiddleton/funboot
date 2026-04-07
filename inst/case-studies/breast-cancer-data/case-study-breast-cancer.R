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
library(grid)

#------------ Image plots ---------------------------------------------------------------------------
library(funboot)
data('breastcancer_data')

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


#------------ Interaction CB -------------------------------------------------------------------------------------
set.seed(12345)

sumfun_data_interaction <- preprocess_data(data=breastcancer_data,
                               from_cell=7,
                               to_cell=3,
                               qc_cellcount_cutoff=20,
                               n_perm=100,
                               perm_yn=TRUE,
                               r_max=500,
                               inc=1,
                               image_dims=c(0,1000,0,1000),
                               summary_function='L',
                               verbose=TRUE)
sumfun_data_interaction$patient_id <- as.factor(sumfun_data_interaction$patient_id)
save(sumfun_data_interaction, file='sumfun_data_interaction.RData')

CB_interaction <- wildBS_CB(formula=outcome ~ 1,
                        data=sumfun_data_interaction,
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
                        re=NULL,
                        nthreads = 1,
                        id='patient_id')
save(CB_interaction, file='CB_interaction.RData')
plot_wildBS_CB(CB_interaction)

#------------ Avoidance CB -------------------------------------------------------------------------------------
set.seed(12345)

sumfun_data_avoidance <- preprocess_data(data=breastcancer_data,
                               from_cell=9,
                               to_cell=16,
                               qc_cellcount_cutoff=20,
                               n_perm=100,
                               perm_yn=TRUE,
                               r_max=500,
                               inc=1,
                               image_dims=c(0,1000,0,1000),
                               summary_function='L',
                               verbose=TRUE)
sumfun_data_avoidance$patient_id <- as.factor(sumfun_data_avoidance$patient_id)
save(sumfun_data_avoidance, file='sumfun_data_avoidance.RData')

CB_avoidance <- wildBS_CB(formula=outcome ~ 1,
                            data=sumfun_data_avoidance,
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
                            re=NULL,
                            nthreads = 1,
                            id='patient_id')
save(CB_avoidance, file='CB_avoidance.RData')
plot_wildBS_CB(CB_avoidance)

# #-------------- Plot results---------------------------------------------------------------------------------------------
lw <- 1.3

pInt <- ggplot() +
  theme_bw() +
  geom_line(
    data = sumfun_data_interaction[order(sumfun_data_interaction$image_number,
                                         sumfun_data_interaction$r), ],
    aes(x = r, y = outcome, group = image_number, col = "Observed Y(r)"),
    linewidth = .7,
    alpha = 1
  ) +
  geom_line(
    aes(x = CB_interaction$grid,
        y = CB_interaction$target_curve_estimates,
        col = "Estimated E[Y(r)]"),
    linewidth = lw
  ) +
  geom_line(
    aes(x = CB_interaction$grid,
        y = CB_interaction$CB_lower,
        col = "95% Confidence Band"),
    linewidth = lw
  ) +
  geom_line(
    aes(x = CB_interaction$grid,
        y = CB_interaction$CB_upper,
        col = "95% Confidence Band"),
    linewidth = lw
  ) +
  labs(
    x = "r",
    y = "Y(r)",
    col = " "
  ) +
  theme(
    axis.title.y = element_text(size = 20),
    axis.title.x = element_text(size = 20)
  ) +
  scale_color_manual(values = c(
    "Observed Y(r)" = "grey",
    "Estimated E[Y(r)]" = "#F8766D",
    "95% Confidence Band" = "#00BA38"
  )) +
  geom_hline(yintercept = 0, lty = 3) +
  ylim(c(-155, 125)) +
  theme(
    title = element_text(size = 18),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16)
  ) +
  ggtitle("Colocalization Curve")

pAvoid <- ggplot() +
  theme_bw() +
  geom_line(
    data = sumfun_data_avoidance[order(sumfun_data_avoidance$image_number,
                                       sumfun_data_avoidance$r), ],
    aes(x = r, y = outcome, group = image_number, col = "Observed Y(r)"),
    linewidth = .7,
    alpha = 1
  ) +
  geom_line(
    aes(x = CB_avoidance$grid,
        y = CB_avoidance$target_curve_estimates,
        col = "Estimated E[Y(r)]"),
    linewidth = lw
  ) +
  geom_line(
    aes(x = CB_avoidance$grid,
        y = CB_avoidance$CB_lower,
        col = "95% Confidence Band"),
    linewidth = lw
  ) +
  geom_line(
    aes(x = CB_avoidance$grid,
        y = CB_avoidance$CB_upper,
        col = "95% Confidence Band"),
    linewidth = lw
  ) +
  labs(
    x = "r",
    y = "Y(r)",
    col = " "
  ) +
  theme(
    axis.title.y = element_text(size = 20),
    axis.title.x = element_text(size = 20)
  ) +
  scale_color_manual(values = c(
    "Observed Y(r)" = "grey",
    "Estimated E[Y(r)]" = "#F8766D",
    "95% Confidence Band" = "#00BA38"
  )) +
  geom_hline(yintercept = 0, lty = 3) +
  ylim(c(-150, 125)) +
  theme(
    title = element_text(size = 18),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16)
  ) +
  ggtitle("Colocalization Curve")

##Image plots
data <- breastcancer_data
cols <- c('0' = "gray", '7' = 'gold', '3' = "darkgreen", '9' = 'skyblue', '16'='cyan4')
data$cells_interested <- ifelse(data$cell_type %in% c(7,3), data$cell_type, 0)
for(i in c(362, 181, 246, 13)){
  plot.data <- data[data$image_number == i,]
  plot.data$cells_interested <- factor(plot.data$cells_interested)
  assign(paste0('p',i), ggplot(data=plot.data, aes(x = cell_x, y = cell_y, col=cells_interested)) +
           geom_point() + theme_bw() + ggtitle(paste0('Image ',i)) +
           labs(x=expression('X Coordinate (' ~ mu ~ 'm)'), y = expression('Y Coordinate (' ~ mu ~ 'm)')) +
           scale_color_manual('Cell Type', values=cols, labels = c('Other', 'Endothelial \nCells', 'T Cells') ) +
           theme(title = element_text(size = 18),
                 axis.title.x = element_text(size = 14),
                 axis.title.y = element_text(size = 14),
                 legend.key.height = unit(1, 'cm'),
                 legend.key.width = unit(1, 'cm'),
                 legend.title = element_text(size=16),
                 legend.text = element_text(size=16)) +
           guides(colour = guide_legend(override.aes = list(size=5))) )
}
p13 <- p13 + theme(legend.position = 'none')
p181 <- p181 + theme(legend.position = 'none')
p246_top <- p246 + theme(legend.position = 'none')

data$cells_interested <- ifelse(data$cell_type %in% c(9,16), data$cell_type, 0)
for(i in c(219, 244, 246, 254)){
  plot.data <- data[data$image_number == i,]
  plot.data$cells_interested <- factor(plot.data$cells_interested)
  assign(paste0('p',i), ggplot(data=plot.data, aes(x = cell_x, y = cell_y, col=cells_interested)) +
           geom_point() + theme_bw() + ggtitle(paste0('Image ',i)) +
           labs(x=expression('X Coordinate (' ~ mu ~ 'm)'), y = expression('Y Coordinate (' ~ mu ~ 'm)')) +
           scale_color_manual('Cell Type', values=cols, labels = c('Other', 'Small Circular \nStromal Cells', 'Proliferative \nEpithelial Cells') ) +
           theme(title = element_text(size = 18),
                 axis.title.x = element_text(size = 14),
                 axis.title.y = element_text(size = 14),
                 legend.key.height = unit(1, 'cm'), #change legend key height
                 legend.key.width = unit(1, 'cm'), #change legend key width
                 legend.title = element_text(size=16), #change legend title font size
                 legend.text = element_text(size=16)) +
           guides(colour = guide_legend(override.aes = list(size=5))) )
}
 p219 <- p219 + theme(legend.position = 'none')
 p244 <- p244 + theme(legend.position = 'none')
 p246_bottom <- p246 + theme(legend.position = 'none')

top_images <- p13 + p181 + p246_top + p362 +
  plot_layout(ncol = 4)

bottom_images <- p219 + p244 + p246_bottom + p254 +
  plot_layout(ncol = 4)

top_row <- (top_images | pInt) +
  plot_layout(widths = c(5, 1)) +
  plot_annotation(
    title = "Interaction Between Endothelial Cells and T Cells",
    theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0))
  )

bottom_row <- (bottom_images | pAvoid) +
  plot_layout(widths = c(5, 1)) +
  plot_annotation(
    title = "Avoidance Between Small Circular Stromal Cells and Proliferative Epithelial Cells",
    theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0))
  )

png('case-study-breast-cancer_top.png', width=5500*1.3, height=1000*1.3, res=300)
top_row
dev.off()
png('case-study-breast-cancer_bottom.png', width=5500*1.3, height=1000*1.3, res=300)
bottom_row
dev.off()


#-------------- Model with covariates ---------------------------------------------------------------------------------------------
set.seed(12345)
##Compute \hat{L}_{Endo,T}(r), obtaining data from 76 images (image 249 has 0 type 3 cells, and 23 images were removed during qc)
sumfun_data_covariates <- funboot::preprocess_data(data=breastcancer_data,
                                                   from_cell=7, #Endothelial cells
                                                   to_cell=3,   #T cells
                                                   qc_cellcount_cutoff=20,
                                                   n_perm=100,
                                                   perm_yn=TRUE,
                                                   r_max=500,
                                                   inc=1,
                                                   image_dims=c(0,1000,0,1000),
                                                   summary_function='L',
                                                   verbose=TRUE)
##Compute \hat{L}_{Endo,Endo}(r), obtaining data from 48 images (52 images were removed during qc)
spatialcov_data <- funboot::preprocess_data(data=breastcancer_data,
                                            from_cell=7, #Endothelial cells
                                            to_cell=7,   #Endothelial cells
                                            qc_cellcount_cutoff=20,
                                            n_perm=100,
                                            perm_yn=TRUE,
                                            r_max=500,
                                            inc=1,
                                            image_dims=c(0,1000,0,1000),
                                            summary_function='L',
                                            verbose=TRUE)
#Retain just \hat{L}_{Endo,Endo}(r) and no other L curves
names(spatialcov_data)[names(spatialcov_data)=='L_obs'] <- 'spatial_cov'
spatialcov_data$L_expect <- spatialcov_data$L_pmean <- spatialcov_data$outcome <- NULL
##Merge the two datasets together, resulting in data from 48 images
sumfun_data_full_covariates <- merge(sumfun_data_covariates,
                                     spatialcov_data,
                                     by=c('image_number','r','patient_id','patient_age','tumor_grade'),
                                     all=FALSE)
##Add in indicator variables for grade
sumfun_data_full_covariates$grade2 <- ifelse(sumfun_data_full_covariates$tumor_grade == '2', 1,
                                             ifelse(is.na(sumfun_data_full_covariates$tumor_grade), NA,
                                                    0))
sumfun_data_full_covariates$grade3 <- ifelse(sumfun_data_full_covariates$tumor_grade == '3', 1,
                                             ifelse(is.na(sumfun_data_full_covariates$tumor_grade), NA,
                                                    0))

# set.seed(12345)
#
# sumfun_data_covariates <- preprocess_data(data=breastcancer_data,
#                                            from_cell=7,
#                                            to_cell=3,
#                                            qc_cellcount_cutoff=20,
#                                            n_perm=100,
#                                            perm_yn=T,
#                                            r_max=500,
#                                            inc=1,
#                                            image_dims=c(0,1000,0,1000),
#                                            summary_function='L',
#                                            verbose=TRUE)
# #combine cell types 3,4,5,6 into 1 category
# #breastcancer_data$cell_type <- ifelse(breastcancer_data$cell_type %in% c(3:6), 3,
# #                                                       breastcancer_data$cell_type)
# spatialcov_data <- preprocess_data(data=breastcancer_data,
#                                    from_cell=7,
#                                    to_cell=7,
#                                    qc_cellcount_cutoff=20,
#                                    n_perm=100,
#                                    perm_yn=T,
#                                    r_max=500,
#                                    inc=1,
#                                    image_dims=c(0,1000,0,1000),
#                                    summary_function='L',
#                                    verbose=TRUE)
# names(spatialcov_data)[names(spatialcov_data)=='L_obs'] <- 'spatial_cov'
# spatialcov_data$L_expect <- spatialcov_data$L_pmean <- spatialcov_data$outcome <- NULL
# sumfun_data_full_covariates <- merge(sumfun_data_covariates,
#                                      spatialcov_data,
#                                      by=c('image_number','r','patient_id','patient_age','tumor_grade'),
#                                      all=T)
# df <- data.frame(image_number = rep(unique(breastcancer_data$image_number), each=501),
#                  r = rep(0:500, length(unique(breastcancer_data$image_number)))
#                  )
# sumfun_data_full_covariates <- merge(df,
#                                      sumfun_data_full_covariates,
#                                      by=c('image_number','r'),
#                                      all=T)
# sumfun_data_full_covariates <- sumfun_data_full_covariates[!is.na(sumfun_data_full_covariates$outcome),] #Drop obs. with missing outcomes
# sumfun_data_full_covariates$spatial_cov <- ifelse(is.na(sumfun_data_full_covariates$spatial_cov), 0,
#                                                   sumfun_data_full_covariates$spatial_cov) #set L^{3,16}(r) = 0 when there are no type 16 cells
# sumfun_data_full_covariates$grade2 <- ifelse(sumfun_data_full_covariates$tumor_grade == '2', 1,
#                                              ifelse(is.na(sumfun_data_full_covariates$tumor_grade), NA,
#                                                     0))
# sumfun_data_full_covariates$grade3 <- ifelse(sumfun_data_full_covariates$tumor_grade == '3', 1,
#                                              ifelse(is.na(sumfun_data_full_covariates$tumor_grade), NA,
#                                                     0))
# save(sumfun_data_full_covariates, file='sumfun_data_full_covariates.RData')

##Get CBs for each coefficient
CB_beta0 <- wildBS_CB(formula=outcome ~ patient_age + grade2 + grade3 + spatial_cov,
                      data=sumfun_data_full_covariates,
                      spatial_covars = c("spatial_cov"),
                      target = c("Intercept"),
                      form   = c("coef"),
                      covar_list = list(),
                      func = function(Intercept){
                        Intercept
                      },
                      B1=50,
                      B2=50,
                      alpha=.05,
                      re=NULL,
                      nthreads = 1,
                      id='patient_id')
save(CB_beta0, file='CB_beta0.RData')
M <- max(abs((CB_beta0$target_curve_estimates - 0)/CB_beta0$target_curve_SEs))
mean(CB_beta0$M>M) #p-value for the test of beta0(r)=0

CB_beta1 <- wildBS_CB(formula=outcome ~ patient_age + grade2 + grade3 + spatial_cov,
                      data=sumfun_data_full_covariates,
                      spatial_covars = c("spatial_cov"),
                      target = c("patient_age"),
                      form   = c("coef"),
                      covar_list = list(),
                      func = function(patient_age){
                        patient_age
                      },
                      B1=50,
                      B2=50,
                      alpha=.05,
                      re=NULL,
                      nthreads = 1,
                      id='patient_id')
save(CB_beta1, file='CB_beta1.RData')
M <- max(abs((CB_beta1$target_curve_estimates - 0)/CB_beta1$target_curve_SEs))
mean(CB_beta1$M>M) #p-value for the test of beta1(r)=0

CB_beta2 <- wildBS_CB(formula=outcome ~ patient_age + grade2 + grade3 + spatial_cov,
                      data=sumfun_data_full_covariates,
                      spatial_covars = c("spatial_cov"),
                      target = c("grade2"),
                      form   = c("coef"),
                      covar_list = list(),
                      func = function(grade2){
                        grade2
                      },
                      B1=50,
                      B2=50,
                      alpha=.05,
                      re=NULL,
                      nthreads = 1,
                      id='patient_id')
save(CB_beta2, file='CB_beta2.RData')
M <- max(abs((CB_beta2$target_curve_estimates - 0)/CB_beta2$target_curve_SEs))
mean(CB_beta2$M>M) #p-value for the test of beta2(r)=0

CB_beta3 <- wildBS_CB(formula=outcome ~ patient_age + grade2 + grade3 + spatial_cov,
                      data=sumfun_data_full_covariates,
                      spatial_covars = c("spatial_cov"),
                      target = c("grade3"),
                      form   = c("coef"),
                      covar_list = list(),
                      func = function(grade3){
                        grade3
                      },
                      B1=50,
                      B2=50,
                      alpha=.05,
                      re=NULL,
                      nthreads = 1,
                      id='patient_id')
save(CB_beta3, file='CB_beta3.RData')
M <- max(abs((CB_beta3$target_curve_estimates - 0)/CB_beta3$target_curve_SEs))
mean(CB_beta3$M>M) #p-value for the test of beta3(r)=0

CB_beta4 <- wildBS_CB(formula=outcome ~ patient_age + grade2 + grade3 + spatial_cov,
                      data=sumfun_data_full_covariates,
                      spatial_covars = c("spatial_cov"),
                      target = c("spatial_cov"),
                      form   = c("coef"),
                      covar_list = list(),
                      func = function(spatial_cov){
                        spatial_cov
                      },
                      B1=50,
                      B2=50,
                      alpha=.05,
                      re=NULL,
                      nthreads = 1,
                      id='patient_id')
save(CB_beta4, file='CB_beta4.RData')
M <- max(abs((CB_beta4$target_curve_estimates - 0)/CB_beta4$target_curve_SEs))
mean(CB_beta4$M>M) #p-value for the test of beta4(r)=0
