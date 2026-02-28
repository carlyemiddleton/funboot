#' @title Calculates wild bootstrap confidence bands for *pffr()* model coefficients
#'
#' @param formula An object of class *formula*.  The formula for the *pffr()* model with special terms as in *mgcv's* GAM
#' @param data Data frame obtained from the output of *funboot::preprocess_data()*
#' @param spatial_covars Vector containing the names of the covariate(s) to be treated as spatially-varying.  For example, c('var1','var2')
#' @param target Number of outer bootstrap samples
#' @param form Number of outer bootstrap samples
#' @param covar_list Number of outer bootstrap samples
#' @param func Number of outer bootstrap samples
#' @param B1 Number of outer bootstrap samples
#' @param B2 Number of inner bootstrap samples
#' @param alpha Desired significance level for the confidence band. For example, 0.05
#' @param re Vector containing names of variable(s) for which to fit a random intercept.  For example, c('patient_id')
#' @param id Vector containing names of variable(s) for which to fit a random intercept.  For example, c('patient_id')
#' @param nthreads Vector containing names of variable(s) for which to fit a random intercept.  For example, c('patient_id')
#' @param inner_SEs_from_package Whether or not to replace the inner bootstrap loop standard errors with *pffr()*'s estimated standard errors of the coefficient function.  Should only be used if the target quantity for the confidence band is a single model coefficient and the data contain one image per patient.
#'
#' @return A list containing the lower and upper bounds of the confidence band, as well as the estimated model coefficients and standard errors from the *pffr()* output and additional intermittent variables used in the calculation of the confidence band.
#'
#' @export

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
