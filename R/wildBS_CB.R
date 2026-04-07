#' @title Calculates wild bootstrap-based confidence bands for linear combinations of functional regression coefficients
#'
#' @param formula An object of class `formula`.  The formula passed to `refund::pffr()` with special terms as in `mgcv::gam()`
#' @param data A data frame obtained from the output of `funboot::preprocess_data()`
#' @param spatial_covars A character vector containing the names of the covariate(s) to be treated as spatially-varying.  For example, `c("var1","var2")`
#' @param target A character vector containing the names of the model terms to be combined to form the target curve for which to construct the confidence band around.  For example, `c("Intercept", "var1", "var2")`
#' @param form A character vector of equal length as `target`, with each entry corresponding to an entry in `target`.  Each entry must be either `"coef"` or `"term"`, where `"coef"` uses the estimated coefficient function directly and `"term"` multiplies the coefficient function by the corresponding value in `covar_list`.
#' @param covar_list A named list giving the covariate values used when constructing target-curve components with `form = "term"`. The names of `covar_list` must match `target`.  For example, if `target = c("Intercept", "var1", "var2")`, then `covar_list` could be `list(Intercept = rep(1, length(v1)), var1 = v1, var2 = v2)`, where `v1` and `v2` each have the length of the grid of evaluation radii.
#' @param func A function specifying how the terms in `target` are combined to form the target curve. The argument names of `func` must match `target`. For example, if `target = c("Intercept", "var1", "var2")`, then `func` might be `function(Intercept, var1, var2){Intercept + var1 + var2}`.
#' @param B1 Number of outer bootstrap samples
#' @param B2 Number of inner bootstrap samples
#' @param alpha Desired significance level for the confidence band. e.g., `0.05`.
#' @param re Name of variable for which to fit a random intercept.  For example, `"patient_id"`.
#' @param id Name of the patient ID variable. e.g., `"patient_id"`.
#' @param nthreads Number of parallel threads to use in the call to `refund::pffr()`
#'
#' @return A list containing, for each evaluation radius: the lower and upper bounds of the confidence band, estimated values of the target curve, and the standard deviation of outer bootstrap target curve estimates.  Also returns the grid of evaluation radii, the critical value, and the bootstrap distribution of the statistic `M`.
#'
#' @seealso [refund::pffr()], [mgcv::gam()]
#' @export

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
