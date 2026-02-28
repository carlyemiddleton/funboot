#' @title Transforms a spatial variable in long format and into wide format for input into pffr()
#'
#' @param data Data frame obtained from the output of *funboot::preprocess_data()*
#' @param spatial_variable The spatial variable to perform the transformation on.  For example, 'outcome' or 'L_obs'
#' @param grid The grid on which the spatial summary function is evaluated.  For example, 0:200
#' @param id The 'by' variable for the transformation, typically the image ID number.  For example, 'image_id'
#'
#' @return A wide-format object compatible with the **data** argument in *pffr()*

format_spatial_variable <- function(data, spatial_variable, grid, id){
  column <- function(i){data[data[id] == i,][spatial_variable][1] } #calling the outcome function Y
  mat_transpose <- data.frame(grid = grid)
  for(i in sort(unique(data[id])[[1]])){
    mat_transpose <- cbind(mat_transpose, column(i))
  }
  mat <- t(mat_transpose[,-c(1)])
  return(mat)
}




#' @title Fits a *pffr()* model to spatial summary function data
#'
#' @param formula An object of class *formula*.  The formula for the *pffr()* model with special terms as in *mgcv's* "bam"
#' @param data Data frame obtained from the output of *funboot::preprocess_data()*
#' @param spatial_covars Vector containing the names of the covariate(s) to be treated as spatially-varying.  For example, c('var1','var2')
#' @param nthreads Number of threads to use in the model fit
#'
#' @return A *pffr()* model object
#'
#' @export

fit_model <- function(formula, data, spatial_covars = NULL, nthreads = 1){

  n <- dim(unique(data['patient_id']))[1]
  n_Im <- dim(unique(data['image_number']))[1]
  grid <- sort(unique(data$r))

  #format image-level covariates
  image_level_covars <- unique(data.frame(data[['image_number']], data[all.vars(formula)[c(-1,-which(all.vars(formula)%in%spatial_covars))]]))
  names(image_level_covars)[1] <- 'image_number'
  pffr_data <- image_level_covars
  if(dim(pffr_data)[1]!=n_Im){stop('spatial_covars argument may be invalid')}

  #format outcome function
  outcome <- format_spatial_variable(data=data, spatial_variable = 'outcome', grid=grid, id='image_number')
  pffr_data$outcome <- outcome

  #format spatial covariates
  if(!is.null(spatial_covars)){
    for(i in 1:length(spatial_covars)){
      pffr_data$temp <- format_spatial_variable(data=data, spatial_variable = spatial_covars[i], grid=grid, id='image_number')
      names(pffr_data)[names(pffr_data) == 'temp'] <- paste0(spatial_covars[i])
    }
  }

  pffr_data <- as.data.frame(pffr_data)
  #Fit the functional regression model
  pffrmodel <- refund::pffr(
    formula   = formula,
    data      = pffr_data,
    yind      = grid,
    algorithm = "bam",
    discrete  = TRUE,
    nthreads  = nthreads
  )
  return(pffrmodel)

}



#' @title Extracts the target curve (e.g. the curve to create the wild bootstrap confidence band around) from a *pffr()* model object
#'
#' @param model A *pffr()* model object
#' @param grid The grid on which the spatial summary function is evaluated.  For example, 0:200
#' @param target A vector containing the names of variables and/or coefficients used to calculate the target curve.  For example, c("var1","var2").  If the target curve includes the intercept curve, specify it as "Intercept".
#' @param form A vector of equal length as *target*.  Each element should be be either "coef" or "term", depending on whether the the corresponding element of *target* represents a coefficient curve (e.g. \eqn{\beta_1(r)}) or full term (e.g. \eqn{\beta_1(r)X_1(r)}).
#' @param covar_list a list containing the values of the covariates across radii.
#'                 For example, if we want a confidence band for \eqn{\beta_0(r) + \beta_1(r)Group + \beta_2(r)L_{AB}(r)} for patient 1, use **covar.list** = list(group = vec1, Lab = vec2), where vec1 and vec2 are R-length vectors containing covariate values for patient 1 at each radii, in ascending order of radii
#'                 The names of *covar.list* must match the names of *target*.
#'                 To calculate a confidence band for \eqn{\beta_0(r)} only, use **covar.list**=list().
#' @param func The function of *form* which produces the target curve.  Its arguments should match *target*.
#'
#' @return A vector containing the target curve.
#'
#' @export

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


#' @title Extracts package SE's
#'
#' @param name A *pffr()* model object
#' @param model A *pffr()* model object
#' @param grid The grid on which the spatial summary function is evaluated.  For example, 0:200
#'
#' @return A vector containing the target curve.
#'
#' @export

extract_package_SEs <- function(name, model, grid) {
  sm <- stats::coef(model, n1 = length(grid))$smterms
  if (name %in% names(sm)) result <- sm[[name]]$coef$se
  ind <- grep(name, names(sm), value = TRUE)
  if (length(ind) == 1) result <- sm[[ind]]$coef$se
  if (length(ind) > 1) stop("Multiple smooth terms match for ", name)
  if (length(ind) == 0)stop("No smooth term found matching ", name)
  return(result)
}



#' @title Fits a *pffr()* model to spatial summary function data
#'
#' @param model An object of class *formula*.  The formula for the *pffr()* model with special terms as in *mgcv's* BAM
#' @param n_Im_total Data frame obtained from the output of *funboot::preprocess_data()*
#' @param grid Vector containing the names of the covariate(s) to be treated as spatially-varying.  For example, c('var1','var2')
#' @param id Vector containing the names of the covariate(s) to be treated as spatially-varying.  For example, c('var1','var2')
#'
#' @return A *pffr()* model object
#'
#' @export

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

