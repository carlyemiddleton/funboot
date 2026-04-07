#' @title Transforms a spatial variable in long format into wide format for input into `refund::pffr()`
#'
#' @param data A data frame obtained from the output of `funboot::preprocess_data()`
#' @param spatial_variable The spatial variable to perform the transformation on.  For example, `outcome`.
#' @param grid The grid of evaluation radii on which the spatial summary function is calculated.  For example, `0:200`.
#' @param id The grouping variable whose name defines the rows of the output matrix (typically the image ID number).  For example, `"image_number"`.
#'
#' @return A wide-format matrix compatible with the `data` argument in `refund::pffr()`, where rows represent images and columns represent radii.

format_spatial_variable <- function(data,
                                    spatial_variable,
                                    grid,
                                    id){
  column <- function(i){data[data[id] == i,][spatial_variable][1] }
  mat_transpose <- data.frame(grid = grid)
  for(i in sort(unique(data[id])[[1]])){
    mat_transpose <- cbind(mat_transpose, column(i))
  }
  mat <- t(mat_transpose[,-c(1)])
  return(mat)
}



#' @title Extracts a target curve from a `refund::pffr()` model object
#'
#' @param model A fitted model object obtained from the output of `refund::pffr()`
#' @param grid The grid of evaluation radii on which the spatial summary function is calculated.  For example, `0:200`
#' @param target A character vector containing the names of the model terms to be combined to form the target curve for which to construct the confidence band around.  For example, `c("Intercept", "var1", "var2")`
#' @param form A character vector of equal length as `target`, with each entry corresponding to an entry in `target`.  Each entry must be either `"coef"` or `"term"`, where `"coef"` uses the estimated coefficient function directly and `"term"` multiplies the coefficient function by the corresponding value in `covar_list`.
#' @param covar_list A named list giving the covariate values used when constructing target-curve components with `form = "term"`. The names of `covar_list` must match `target`.  For example, if `target = c("Intercept", "var1", "var2")`, then `covar_list` could be `list(Intercept = rep(1, length(v1)), var1 = v1, var2 = v2)`, where `v1` and `v2` each have the length of the grid of evaluation radii.
#' @param func A function specifying how the terms in `target` are combined to form the target curve. The argument names of `func` must match `target`. For example, if `target = c("Intercept", "var1", "var2")`, then `func` might be `function(Intercept, var1, var2){Intercept + var1 + var2}`.
#'
#' @return A vector containing estimated values of the target curve for each evaluation radius.
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



#' @title Extracts estimated fixed effects and random effects from a `refund::pffr()` model object
#'
#' @param model A fitted model object obtained from the output of `refund::pffr()`
#' @param n_Im_total The total number of observed images
#' @param grid The grid of evaluation radii on which the spatial summary function is calculated.  For example, `0:200`
#' @param id The grouping variable whose name defines the random effect in the fitted model.  For example, `"patient_id"`
#'
#' @return A list containing estimated fixed effects and random effects at each evaluation radius
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


