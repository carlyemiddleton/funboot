#' @title format.spatial.variable takes a spatial variable in long format and puts it into matrix format for input into pffr()
#' @param sumfun.data is
#' @param covars is a vector of model covariate names found in sumfun.data
#' @return It returns res.wildBS
#' @export

format.spatial.variable <- function(data, spatial.variable, grid, id){
    column <- function(i){data[data[id] == i,][spatial.variable][1] } #calling the outcome function Y
    mat.transpose <- data.frame(grid = grid)
    for(i in sort(unique(data[id])[[1]])){
      mat.transpose <- cbind(mat.transpose, column(i))
    }
    mat <- t(mat.transpose[,-c(1)])
    return(mat)
}





#' @title Preprocesses data
#' @param sumfun.data is
#' @param covars is a vector of model covariate names found in sumfun.data
#' @return It returns res.wildBS
#' @export

#formula <- outcome ~ patient_age + tumor_grade + L.obs #+ s(image_number, bs="re")
#data <- sumfun.data
#data$tumor_grade <- ifelse(data$tumor_grade=='1', 1, 0) #omit this to see if our package returns pffr() errors
#spatial.covars = c('L.obs')

fit_model <- function(formula, data, spatial.covars = NULL){

  #format outcome function
  n <- length(unique(data$patient_id))
  grid <- sort(unique(data$r))
  outcome <- format.spatial.variable(data=data, spatial.variable = 'outcome', grid=grid, id='image_number')
  #format spatial covariates
  #Note:  make an error message checking if spatial.covars match the formula, which matches the dataset)
  if(!is.null(spatial.covars)){
    for(i in 1:length(spatial.covars)){
      assign(paste0(spatial.covars[i]),
            format.spatial.variable(data=data, spatial.variable = spatial.covars[i], grid=grid, id='image_number') )
    }
  }
  #format image-level covariates
  image.level.covars <- unique(data.frame(data['image_number'], data[all.vars(formula)[c(-1,-which(all.vars(formula)%in%spatial.covars))]]))


  #Fit the functional regression model
  pffrmodel <- refund::pffr(formula, data=image.level.covars, yind=grid)
  return(pffrmodel)

}








  #put in the documentation:  this package uses one constant and one functional intercept, and does not allow suppression of either
  # beta0_constant <- rep(summary(pffrmodel)$p.coeff, length(grid))
  # beta0_hat <- coef(pffrmodel,n1=length(grid))$smterms$`Intercept(grid)`$coef$value
  # for(i in 1:length(formula)){
  #   list <- coef(pffrmodel,n1=length(grid))$smterms[
  #     which(gsub('.{6}$', '', names(coef(pffrmodel,n1=length(grid))$smterms))==all.vars(formula)[i+1])]
  #   assign(paste0('beta_hat_',all.vars(formula)[i+1]),
  #          list[[1]]$coef$value)
  # }




























