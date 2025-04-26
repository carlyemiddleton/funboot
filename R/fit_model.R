#' @title Takes a spatial variable in long format and puts it into wide format for input into pffr()
#'
#' @param data Data frame obtained from the output of *funboot::preprocess_data()*
#' @param spatial.variable The spatial variable to perform the transformation on.  For example, 'outcome' or 'L.obs'
#' @param grid The grid on which the spatial summary function is evaluated.  For example, 0:200
#' @param id The 'by' variable for the transformation, typically the image ID number.  For example, 'image_id'
#'
#' @return A wide-format object compatible with the **data** argument in *pffr()*

format.spatial.variable <- function(data, spatial.variable, grid, id){
  column <- function(i){data[data[id] == i,][spatial.variable][1] } #calling the outcome function Y
  mat.transpose <- data.frame(grid = grid)
  for(i in sort(unique(data[id])[[1]])){
    mat.transpose <- cbind(mat.transpose, column(i))
  }
  mat <- t(mat.transpose[,-c(1)])
  return(mat)
}




#' @title Fits a the *pffr()* model to spatial summary function data
#'
#' @param formula An object of class *formula*.  The formula for the *pffr()* model
#' @param data Data frame obtained from the output of *funboot::preprocess_data()*
#' @param spatial.covars Vector containing the names of the covariate(s) to be specified as spatially-varying.  For example, c('var1','var2')
#'
#' @return A *pffr()* model object.  See *refund::pffr()*
#'
#' @export

fit_model <- function(formula, data, spatial.covars = NULL){

  n <- dim(unique(data['patient_id']))[1]
  n.Im <- dim(unique(data['image_number']))[1]
  grid <- sort(unique(data$r))

  #format image-level covariates
  image.level.covars <- unique(data.frame(data[['image_number']], data[all.vars(formula)[c(-1,-which(all.vars(formula)%in%spatial.covars))]]))
  names(image.level.covars)[1] <- 'image_number'
  pffr.data <- image.level.covars
  if(dim(pffr.data)[1]!=n.Im){stop('spatial.covars argument may be invalid')}

  #format outcome function
  outcome <- format.spatial.variable(data=data, spatial.variable = 'outcome', grid=grid, id='image_number')
  pffr.data$outcome <- outcome

  #format spatial covariates
  if(!is.null(spatial.covars)){
    for(i in 1:length(spatial.covars)){
      pffr.data$temp <- format.spatial.variable(data=data, spatial.variable = spatial.covars[i], grid=grid, id='image_number')
      names(pffr.data)[names(pffr.data) == 'temp'] <- paste0(spatial.covars[i])
    }
  }

  pffr.data <- as.data.frame(pffr.data)
  #Fit the functional regression model
  pffrmodel <- refund:::pffr(formula=formula, data=pffr.data, yind=grid)
  return(pffrmodel)

}




































