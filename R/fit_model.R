#' @title format.spatial.variable takes a spatial variable in long format and puts it into matrix format for input into pffr()
#' @param sumfun.data is
#' @param covars is a vector of model covariate names found in sumfun.data
#' @return It returns res.wildBS

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




































