#' @title Preprocesses data
#'
#' @param sumfun.data is
#' @return It returns res.wildBS
#'
#' @export

pointwise_test <- function(data, formula = formula, r.star = 100, re=c('patient_id'),alpha=.05,
                           spatial.covars = NULL,method = c('Bonferroni')){

  n <- dim(unique(data['patient_id']))[1]
  n.Im <- dim(unique(data['image_number']))[1]
  grid <- sort(unique(data$r))
  if(sum(r.star %in% grid)!=length(r.star)){stop('each r.star must be within the grid of radii used to compute the spatial summary function')}

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

  ##Perform the pointwise test
  #Estimates
  beta0_hat_constant <- rep(summary(pffrmodel)$p.coeff, length(grid))
  invisible(capture.output(beta0_hat <- coef(pffrmodel,n1=length(grid))$smterms$`Intercept(grid)`$coef$value))
  beta_hat <- data.frame(beta0_hat_constant = beta0_hat_constant,
                         beta0_hat = beta0_hat)
  if(length(attr(terms(formula), "term.labels")) != 0){
    for(i in names(coef(pffrmodel,n1=length(grid))$smterms[-1])){
      invisible(capture.output(list <- coef(pffrmodel,n1=length(grid))$smterms[
        which(names(coef(pffrmodel,n1=length(grid))$smterms) == i )]))
      assign(paste0('beta_hat_',i),
             list[[1]]$coef$value)
    }
    beta_hat <- cbind(beta_hat, bind_cols(mget(paste0('beta_hat_',
                names(coef(pffrmodel,n1=length(grid))$smterms)[!grepl(re, names(coef(pffrmodel,n1=length(grid))$smterms)) &
                                                               !grepl("Intercept", names(coef(pffrmodel,n1=length(grid))$smterms)) ]  ))))
  }
  #SEs
  beta0_hat_constant.se <- rep(summary(pffrmodel)$p.table[2], length(grid))
  invisible(capture.output(beta0_hat.se <- coef(pffrmodel,n1=length(grid))$smterms$`Intercept(grid)`$coef$se))
  beta_hat.se <- data.frame(beta0_hat_constant = beta0_hat_constant.se,
                         beta0_hat = beta0_hat.se)
  if(length(attr(terms(formula), "term.labels")) != 0){
    for(i in names(coef(pffrmodel,n1=length(grid))$smterms[-1])){
      invisible(capture.output(list <- coef(pffrmodel,n1=length(grid))$smterms[
        which(names(coef(pffrmodel,n1=length(grid))$smterms) == i )]))
      assign(paste0('beta_hat_',i),
             list[[1]]$coef$se)
    }
    beta_hat.se <- cbind(beta_hat.se, bind_cols(mget(paste0('beta_hat_',
                                                      names(coef(pffrmodel,n1=length(grid))$smterms)[!grepl(re, names(coef(pffrmodel,n1=length(grid))$smterms)) &
                                                                                                       !grepl("Intercept", names(coef(pffrmodel,n1=length(grid))$smterms)) ]  ))))
  }
  #Adjusted p values
  Ests <- beta_hat[grid %in% r.star,]
  SEs <- beta_hat.se[grid %in% r.star,]
  Zs <- Ests/SEs
  p.unadjusted <- sapply(abs(Zs), function(x){2*pnorm(x, lower.tail=F)})
  if(length(r.star)>1){ p.adjusted <- apply(p.unadjusted, 2, p.adjust, method=paste0(method))}
  if(length(r.star)==1){ p.adjusted <- p.unadjusted}
  p.adjusted <- cbind(r.star, p.adjusted)
  if(method=='bonferroni'){
    LB <- cbind(r.star, Ests - qnorm(alpha/2/length(r.star), lower.tail=F)*SEs)
    UB <- cbind(r.star, Ests + qnorm(alpha/2/length(r.star), lower.tail=F)*SEs)
    return(list(p.adjusted=p.adjusted, LB = LB, UB = UB))
  }else if (method == 'none'){
    LB <- cbind(r.star, Ests - qnorm(alpha/2, lower.tail=F)*SEs)
    UB <- cbind(r.star, Ests + qnorm(alpha/2, lower.tail=F)*SEs)
    return(list(p.adjusted=p.adjusted, LB = LB, UB = UB))
  }else{
    return(list(p.adjusted=p.adjusted))
  }

}
