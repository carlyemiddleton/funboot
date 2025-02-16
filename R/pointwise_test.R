#' @title Preprocesses data
#'
#' @param sumfun.data is
#' @return It returns res.wildBS
#'
#' @export

pointwise_test <- function(data=sumfun.data, formula = formula, r.star = 100, re=c('patient_id'), patient.id=c('patient_id'),alpha=.05,
                           image.id=c('image_number'),spatial.covars = c('g.obs'),method = c('Bonferroni')){

  n <- dim(unique(data[paste0(patient.id)]))[1]
  n.Im <- dim(unique(data[paste0(image.id)]))[1]
  grid <- sort(unique(data$r))
  if(sum(r.star %in% grid)!=length(r.star)){stop('each r.star must be within the grid of radii used to compute the spatial summary function')}

  #format image-level covariates
  image.level.covars <- unique(data.frame(data[[paste0(image.id)]], data[all.vars(formula)[c(-1,-which(all.vars(formula)%in%spatial.covars))]]))
  names(image.level.covars)[1] <- paste0(image.id)
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
    for(i in 1:length(attr(terms(formula), "term.labels"))){
      invisible(capture.output(list <- coef(pffrmodel,n1=length(grid))$smterms[
        #which(names(coef(pffrmodel,n1=length(grid))$smterms) %in% attr(terms(formula), "term.labels")[i])]
        which(sapply(names(coef(pffrmodel,n1=length(grid))$smterms), function(x){grepl(paste0(all.vars(formula)[i+1]), x)}) )]))
      assign(paste0('beta_hat_',all.vars(formula)[i+1]),
             list[[1]]$coef$value)
    }
    beta_hat <- cbind(beta_hat, bind_cols(mget(paste0('beta_hat_',
                                                      all.vars(formula)[!all.vars(formula) %in% re][-1]  ))))
  }
  #SEs
  beta0_hat_constant.se <- rep(summary(pffrmodel)$p.table[2], length(grid))
  invisible(capture.output(beta0_hat.se <- coef(pffrmodel,n1=length(grid))$smterms$`Intercept(grid)`$coef$se))
  beta_hat.se <- data.frame(beta0_hat_constant = beta0_hat_constant.se,
                            beta0_hat = beta0_hat.se)
  if(length(attr(terms(formula), "term.labels")) != 0){
    for(i in 1:length(attr(terms(formula), "term.labels"))){
      invisible(capture.output(list <- coef(pffrmodel,n1=length(grid))$smterms[
        which(sapply(names(coef(pffrmodel,n1=length(grid))$smterms), function(x){grepl(paste0(all.vars(formula)[i+1]), x)}) )]))
      assign(paste0('beta_hat_',all.vars(formula)[i+1]),
             list[[1]]$coef$se)
    }
    beta_hat.se <- cbind(beta_hat.se, bind_cols(mget(paste0('beta_hat_',
                                                            all.vars(formula)[!all.vars(formula) %in% re][-1]  ))))
  }
  #Adjusted p values
  Ests <- beta_hat[grid %in% r.star,]
  SEs <- beta_hat.se[grid %in% r.star,]
  Zs <- Ests/SEs
  p.unadjusted <- sapply(Zs, function(x){pnorm(x, alpha/2, lower.tail=T)})
  if(length(r.star)>1){ p.adjusted <- apply(p.unadjusted, 2, p.adjust, method=paste0(method))}
  if(length(r.star)==1){ p.adjusted <- p.unadjusted}
  p.adjusted <- cbind(r.star, p.adjusted)
  return(list(p.adjusted=p.adjusted))

}
