#' @title Preprocesses data
#'
#' @param sumfun.data is
#' @return It returns res.wildBS
#'
#' @export

wildBS_CB <- function(formula, data, spatial.covars = NULL, B=1000,alpha=.05,re=NULL,seed=456){

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

  ###################
  ## Wild bootstrap #
  ###################

  ##Step B1
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
  e <- function(i){outcome[i,] - predict(pffrmodel)[i,]}
  ##Step B2.1
  c <- matrix(NA, nrow=n.Im, ncol=B)
  for(b in 1:B){
    for(i in 1:n.Im){
      c[i,b] <- ifelse(rbinom(1, size=1, prob=.5) == 1, 1, -1)
    }
  }
  ##Step B2.2
  M <- matrix(NA, ncol=c(sum(!grepl(re, names(coef(pffrmodel,n1=length(grid))$smterms)))+1), nrow=B)
  for(b in 1:B){
    Y.bs <- matrix(NA, nrow=n.Im, ncol=length(grid))
    for(i in 1:n.Im){
      Y.bs[i,] <- predict(pffrmodel)[i,] + c[i,b]*e(i)
    }
    data.bs <- data.frame(grid = grid)
    for(i in 1:n.Im){
      data.bs <- cbind(data.bs, Y.bs[i,])
    }
    Y.mat.bs <- t(data.bs[,-c(1)])
    pffr.data.bs <- pffr.data
    pffr.data.bs$Y.mat.bs <- Y.mat.bs
    formula.bs <- update.formula(formula, as.formula('Y.mat.bs ~ .'))
    model.bs <- refund:::pffr(formula=formula.bs , data=pffr.data.bs, yind=grid)
    #retrieve the coefficient estimates from the pffr() output
    beta0_hat_constant.bs <- rep(summary(model.bs)$p.coeff, length(grid))
    invisible(capture.output(beta0_hat.bs <- coef(model.bs,n1=length(grid))$smterms$`Intercept(grid)`$coef$value))
    beta_hat.bs <- data.frame(beta0_hat_constant = beta0_hat_constant.bs,
                           beta0_hat = beta0_hat.bs)
    if(length(attr(terms(formula), "term.labels")) != 0){
      for(i in names(coef(model.bs,n1=length(grid))$smterms[-1])){
        invisible(capture.output(list <- coef(model.bs,n1=length(grid))$smterms[
          which(names(coef(model.bs,n1=length(grid))$smterms) == i )]))
        assign(paste0('beta_hat_',i),
               list[[1]]$coef$value)
      }
      beta_hat.bs <- cbind(beta_hat.bs, bind_cols(mget(paste0('beta_hat_',
                                                        names(coef(model.bs,n1=length(grid))$smterms)[!grepl(re, names(coef(model.bs,n1=length(grid))$smterms)) &
                                                                                                         !grepl("Intercept", names(coef(model.bs,n1=length(grid))$smterms)) ]  ))))
    }
    #retrieve the SEs from the pffr() output
    beta0_hat_constant.se.bs <- rep(summary(model.bs)$p.table[2], length(grid))
    invisible(capture.output(beta0_hat.se.bs <- coef(model.bs,n1=length(grid))$smterms$`Intercept(grid)`$coef$se))
    beta_hat.se.bs <- data.frame(beta0_hat_constant = beta0_hat_constant.se.bs,
                              beta0_hat = beta0_hat.se.bs)
    if(length(attr(terms(formula), "term.labels")) != 0){
      for(i in names(coef(model.bs,n1=length(grid))$smterms[-1])){
        invisible(capture.output(list <- coef(model.bs,n1=length(grid))$smterms[
          which(names(coef(model.bs,n1=length(grid))$smterms) == i )]))
        assign(paste0('beta_hat_',i),
               list[[1]]$coef$se)
      }
      beta_hat.se.bs <- cbind(beta_hat.se.bs, bind_cols(mget(paste0('beta_hat_',
                                                              names(coef(model.bs,n1=length(grid))$smterms)[!grepl(re, names(coef(model.bs,n1=length(grid))$smterms)) &
                                                                                                               !grepl("Intercept", names(coef(model.bs,n1=length(grid))$smterms)) ]  ))))
    }
    M[b,] <- apply(abs((beta_hat.bs - beta_hat)/beta_hat.se.bs), 2, max)
    print(paste0('completed bootstrap sample ',b))
  }
  #apply(M, 2, hist)
  #retrieve the SEs from the pffr() output
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
  #get q and calculate the confidence bands
  q <- diag(apply(M, 2, function(x){quantile(x, 1-alpha, na.rm=T)}))
  MoE <- data.frame(as.matrix(beta_hat.se) %*% q)
  names(MoE) <- names(beta_hat)
  CB.lower <- beta_hat - MoE
  CB.upper <- beta_hat + MoE
  CBs <- list(CB.lower, CB.upper, grid, beta_hat, beta_hat.se, q, re)
  names(CBs) <- c('CB.lower','CB.upper', 'grid', 'pffr.Estimates', 'pffr.SEs', 'q', 're')
  return(CBs)
}
