#' @title Calculates predicted values for linear combinations of *pffr()* model coefficient functions
#'
#' @param beta_hat An R x p matrix of values containing the estimated model coefficient functions from *pffr()*
#' @param lin.comb a list containing the values of the covariates across radii in the linear combination.
#'                 For example, if we want a confidence band for the linear combination \eqn{\beta_0(r) + \beta_1(r)Group + \beta_2(r)L_{AB}(r)} for patient 1, use **lin.comb** = list(group = data[data$patient_id=='1',]$group, Lab = data[data$patient_id=='1',]$Lab).
#'                 To calculate a confidence band for \eqn{\beta_0(r)} only, use **lin.comb**=list()
#' @return A length-R vector containing predicted values for the linear combination across radii

calc.preds <- function(beta_hat, lin.comb){
  newdata <- beta_hat
  newdata[,1:2] <- 1
  if(length(lin.comb)>0){
    for(i in 1:length(lin.comb)){newdata[,i+2] <- lin.comb[[i]]}
  }
  newdata <- newdata[,1:(length(lin.comb)+2)]
  preds <- as.matrix(newdata)*as.matrix(beta_hat[,1:(length(lin.comb)+2)])
  preds <- apply(preds, 1, sum)
  return(preds)
}


#' @title Calculates wild bootstrap confidence bands for linear combinations of *pffr()* model coefficient functions
#'
#' @param formula An object of class *formula*.  The formula for the *pffr()* model
#' @param lin.comb a list containing the values of the covariates across radii in the linear combination.
#'                 For example, if we want a confidence band for the linear combination \eqn{\beta_0(r) + \beta_1(r)Group + \beta_2(r)L_{AB}(r)} for patient 1, use **lin.comb** = list(group = data[data$patient_id=='1',]$group, Lab = data[data$patient_id=='1',]$Lab).
#'                 To calculate a confidence band for \eqn{\beta_0(r)} only, use **lin.comb**=list()
#' @param data Data frame obtained from the output of *funboot::preprocess_data()*
#' @param spatial.covars Vector containing the names of the covariate(s) to be specified as spatially-varying.  For example, c('var1','var2')
#' @param B1 Number of bootstrap samples in the outside layer
#' @param B2 Number of bootstrap samples in the inside layer
#' @param alpha Desired significance level for the confidence band. For example, 0.05
#' @param re Vector containing names of variable(s) for which to fit a random intercept.  For example, c('patient_id')
#' @param n_cores Number of cores to use during parallel processing.  If NULL, parallel::detectCores(logical=FALSE)-1 cores are used
#' @param seed Random seed to be used during the bootstrap resampling (optional).
#'
#' @return A list containing the lower and upper bounds of the confidence band, as well as estimated model coefficients from the *pffr()* output and
#'         additional intermittent variables used in the calculation of the confidence band.  'p.vals' are p-values for each model coefficient function yielded by the wild boostrap test of the model coefficient against 0, and 'test.stats' are the associated test statistics.
#'         'M1' is the bootstrap distribution of M from the outer layer.  'bs.sd' are the bootstrap standard deviations for each radius, yielded by the inner layer of resampling
#' @export

lin_comb_CB <- function(formula,lin.comb=list(), data, spatial.covars = NULL, B1=500, B2=100,
                        alpha=.05,re=NULL,n_cores=NULL,seed=NULL){
  first_layer <- function(b1, B2=B2,preds=preds,pffrmodel=pffrmodel,c=c,e=function(i){outcome[i,] - predict(pffrmodel)[i,]},
                          n.Im=100,grid=grid,model.bs=model.bs,c.bs=c.bs,e.bs=e.bs,pffr.data.bs=pffr.data.bs,pffr.data=pffr.data,
                          formula=formula,re=re,outcome=outcome,lin.comb=lin.comb,calc.preds=calc.preds,PREDS.BS.BS=PREDS.BS.BS, seed=seed
  ){if(!is.null(seed)){set.seed(seed)}
    Y.bs <- matrix(NA, nrow=n.Im, ncol=length(grid))
    for(i in 1:n.Im){
      Y.bs[i,] <- predict(pffrmodel)[i,] + c[i,b1]*e(i)
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
      if(is.null(re)){
        beta_hat.bs <- cbind(beta_hat.bs, do.call('cbind', mget(paste0('beta_hat_',
                                                                       names(coef(model.bs,n1=length(grid))$smterms)[!grepl("Intercept", names(coef(model.bs,n1=length(grid))$smterms)) ]  ))))
      }else{
        beta_hat.bs <- cbind(beta_hat.bs, do.call('cbind', mget(paste0('beta_hat_',
                                                                       names(coef(model.bs,n1=length(grid))$smterms)[!grepl(re, names(coef(model.bs,n1=length(grid))$smterms)) &
                                                                                                                       !grepl("Intercept", names(coef(model.bs,n1=length(grid))$smterms)) ]  ))))
      }
    }
    #bootstrap layer 2:  estimate the bootstrap standard error of lin.comb
    e.bs <- function(i){outcome[i,] - predict(model.bs)[i,]}
    ##Step B2.1
    c.bs <- matrix(NA, nrow=n.Im, ncol=B2)
    for(b in 1:B2){
      for(i in 1:n.Im){
        c.bs[i,b] <- ifelse(rbinom(1, size=1, prob=.5) == 1, 1, -1)
      }
    }
    M2 <- matrix(NA, ncol=1, nrow=B2)
    preds.bs <- calc.preds(beta_hat = beta_hat.bs, lin.comb = lin.comb)
    PREDS.BS[,b1] <- preds.bs
    PREDS.BS.BS <- matrix(NA, ncol=B2, nrow=length(preds))
    for(b2 in 1:B2){
      Y.bs.bs <- matrix(NA, nrow=n.Im, ncol=length(grid))
      for(i in 1:n.Im){
        Y.bs.bs[i,] <- predict(model.bs)[i,] + c.bs[i,b2]*e.bs(i)
      }
      data.bs.bs <- data.frame(grid = grid)
      for(i in 1:n.Im){
        data.bs.bs <- cbind(data.bs.bs, Y.bs.bs[i,])
      }
      Y.mat.bs.bs <- t(data.bs.bs[,-c(1)])
      pffr.data.bs.bs <- pffr.data.bs
      pffr.data.bs.bs$Y.mat.bs.bs <- Y.mat.bs.bs
      formula.bs.bs <- update.formula(formula, as.formula('Y.mat.bs.bs ~ .'))
      model.bs.bs <- refund:::pffr(formula=formula.bs.bs , data=pffr.data.bs.bs, yind=grid)
      #retrieve the coefficient estimates from the pffr() output
      beta0_hat_constant.bs.bs <- rep(summary(model.bs.bs)$p.coeff, length(grid))
      invisible(capture.output(beta0_hat.bs.bs <- coef(model.bs.bs,n1=length(grid))$smterms$`Intercept(grid)`$coef$value))
      beta_hat.bs.bs <- data.frame(beta0_hat_constant = beta0_hat_constant.bs.bs,
                                   beta0_hat = beta0_hat.bs.bs)
      if(length(attr(terms(formula), "term.labels")) != 0){
        for(i in names(coef(model.bs.bs,n1=length(grid))$smterms[-1])){
          invisible(capture.output(list <- coef(model.bs.bs,n1=length(grid))$smterms[
            which(names(coef(model.bs.bs,n1=length(grid))$smterms) == i )]))
          assign(paste0('beta_hat_',i),
                 list[[1]]$coef$value)
        }
        if(is.null(re)){
          beta_hat.bs.bs <- cbind(beta_hat.bs.bs, do.call('cbind', mget(paste0('beta_hat_',
                                                                               names(coef(model.bs.bs,n1=length(grid))$smterms)[!grepl("Intercept", names(coef(model.bs.bs,n1=length(grid))$smterms)) ]  ))))
        }else{
          beta_hat.bs.bs <- cbind(beta_hat.bs.bs, do.call('cbind', mget(paste0('beta_hat_',
                                                                               names(coef(model.bs.bs,n1=length(grid))$smterms)[!grepl(re, names(coef(model.bs.bs,n1=length(grid))$smterms)) &
                                                                                                                                  !grepl("Intercept", names(coef(model.bs.bs,n1=length(grid))$smterms)) ]  ))))
        }
      }
      PREDS.BS.BS[,b2] <- calc.preds(beta_hat = beta_hat.bs.bs, lin.comb = lin.comb)
      print(paste0('completed second layer bootstrap sample ',b2))
    }
    load(file='M1.RData')
    load(file='PREDS.BS.RData')
    m1 <- max(abs(preds.bs - preds)/apply(PREDS.BS.BS, 1, sd))
    M1[,b1] <- m1
    PREDS.BS[,b1] <- preds.bs
    save(M1, file='M1.RData')
    save(PREDS.BS, file='PREDS.BS.RData')
    print(paste0('completed first layer bootstrap sample ',b1))
    return(list(m1=m1, preds.bs=preds.bs))
  }
  #########################################
  if(!is.null(seed)){set.seed(seed)}
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


  ###########################
  ## Wild bootstrap layer 1 #
  ###########################

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
    if(is.null(re)){
      beta_hat <- cbind(beta_hat, do.call('cbind', mget(paste0('beta_hat_',
                                                               names(coef(pffrmodel,n1=length(grid))$smterms)[!grepl("Intercept", names(coef(pffrmodel,n1=length(grid))$smterms)) ]  ))  ))
    }else{
      beta_hat <- cbind(beta_hat, do.call('cbind', mget(paste0('beta_hat_',
                                                               names(coef(pffrmodel,n1=length(grid))$smterms)[!grepl(re, names(coef(pffrmodel,n1=length(grid))$smterms)) &
                                                                                                                !grepl("Intercept", names(coef(pffrmodel,n1=length(grid))$smterms)) ]  ))))
    }

  }
  e <- function(i){outcome[i,] - predict(pffrmodel)[i,]}
  ##Step B2.1
  c <- matrix(NA, nrow=n.Im, ncol=B1)
  for(b in 1:B1){
    for(i in 1:n.Im){
      c[i,b] <- ifelse(rbinom(1, size=1, prob=.5) == 1, 1, -1)
    }
  }
  ##Step B2.2
  M1 <- matrix(NA, nrow=1, ncol=B1)
  #calculate preds (the value we want to bootstrap)
  preds <- calc.preds(beta_hat = beta_hat, lin.comb = lin.comb)
  PREDS.BS <- matrix(NA, ncol=B1, nrow=length(preds))
  save(M1, file='M1.RData')
  save(PREDS.BS, file='PREDS.BS.RData')

  if(is.null(n_cores)){n_cores <- detectCores(logical=FALSE)-1}
  cl <- makeCluster(n_cores, outfile='cluster-output.txt')
  first.layer.out <- parLapply(cl=cl, as.list(1:B1), first_layer, B2=B2, preds=preds, pffrmodel=pffrmodel,c=c,pffr.data=pffr.data,
                               n.Im=n.Im,grid=grid,
                               formula=formula,re=re,outcome=outcome,lin.comb=lin.comb,calc.preds=calc.preds, seed=456)
  M1 <- do.call('cbind', do.call('cbind',first.layer.out)[1,])
  PREDS.BS <- do.call('cbind', do.call('cbind',first.layer.out)[2,])
  save(M1, file='M1.RData')
  save(PREDS.BS, file='PREDS.BS.RData')

  #get q and calculate the confidence bands
  q <- quantile(M1, 1-alpha, na.rm=T)
  MoE <- q*apply(PREDS.BS, 1, sd)
  CB.lower <- preds - MoE
  CB.upper <- preds + MoE
  test.stat <- max(abs(preds - 0)/apply(PREDS.BS, 1, sd))
  p.val <- mean(M1>test.stat)
  CBs <- list(CB.lower, CB.upper, grid, preds, apply(PREDS.BS, 1, sd), q, re, M1, test.stat, p.val)
  names(CBs) <- c('CB.lower','CB.upper', 'grid', 'estimate', 'bs.sd', 'q', 're', 'M1', 'test.stat', 'p.val')
  return(CBs)
}
