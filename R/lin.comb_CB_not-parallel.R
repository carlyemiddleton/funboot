calc.preds <- function(beta_hat, lin.comb){
  newdata <- beta_hat
  newdata[,1:2] <- 1
  for(i in 1:length(lin.comb)){
    newdata[,i+2] <- lin.comb[[i]]
  }
  newdata <- newdata[,1:(length(lin.comb)+2)]
  preds <- as.matrix(newdata)*as.matrix(beta_hat[,1:(length(lin.comb)+2)])
  preds <- apply(preds, 1, sum)
  return(preds)
}

lin.comb_CB <- function(formula,lin.comb=list(), data, spatial.covars = NULL, B1=500, B2=500,
                        alpha=.05,re=NULL,seed=NULL){
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
  M1 <- matrix(NA, ncol=1, nrow=B1)
  #calculate preds (the value we want to bootstrap)
  preds <- calc.preds(beta_hat = beta_hat, lin.comb = lin.comb)
  PREDS.BS <- matrix(NA, ncol=B1, nrow=length(preds))
  for(b1 in 1:B1){
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
      preds.bs.bs <- calc.preds(beta_hat = beta_hat.bs.bs, lin.comb = lin.comb)
      PREDS.BS.BS[,b2] <- preds.bs.bs
      print(paste0('completed second layer bootstrap sample ',b2))
    }
    M1[b1,] <- max(abs(preds.bs - preds)/apply(PREDS.BS.BS, 1, sd))
    print(paste0('completed first layer bootstrap sample ',b1))
  }
  
  #get q and calculate the confidence bands
  q <- quantile(M1, 1-alpha, na.rm=T)
  MoE <- q*apply(PREDS.BS, 1, sd)
  CB.lower <- preds - MoE
  CB.upper <- preds + MoE
  test.stat <- max(abs(preds - 0)/apply(PREDS.BS, 1, sd))
  p.val <- mean(M1>test.stat)
  CBs <- list(CB.lower, CB.upper, grid, preds, apply(PREDS.BS, 1, sd), q, re, M1, M2, test.stat, p.val)
  names(CBs) <- c('CB.lower','CB.upper', 'grid', 'estimate', 'bs.SE', 'q', 're', 'M1', 'M2', 'test.stat', 'p.val')
  return(CBs)
}


n_cores <- detectCores(logical=FALSE)
cl <- makeCluster(n_cores-1)
Sys.time()
parLapply(cl , as.list(seq_len(500)+seed), sim,
          counts = counts, nPatients = nPatients, nIm = nIm, window = window, sigma = sigma)
Sys.time()


formula=outcome ~ group + var.mat + var.group.mat;
lin.comb=list(group = 1); data=sumfun.data.full; spatial.covars = c('var.mat','var.group.mat');
B1=2;B2=2;alpha=.05;re=NULL;seed=NULL;
CBs <- lin.comb_CB(formula=outcome ~ group + var.mat + var.group.mat, lin.comb=list(group = 1), 
                   data=sumfun.data.full, spatial.covars = c('var.mat','var.group.mat'),
                   B1=2,B2=2,alpha=.05,re=NULL,seed=NULL,cores=NULL)