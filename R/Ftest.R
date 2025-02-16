#' @title Preprocesses data
#'
#' @param grid the grid that the CBs were evaluated over.  eg. 0:200
#' @return It returns res.wildBS
#'
#' @export

Ftest <- function(formula.full, formula.red, image.id='image_number',patient.id='patient_id',
                  data, spatial.covars = NULL, B=1000,alpha=.05,re=NULL,seed=456){

  n <- dim(unique(data[paste0(patient.id)]))[1]
  n.Im <- dim(unique(data[paste0(image.id)]))[1]
  grid <- sort(unique(data$r))

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
  ##calculate the observed F
  model.full <- refund:::pffr(formula=formula.full, data=pffr.data, yind=grid)
  model.red <- refund:::pffr(formula=formula.red, data=pffr.data, yind=grid)
  rss.full <- apply((pffr.data$outcome - predict(model.full))^2, 2, sum)
  rss.red <- apply((pffr.data$outcome - predict(model.red))^2, 2, sum)
  (F.test <- sum((rss.red-rss.full)/(length(attr(terms(formula.full), "term.labels"))-
                                       length(attr(terms(formula.red), "term.labels")))/
                   (rss.full/(n.Im-length(attr(terms(formula.full), "term.labels")))) ) )

  ##########################
  ## Wild bootstrap F test #
  ##########################
  e <- function(i){pffr.data$outcome[i,] - predict(model.red)[i,] }
  ##Step B2.1
  c <- matrix(NA, nrow=n.Im, ncol=B)
  for(b in 1:B){
    for(i in 1:n.Im){
      c[i,b] <- ifelse(rbinom(1, size=1, prob=.5) == 1, 1, -1)
    }
  }
  ##Step B2.2
  M <- rep(NA, B)
  for(b in 1:B){
    Y.bs <- matrix(NA, nrow=n.Im, ncol=length(grid))
    for(i in 1:n.Im){
      Y.bs[i,] <- predict(model.red)[i,] + c[i,b]*e(i) #resample from the null assumption model
    }
    data.bs <- data.frame(grid = grid)
    for(i in 1:n.Im){
      data.bs <- cbind(data.bs, Y.bs[i,])
    }
    Y.mat.bs <- t(data.bs[,-c(1)])
    pffr.data.bs <- pffr.data
    pffr.data.bs$Y.mat.bs <- Y.mat.bs
    formula.full.bs <- update.formula(formula.full, as.formula('Y.mat.bs ~ .'))
    model.full.bs <- refund:::pffr(formula=formula.full.bs , data=pffr.data.bs, yind=grid)
    formula.red.bs <- update.formula(formula.red, as.formula('Y.mat.bs ~ .'))
    model.red.bs <- refund:::pffr(formula=formula.red.bs , data=pffr.data.bs, yind=grid)
    rss.full.bs <- apply((Y.mat.bs - predict(model.full.bs))^2, 2, sum)
    rss.red.bs <- apply((Y.mat.bs - predict(model.red.bs))^2, 2, sum)
    (F.b <- sum((rss.red.bs-rss.full.bs)/(length(attr(terms(formula.full), "term.labels"))-
                                           length(attr(terms(formula.red), "term.labels")))/
                  (rss.full.bs/(n.Im-length(attr(terms(formula.full), "term.labels")))) ) )
    M[b] <- F.b
    print(paste0('completed bootstrap sample ',b))
  }
  #hist(M)
  p.value <- mean(M > F.test, na.rm=T)
  return(p.value)

}

