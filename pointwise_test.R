library(spatstat.data)
library(spatstat.univar)
library(spatstat.geom)
library(spatstat.random)
library(spatstat.explore)
library(refund)
data('data_example')
data <- data_example; rm(data_example)
summary.function <- 'L' #can also be 'K' or 'g'
from.cell <- 7
to.cell <- 3
R <- 200
inc <- 1
qc.cellcount.cutoff <- 20
P <- 50 #number of permutations to use
perm.yn <- F #whether or not to permute the outcome
image.dims=c(0,1000,0,1000) #xmin, xmax, ymin, ymax
image.xmax <- image.dims[2]
image.ymax <- image.dims[4]
image.xmin <- image.dims[1]
image.ymin <- image.dims[3]
W <- owin(c(image.xmin,image.xmax), c(image.ymin,image.ymax))

#########################
##Calculate model.data ## (it's Kdata, but for any summary function)
#########################
if(!(summary.function %in% c('K', 'L', 'g'))){stop("summary.function must be one of 'K', 'L', or g")}
Kdata <- K.pmean.vec <- NULL
if(summary.function=='L'){Ldata<- L.pmean.vec <- NULL}
if(summary.function=='g'){gdata<- g.pmean.vec <- NULL}
for(i in unique(data$image_number)){
  if(sum(data$image_number == i & data$cell_metacluster == from.cell)>qc.cellcount.cutoff & #quality control criteria
     sum(data$image_number == i & data$cell_metacluster == to.cell)>qc.cellcount.cutoff){
    qc.data <- data[data$image_number == i,]
    if(perm.yn==T){
      permuted.K <- NULL
      if(summary.function=='L'){permuted.L <- NULL}
      if(summary.function=='g'){permuted.g <- NULL}
      for(p in 1:P){
        ppp <- ppp(x = qc.data$cell_x, y = qc.data$cell_y,
                   marks = factor(sample(qc.data$cell_metacluster,
                                         size=length(qc.data$cell_metacluster), replace = F)), window = W)
        Kdata.temp <- data.frame(Kcross(ppp, i = from.cell, j = to.cell, correction='Ripley', r=seq(0,R,by=inc) ), image=i)
        permuted.K <- cbind(permuted.K, Kdata.temp$iso)
        if(summary.function=='L'){permuted.L <- sqrt(permuted.K/pi)}
        if(summary.function=='g'){
        gdata.temp <- data.frame(pcf(Kcross(ppp, i = from.cell, j = to.cell, correction='Ripley', r=seq(0,R,by=inc) ),
                                       method='c', divisor ="d"), image=i)
        permuted.g <- cbind(permuted.g, gdata.temp$pcf)
        }
      }
      K.pmean <- apply(permuted.K, 1, mean);  K.pmean.vec <- c(K.pmean.vec, K.pmean)
      if(summary.function=='L'){
      L.pmean <- apply(permuted.L, 1, mean); L.pmean.vec <- c(L.pmean.vec, L.pmean)
      }
      if(summary.function=='g'){
      g.pmean <- apply(permuted.g, 1, mean); g.pmean.vec <- c(g.pmean.vec, g.pmean)
      }
      print(paste0('permuted outcome for image ',i,' calculated'))
    }
    ppp <- ppp(x = qc.data$cell_x, y = qc.data$cell_y,
               marks = factor(qc.data$cell_metacluster), window = W)
    Kdata.temp <- data.frame(Kcross(ppp, i = from.cell, j = to.cell, correction='Ripley', r=seq(0,R,by=inc) ), image=i)
    names(Kdata.temp) <- c('r','K.expect','K.obs','image_number')
    Kdata <- rbind(Kdata, Kdata.temp)
    if(summary.function=='L'){
    Ldata <- Kdata
    Ldata$K.expect <- sqrt(Ldata$K.expect/pi)
    Ldata$K.obs <- sqrt(Ldata$K.obs/pi)
    names(Ldata) <- c('r','L.expect','L.obs','image_number')
    }
    if(summary.function=='g'){
    gdata.temp <- data.frame(pcf(Kcross(ppp, i = from.cell, j = to.cell, correction='Ripley', r=seq(0,R,by=inc) ),
                                   method='c', divisor ="d"), image=i)
    names(gdata.temp) <- c('r','g.expect','g.obs','image_number')
    gdata <- rbind(gdata, gdata.temp)
    }
  }
}
#Calculate the outcome variable
if(summary.function=='L'){
  model.data <- Ldata
  model.data$L.pmean <- ifelse(perm.yn==T, L.pmean.vec, NA)
  names(model.data) <- c('r','L.expect','L.obs','image_number','L.pmean.vec')
}else if(summary.function=='g'){
  model.data <- gdata
  model.data$g.pmean <- ifelse(perm.yn==T, g.pmean.vec, NA)
  names(model.data) <- c('r','g.expect','g.obs','image_number','g.pmean.vec')
}else{
  model.data <- Kdata
  model.data$K.pmean <- ifelse(perm.yn==T, K.pmean.vec, NA)
  names(model.data) <- c('r','K.expect','K.obs','image_number','K.pmean.vec')
}
covariate.df <- data
covariate.df$cell_id <- covariate.df$cell_x <- covariate.df$cell_y <- covariate.df$cell_metacluster <- NULL
covariate.df <- unique(covariate.df)
model.data <- merge(model.data, covariate.df, by='image_number'
                     , all=F) #all=F:  if an image doesn't meet the qc.cutoff, don't include its covariates

















########################################################################################################
##look at colocalization between metaclusters from=7 (Endothelial cells) and to=3 (T cells). We expect strong interaction
library(spatstat.data)
library(spatstat.univar)
library(spatstat.geom)
library(spatstat.random)
library(spatstat.explore) #go with 9,15 or 16
load(file='data.RData')
from.cell <- 7
to.cell <- 3
##Prepare the data for TVEM
##look at colocalization between metaclusters from=9 (Small circular cells) and
image.xmax <- 1000
image.ymax <- 1000
image.xmin <- 0
image.ymin <- 0
W <- owin(c(image.xmin,image.xmax), c(image.ymin,image.ymax))
Kdata <- NULL
qc.cellcount.cutoff <- 20
P <- 50 #number of permutations to use
for(i in unique(data$image_number)){
  if(sum(data$image_number == i & data$cell_metacluster == from.cell)>qc.cellcount.cutoff &
     sum(data$image_number == i & data$cell_metacluster == to.cell)>qc.cellcount.cutoff){
    qc.data <- data[data$image_number == i,]
    permuted.k <- NULL
    for(p in 1:P){
      ppp <- ppp(x = qc.data$cell_x, y = qc.data$cell_y,
                 marks = factor(sample(qc.data$cell_metacluster,
                                       size=length(qc.data$cell_metacluster), replace = F)), window = W)
      Kdata.temp <- data.frame(Kcross(ppp, i = from.cell, j = to.cell, correction='Ripley', r=seq(0,500,by=1) ), image=i)
      permuted.k <- cbind(permuted.k, Kdata.temp[,3])
    }
    k.pmean <- apply(permuted.k, 1, mean)
    ppp <- ppp(x = qc.data$cell_x, y = qc.data$cell_y,
               marks = factor(qc.data$cell_metacluster), window = W)
    Kdata.temp <- data.frame(Kcross(ppp, i = from.cell, j = to.cell, correction='Ripley', r=seq(0,500,by=1) ), image=i)
    Kdata.temp <- cbind(Kdata.temp, k.pmean)
    names(Kdata.temp) <- c('r','k.expect','k.obs','image_number','k.pmean')
    Kdata <- rbind(Kdata, Kdata.temp)
  }
}
Kdata$l.obs <- sqrt(Kdata$k.obs/pi)
Kdata$outcome <- Kdata$l.obs - Kdata$r
Kdata$l.pmean <- sqrt(Kdata$k.pmean/pi)
Kdata$poutcome <- Kdata$l.obs - Kdata$l.pmean
image.level.data <- data.frame(data$image_number, data$cells_per_image,
                               data$patient_age, data$patient_status, data$patient_DFS_months, data$patient_OS_months,
                               data$tumor_grade, data$tumor_clinical_type, data$tumor_size )
names(image.level.data) <- c('image_number', 'cells_per_image',
                             'patient_age', 'patient_status', 'patient_DFS_months', 'patient_OS_months',
                             'tumor_grade', 'tumor_clinical_type', 'tumor_size' )
Kdata <- merge(Kdata, unique(image.level.data), by='image_number', all=T)
Kdata <- Kdata[!is.na(Kdata$r),]



#############################################################
##calculate the confidence band and overlay it on the plot ##
#############################################################
##format data
n <- nPatients <- length(unique(Kdata$image_number))
nIm <- 1
s <- sort(unique(Kdata$r))
g <- length(s)
##Estimate beta_hat(s) using pffr()
Y <- function(i){Kdata[Kdata$image_number == i,]$poutcome[order(Kdata[Kdata$image_number == 1,]$r )]}
data <- data.frame(r = sort(unique(Kdata$r)))
for(i in unique(Kdata$image_number)){
  data <- cbind(data, Y(i))
}
Y.mat <- t(data[,-c(1)])
library(refund)
Y.mat <- Y.mat[!is.na(Y.mat[,1]),]
model <- pffr(Y.mat ~ 1, yind=sort(unique(Kdata$r)))
plot(coef(model, n1=g)$smterms$`Intercept(yindex)`$coef$yindex.vec,
     coef(model, n1=g)$smterms$`Intercept(yindex)`$coef$value
     ,ylim=c(-10,10))
lines(coef(model, n1=g)$smterms$`Intercept(yindex)`$coef$yindex.vec,
      coef(model, n1=g)$smterms$`Intercept(yindex)`$coef$value - qnorm(.025, lower.tail=F)*coef(model, n1=g)$smterms$`Intercept(yindex)`$coef$se
      , col='blue')
lines(coef(model, n1=g)$smterms$`Intercept(yindex)`$coef$yindex.vec,
      coef(model, n1=g)$smterms$`Intercept(yindex)`$coef$value + qnorm(.025, lower.tail=F)*coef(model, n1=g)$smterms$`Intercept(yindex)`$coef$se
      , col='blue')
abline(h=0, col='grey')




##Step B1
beta0_hat <- coef(model,n1=g)$smterms$`Intercept(yindex)`$coef$value
e <- Y.mat - predict(model)   #theresidual functions
s <- 0:500
n <- dim(Y.mat)[1]
##Step B2.1
B <- 1000
c <- matrix(NA, nrow=n, ncol=B)
for(b in 1:B){
  for(i in 1:n){
    c[i,b] <- ifelse(rbinom(1, size=1, prob=.5) == 1, 1, -1)
  }
}
##Step B2.2
beta_hat.bs <- matrix(NA, nrow=g, ncol=B)
M0 <-  rep(NA, B)
for(b in 1:B){
  Y.bs <- matrix(NA, nrow=n, ncol=g)
  for(i in 1:n){
    Y.bs[i,] <- predict(model)[i,] + c[i,b]*e[i,]
  }
  data.bs <- data.frame(s = s)
  for(i in 1:n){
    data.bs <- cbind(data.bs, Y.bs[i,])
  }
  Y.mat.bs <- t(data.bs[,-c(1)])
  model.bs <- pffr(Y.mat.bs ~  1, yind=sort(unique(Kdata$r)))
  beta0_hat.bs <- coef(model.bs,n1=g)$smterms$`Intercept(yindex)`$coef$value
  SE.beta0_hat.bs <- coef(model.bs,n1=g)$smterms$`Intercept(yindex)`$coef$se
  M0[b] <- max(abs((beta0_hat.bs - beta0_hat)/SE.beta0_hat.bs))
  print(b)
}
save.image(file='myEnvironment.RData')
#hist(M0)
alpha <- .05
q <- quantile(M0,1-alpha, na.rm=T)
CB0.lower <- beta0_hat - q*coef(model,n1=g)$smterms$`Intercept(yindex)`$coef$se
CB0.upper <- beta0_hat + q*coef(model,n1=g)$smterms$`Intercept(yindex)`$coef$se


#point estimate and CIs for the coefficient of var.mat
plot(coef(model, n1=g)$smterms$`Intercept(yindex)`$coef$yindex.vec,
     coef(model, n1=g)$smterms$`Intercept(yindex)`$coef$value, type='l'
     ,ylim=c(-100,100), xlab='r', ylab='L_7-3(r) - r')
lines(coef(model, n1=g)$smterms$`Intercept(yindex)`$coef$yindex.vec,
      coef(model, n1=g)$smterms$`Intercept(yindex)`$coef$value - qnorm(.025, lower.tail=F)*coef(model, n1=g)$smterms$`Intercept(yindex)`$coef$se
      , col='blue')
lines(coef(model, n1=g)$smterms$`Intercept(yindex)`$coef$yindex.vec,
      coef(model, n1=g)$smterms$`Intercept(yindex)`$coef$value + qnorm(.025, lower.tail=F)*coef(model, n1=g)$smterms$`Intercept(yindex)`$coef$se
      , col='blue')
lines(coef(model, n1=g)$smterms$`Intercept(yindex)`$coef$yindex.vec, CB0.lower, col='magenta')
lines(coef(model, n1=g)$smterms$`Intercept(yindex)`$coef$yindex.vec, CB0.upper, col='magenta')
abline(h=0, col='black')
for(i in 1:nPatients){
  lines(coef(model, n1=g)$smterms$`Intercept(yindex)`$coef$yindex.vec, Y.mat[i,], col='grey')
}






##plot images of the observed interaction patterns
#library(ggplot2)
#for(i in unique(data$image_number)){
#  plot.data <- data[data$image_number == i,]
#  plot.data$cell_metacluster <- factor(plot.data$cell_metacluster)
#  png(paste0('image_id_',i,'.png'))
#  print(ggplot(data=plot.data, aes(x = cell_x, y = cell_y, col=cell_metacluster)) +
#    geom_point() + theme_bw() + ggtitle(paste0('image_id ',i)))
#  dev.off()
#}
#data$cells_interested <- ifelse(data$cell_metacluster %in% c(9,15), data$cell_metacluster, 0)
#for(i in unique(data$image_number)){
#  plot.data <- data[data$image_number == i,]
#  plot.data$cells_interested <- factor(plot.data$cells_interested)
#  png(paste0('image_id_',i,'.png'))
#  print(ggplot(data=plot.data, aes(x = cell_x, y = cell_y, col=cells_interested)) +
#          geom_point() + theme_bw() + ggtitle(paste0('image_id ',i)))
#  dev.off()
#}




















