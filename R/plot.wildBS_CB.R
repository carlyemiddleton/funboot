#' @title Preprocesses data
#'
#' @param grid the grid that the CBs were evaluated over.  eg. 0:200
#' @return It returns res.wildBS
#'
#' @export

plot.wildBS_CB <- function(CB.object=CBs, alpha=.05){
  re <- CB.object$re
  for(i in 1:dim(CB.object$CB.lower)[2]){
    plot(CB.object$grid, CB.object[[1]][,i], ylim=c(min(CB.object[[1]][,i]),
                                                    min(CB.object[[1]][,i]) +
                  1.5*abs(min(CB.object[[1]][,i])-max(CB.object[[2]][,i]))), col='magenta', main=paste0(names(CB.object[[1]][i])))
    points(CB.object$grid, CB.object[[2]][,i], col='magenta', main=paste0(names(CB.object[[1]][i])))
    points(CB.object$grid, CB.object[[4]][,i], col='black')
    lines(CB.object$grid, CB.object[[4]][,i]- qnorm(1-alpha/2)*CB.object[[5]][,i], col='blue')
    lines(CB.object$grid, CB.object[[4]][,i]+ qnorm(1-alpha/2)*CB.object[[5]][,i], col='blue')
    legend(x='topleft', legend=c('Estimate','Pointwise CI (unadjusted)','Wild BS CB'), fill=c('black','blue','magenta'))
  }
}

