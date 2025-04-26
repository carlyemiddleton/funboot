#' @title Plots wild bootstrap confidence bands obtained from *funboot::wildBS_CB(),* along with pointwise confidence intervals
#'
#' @param CB.object List obtained from the output of *funboot::wildBS_CB()*
#' @param alpha Desired significance level for the pointwise confidence intervals.  For example, 0.05.  The pointwise confidence intervals are obtained using the asymptotic normality of the statistic \eqn{\frac{\hat{\beta}_{j}(r)}{SE(\hat{\beta}_{j}(r))}} as calculated from *pffr()* output, and are not adjusted for multiplicity.
#'
#' @return A *ggplot* object for each model coefficient function
#'
#' @export

plot_wildBS_CB <- function(CB.object, alpha=.05){

  for(i in 1:dim(CB.object$CB.lower)[2]){
    print(ggplot() + theme_bw() +
            geom_line(aes(x=CB.object$grid,y=CB.object[[4]][,i], col='Estimate')) +
            geom_line(aes(x=CB.object$grid,y=CB.object[[4]][,i]- qnorm(1-alpha/2)*CB.object[[5]][,i], col='Pointwise CI (unadjusted)')) +
            geom_line(aes(x=CB.object$grid,y=CB.object[[4]][,i]+ qnorm(1-alpha/2)*CB.object[[5]][,i], col='Pointwise CI (unadjusted)')) +
            geom_line(aes(x=CB.object$grid,y=CB.object[[1]][,i], col='Wild Bootstrap CB')) +
            geom_line(aes(x=CB.object$grid,y=CB.object[[2]][,i], col='Wild Bootstrap CB')) +
            labs(x='r', y='coefficient', col=' ') + theme(axis.title.y = element_text(size = 20),
                                                          axis.title.x = element_text(size = 20)) +
            geom_hline(yintercept=0, lty=2) + ggtitle(paste0(names(CB.object$CB.lower)[i]))  +
            scale_color_manual(values = c("Estimate" = '#F8766D',
                                          "Pointwise CI (unadjusted)"='#00BFC4',
                                          "Wild Bootstrap CB"='#00BA38'))  )
  }
}

