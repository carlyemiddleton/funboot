#' @title Preprocesses data
#'
#' @param grid the grid that the CBs were evaluated over.  eg. 0:200
#' @return It returns res.wildBS
#'
#' @export

plot.wildBS_CB <- function(CB.object=CBs, alpha=.05){

  plot.list <- list()
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

plot.lin.comb_CB <- function(CB.object=CBs, alpha=.05){
    print(ggplot() + theme_bw() +
            geom_line(aes(x=CB.object$grid,y=CB.object$estimate, col='Estimate')) +
            geom_line(aes(x=CB.object$grid,y=CB.object$estimate- qnorm(1-alpha/2)*CB.object$bs.SE, col='Pointwise CI (unadjusted)')) +
            geom_line(aes(x=CB.object$grid,y=CB.object$estimate+ qnorm(1-alpha/2)*CB.object$bs.SE, col='Pointwise CI (unadjusted)')) +
            geom_line(aes(x=CB.object$grid,y=CB.object$CB.lower, col='Wild Bootstrap CB')) +
            geom_line(aes(x=CB.object$grid,y=CB.object$CB.upper, col='Wild Bootstrap CB')) +
            labs(x='r', y='coefficient', col=' ') + theme(axis.title.y = element_text(size = 20),
                                                          axis.title.x = element_text(size = 20)) +
            geom_hline(yintercept=0, lty=2) + ggtitle('Predicted value of the linear combination of model coefficient functions')  +
            scale_color_manual(values = c("Estimate" = '#F8766D',
                                          "Pointwise CI (unadjusted)"='#00BFC4',
                                          "Wild Bootstrap CB"='#00BA38'))  )

}

