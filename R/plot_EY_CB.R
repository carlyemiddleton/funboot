#' @title Plots wild bootstrap confidence bands obtained from *funboot::EY_CB()*
#'
#' @param CB.object List obtained from the output of *funboot::EY_CB()*
#'
#' @return A *ggplot* object
#'
#' @export

plot_EY_CB <- function(CB.object=CBs){
  print(
    ggplot() + theme_bw() +
          geom_line(aes(x=CB.object$grid,y=CB.object$estimate, col='Estimate')) +
          geom_line(aes(x=CB.object$grid,y=CB.object$CB.lower, col='Wild Bootstrap CB')) +
          geom_line(aes(x=CB.object$grid,y=CB.object$CB.upper, col='Wild Bootstrap CB')) +
          labs(x='r', y=expression(E* "[" * Y(r) * "]"), col=' ') + theme(axis.title.y = element_text(size = 20),
                                                        axis.title.x = element_text(size = 20)) +
          geom_hline(yintercept=0, lty=2) +
          scale_color_manual(values = c("Estimate" = '#F8766D',
                                        "Wild Bootstrap CB"='#00BA38'))
  )

}

