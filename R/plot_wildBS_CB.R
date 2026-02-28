#' @title Plots wild bootstrap confidence bands obtained from *funboot::EY_CB()*
#'
#' @param CB_object List obtained from the output of *funboot::EY_CB()*
#'
#' @return A *ggplot* object
#'
#' @export

plot_wildBS_CB <- function(CB_object){
  print(
    ggplot2::ggplot() +
      ggplot2::theme_bw() +
      ggplot2::geom_line(ggplot2::aes(x=CB_object$grid,
                    y=CB_object$target_curve_estimates,
                    col='Estimate')
                ) +
      ggplot2::geom_line(ggplot2::aes(x=CB_object$grid,
                    y=CB_object$CB_lower,
                    col='Wild Bootstrap CB')
                ) +
      ggplot2::geom_line(ggplot2::aes(x=CB_object$grid,
                    y=CB_object$CB_upper,
                    col='Wild Bootstrap CB')
                ) +
      ggplot2::labs(x='Radius',
           y="Target Curve",
           col=' ') +
      ggplot2::theme(axis.title.y = ggplot2::element_text(size = 20),
            axis.title.x = ggplot2::element_text(size = 20)) +
      ggplot2::geom_hline(yintercept=0,
                     lty=2) +
      ggplot2::scale_color_manual(values = c("Estimate" = '#F8766D',
                                        "Wild Bootstrap CB"='#00BA38'))
  )

}

