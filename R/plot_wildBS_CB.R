#' @title Plots target curve estimates and wild bootstrap confidence bands obtained from `funboot::wildBS_CB()`
#'
#' @param CB_object List obtained from the output of `funboot::wildBS_CB()`
#'
#' @return A `ggplot` object showing the estimated target curve and its wild bootstrap confidence band
#'
#' @export

plot_wildBS_CB <- function(CB_object){
  p <- ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::geom_line(ggplot2::aes(x = CB_object$grid,
                                    y = CB_object$target_curve_estimates,
                                    col = "Estimated Target Curve")) +
    ggplot2::geom_line(ggplot2::aes(x = CB_object$grid,
                                    y = CB_object$CB_lower,
                                    col = "Confidence Band")) +
    ggplot2::geom_line(ggplot2::aes(x = CB_object$grid,
                                    y = CB_object$CB_upper,
                                    col = "Confidence Band")) +
    ggplot2::labs(x = "Radius",
                  y = "Target Curve",
                  col = " ") +
    ggplot2::theme(axis.title.y = ggplot2::element_text(size = 20),
                   axis.title.x = ggplot2::element_text(size = 20)) +
    ggplot2::geom_hline(yintercept = 0, lty = 2) +
    ggplot2::scale_color_manual(values = c("Estimated Target Curve" = "#F8766D",
                                           "Confidence Band" = "#00BA38"))
  return(p)
}

