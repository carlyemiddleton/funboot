#' @title Plots the locations of cells within a set of images, colored by cell phenotype.
#'
#' @param data Data frame containing the cell locations and cell types to be plotted.  It should contain the variables `image_number`, `cell_x`, `cell_y`, and `cell_type` as in the example data.
#' @param images_to_plot Vector containing the values of `image_number` corresponding to the images to be plotted. For example, `1:10`
#' @param cell_types_to_plot Vector containing the values of `cell_type` corresponding to the cell types to be plotted with color.  For example, `c(7,3,9,16)`
#' @param palette Vector containing assignments of `cell_type` values to their desired color.  One element should be titled `"Other"` and specify the desired color for all other cells.  For example, `c("Other" = "gray", "7" = "gold", "3" = "darkgreen", '9' = "skyblue", "16"="cyan4")`.
#'
#' @return A list containing `ggplot` objects for each image to be plotted
#'
#' @export

plot_images <- function(data, images_to_plot, cell_types_to_plot, palette){

  data$cell_graphlabels <- ifelse(data[[paste0("cell_type")]] %in% cell_types_to_plot,
                                  data[[paste0("cell_type")]],
                                  'Other')
  data$cell_graphlabels <- factor(data$cell_graphlabels)
  plot_list <- list()
  for(image in images_to_plot){
    assign(paste0('p',image),
           ggplot2::ggplot(data=data[data['image_number']==image,]) +
             ggplot2::theme_bw() +
             ggplot2::geom_point(ggplot2::aes(x=cell_x,
                            y=cell_y,
                            col=cell_graphlabels),
                        size=.5) +
             ggplot2::labs(title=paste0('Image ', image),
                  x=expression('X Coordinate (' ~ mu ~ 'm)'),
                  y=expression('Y Coordinate (' ~ mu ~ 'm)'),
                  col='Cell Type') +
             ggplot2::scale_color_manual(values = palette) +
             ggplot2::theme(plot.title = ggplot2::element_text(size = 20))
    )
    len <- length(plot_list)
    plot_list[[len+1]] <- get(paste0('p',image))
  }
  return(plot_list)
}


