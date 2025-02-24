#' @title Preprocesses data
#'
#' @param grid the grid that the CBs were evaluated over.  eg. 0:200
#' @return It returns res.wildBS
#'
#' @export

plot_images <- function(data, images.to.plot, cell.types.to.plot, pallete){

  data$cell_graphlabels <- ifelse(data[[paste0("cell_type")]] %in% cell.types.to.plot,
                                  data[[paste0("cell_type")]], 'Other')
  data$cell_graphlabels <- factor(data$cell_graphlabels)
  plot.list <- list()
  for(image in images.to.plot){
    assign(paste0('p',image),
           ggplot(data=data[data['image_number']==image,]) + theme_bw() +
             geom_point(aes(x=cell_x, y=cell_y, col=cell_graphlabels),size=.5) +
             labs(title=paste0('Image ', image),x=expression('X Coordinate (' ~ mu ~ 'm)'),
                  y=expression('Y Coordinate (' ~ mu ~ 'm)'), col='Cell Type') +
             scale_color_manual(values = pallete) +
             theme(plot.title = element_text(size = 20))
    )
    len <- length(plot.list)
    plot.list[[len+1]] <- get(paste0('p',image))
  }
  return(plot.list)
}


