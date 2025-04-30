# Plotting functions

#' Get Raster Image
#'
#' @param image_path File path to the underlying Visium 10X tissue image, which will be plotted.
#'
#' @returns An image raster object which is used to plot the tissue image.
#'
img_get_raster <- function(image_path) {

  if(!file.exists(image_path)){
    stop("Image file cannot be found.")
  }

  hires <- magick::image_read(image_path)
  hires_raster <- as.raster(hires)

  return(hires_raster)
}



# offset

#' Name colours for plotting
#'
#' @param newick_df A data frame generated from the .new phylogenetic tree file.
#' @param palette_name The name of the RColorBrewer palette to be used.
#'
#' @import RColorBrewer
#'
#' @returns A named list with the colour palette associated with the list of clones.
#'
img_name_colours <- function(newick_df, palette_name = "Set2") {
  nclones <- length(unique(newick_df$colour))
  mypal <- RColorBrewer::brewer.pal(n = nclones, name = palette_name)
  names(mypal) <- unique(newick_df$colour)

  return(mypal)
}

#' Plot tree and spots
#'
#' @param raster_img The raster image object containing the tissue image to be plotted.
#' @param coordinates_df_scaled A `data.frame` with the scaled (0-1) spot coordinates for each Visium barcode.
#' @param newick_df The newick data.frame generated from the .new phylogenetic tree file.
#' @param from_to_df A `data.frame` containing the start and end coordinates of each branch on the phylogenetic tree.
#' @param point_alpha The alpha transparency value for points, default is `0.8`.
#' @param hull_alpha The alpha transparency value for polygons/hulls to be plotted, default is `0` (not plotted).
#' @param hull_expansion The hull/polygon extension amount, default is `0.005`.
#' @param centroid_alpha The alpha transparency value for the centroids plotted as part of the phylogenetic tree, default is `0.9`.
#' @param centroid_size The size of the plotted centroid points, default is `8`.
#' @param segment_alpha The alpha transparency value for phylogenetic tree segments, default is `0.8`.
#' @param segment_width The width of the phylogenetic tree segments, default is `2`.
#' @param segment_colour The colour of the phylogenetic tree segments, default is `"grey"`.
#' @param fig_offset_x The offset to adjust the histology image x location, if required. Default is `0.025`, adjustments not recommended.
#' @param fig_offset_y The offset to adjust the histology image y location, if required. Default is `0.025`, adjustments not recommended.
#' @param palette
#'
#' @import ggplot2
#' @import ggforce
#'
#' @returns A ggplot2 object with clones, trees, and tissue image plotted together.
#'
img_plot <- function(raster_img,
                     coordinates_df_scaled,
                     newick_df,
                     from_to_df,
                     point_alpha = 0.8,
                     hull_alpha = 0,
                     hull_expansion = 0.005,
                     centroid_alpha = 0.9,
                     centroid_size = 8,
                     segment_alpha = 0.8,
                     segment_width = 2,
                     segment_colour = "grey",
                     fig_offset_x = 0.025,
                     fig_offset_y = 0.025,
                     palette
                     ) {

  p <- ggplot()+
    annotation_raster(raster_img, xmin = 0,  xmax = 1+fig_offset_x,
                      ymin = 0, ymax = 1+fig_offset_y)+
    geom_point(data = coordinates_df_scaled, aes(x = new_x_scaled, y = new_y_scaled, colour = clone_bar),
               alpha = point_alpha, show.legend = NA)+
    geom_mark_hull(data = coordinates_df_scaled, aes(x = new_x_scaled, y = new_y_scaled, fill = Clone),
                   alpha = hull_alpha, expand=hull_expansion, show.legend = NA)+
    geom_segment(data = from_to_df, aes(x = from_x, y = from_y, xend = to_x, yend = to_y),
                 colour = segment_colour, linewidth = segment_width, alpha = segment_alpha)+
    geom_point(data = subset(newick_df, !is.na(colour)), aes(x = centroid_x, y = centroid_y, colour = colour),
               size = centroid_size, alpha = centroid_alpha, show.legend = NA)+
    xlab("")+
    ylab("")+
    theme_void()+
    coord_cartesian(xlim = c(0,1),
                    ylim = c(0,1))+
    scale_color_manual(values = palette)+
    scale_fill_manual(values = palette)

  return(p)

}


# gif version
