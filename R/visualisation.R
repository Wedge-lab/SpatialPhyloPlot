# Plotting functions

#' Get Raster Image
#'
#' @param image_path File path to the underlying Visium 10X tissue image, which will be plotted.
#'
#' @returns An image magick object which is used to plot the tissue image.
#'
img_get_raster <- function(image_path) {

  if(!file.exists(image_path)){
    stop("Image file cannot be found.")
  }

  hires <- magick::image_read(image_path)
  # hires_raster <- as.raster(hires)

  return(hires)
}

#' Crop raster image
#'
#' @param hires_scale High resolution scale factor extracted from Seurat object.
#' @param tissue_positions A `data.frame` with the pixel coordinates of each spot.
#' @param raster An image magick object that is plotted as the tissue.
#'
#' @returns A magick image object cropped to the appropriate bounding box.

img_crop_raster <- function(raster, hires_scale, tissue_positions) {

  hires_info <- magick::image_info(raster)

  from_left <- min(tissue_positions$pxl_col_in_fullres)*hires_scale
  from_bottom <- (min(tissue_positions$pxl_row_in_fullres)*hires_scale)
  width <- (max(tissue_positions$pxl_col_in_fullres)*hires_scale) - (min(tissue_positions$pxl_col_in_fullres)*hires_scale)
  height <- (max(tissue_positions$pxl_row_in_fullres)*hires_scale) - (min(tissue_positions$pxl_row_in_fullres)*hires_scale)

  hires_crop <- magick::image_crop(raster, geometry =paste0(width,"x",height,"+",from_left,"+",from_bottom))

  return(hires_crop)

}

#' Name colours for plotting
#'
#' @param newick_df A data frame generated from the .new phylogenetic tree file.
#'
#' @import pals
#'
#' @returns A named list with the colour palette associated with the list of clones.
#'
img_name_colours <- function(newick_df) {
  nclones <- length(unique(na.omit(newick_df$colour)))

  long_palette <- c(pals::brewer.set2(n = 8),
                    pals::alphabet(n = 26),
                    pals::alphabet2(n = 26))

  mypal <- long_palette[1:nclones]
  names(mypal) <- unique(na.omit(newick_df$colour))

  return(mypal)
}

#' Plot tree and spots
#'
#' @param raster_img The raster image object containing the tissue image to be plotted.
#' @param coordinates_df_scaled A `data.frame` with the scaled (0-1) spot coordinates for each Visium barcode.
#' @param newick_df The newick `data.frame` generated from the .new phylogenetic tree file.
#' @param from_to_df A `data.frame` containing the start and end coordinates of each branch on the phylogenetic tree.
#' @param point_alpha The alpha transparency value for points, default is `0.8`.
#' @param hull_alpha The alpha transparency value for polygons/hulls to be plotted, default is `0` (not plotted).
#' @param hull_expansion The hull/polygon extension amount, default is `0.005`.
#' @param centroid_alpha The alpha transparency value for the centroids plotted as part of the phylogenetic tree, default is `0.9`.
#' @param centroid_size The size of the plotted centroid points, default is `8`.
#' @param segment_alpha The alpha transparency value for phylogenetic tree segments, default is `0.8`.
#' @param segment_width The width of the phylogenetic tree segments, default is `2`.
#' @param segment_colour The colour of the phylogenetic tree segments, default is `"grey"`.
#' @param fig_offset_x The offset to adjust the histology image x_end location, if required. Default is `0.011`, adjustments not recommended.
#' @param fig_offset_y The offset to adjust the histology image y location, if required. Default is `0.011` adjustments not recommended.
#' @param palette A named list of colours for each clone. Defaults to `RColorBrewer` `Set2` followed by the `pals` `alphabet` sets.
#' @param raster_img_left The raster image object containing the tissue image to be plotted to the left of the main image.
#' @param raster_img_right The raster image object containing the tissue image to be plotted to the right of the main image.
#' @param raster_img_top The raster image object containing the tissue image to be plotted above the main image.
#' @param raster_img_bottom The raster image object containing the tissue image to be plotted below the main image.
#' @param coordinates_df_scaled_left A `data.frame` with the scaled (0-1) spot coordinates for each Visium barcode to be plotted to the left of the main tissue.
#' @param coordinates_df_scaled_right A `data.frame` with the scaled (0-1) spot coordinates for each Visium barcode to be plotted to the right of the main tissue.
#' @param coordinates_df_scaled_top A `data.frame` with the scaled (0-1) spot coordinates for each Visium barcode to be plotted above the main tissue.
#' @param coordinates_df_scaled_bottom A `data.frame` with the scaled (0-1) spot coordinates for each Visium barcode to be plotted below the main tissue.
#' @param newick_df_left The newick `data.frame` generated from the .new phylogenetic tree file with connections to plot to the left of the main tissue.
#' @param newick_df_right The newick `data.frame` generated from the .new phylogenetic tree file with connections to plot to the right of the main tissue.
#' @param newick_df_top The newick `data.frame` generated from the .new phylogenetic tree file with connections to plot above the main tissue.
#' @param newick_df_bottom The newick `data.frame` generated from the .new phylogenetic tree file with connections to plot below the main tissue.
#' @param from_to_df_left A `data.frame` containing the start and end coordinates of each branch on the phylogenetic tree to be plotted to the left of the main tissue.
#' @param from_to_df_right A `data.frame` containing the start and end coordinates of each branch on the phylogenetic tree to be plotted to the right of the main tissue.
#' @param from_to_df_top A `data.frame` containing the start and end coordinates of each branch on the phylogenetic tree to be plotted above the main tissue.
#' @param from_to_df_bottom A `data.frame` containing the start and end coordinates of each branch on the phylogenetic tree to be plotted below the main tissue.
#' @param plot_points Whether individual Visium spots should be plotted. Defaults to `TRUE`.
#' @param plot_polygon Whether a polygon should be plotted for each clone. Defaults to `FALSE`.
#' @param multisample Whether multiple samples are being plotted with the same phylogenetic tree. Defaults to `FALSE`.
#' @param plot_connections Whether connections should be plotted between clones that are present in multiple samples. Defaults to `FALSE`.
#' @param connections_coords A `data.frame` of coordinates linking clones found in multiple samples.
#' @param connection_width The width of the connecting segment linking clones between samples. Defaults to `1`.
#' @param connection_colour The colour of the connecting segment linking clones between samples. Defaults to `grey`.
#'
#' @import ggplot2
#' @import ggforce
#'
#' @returns A ggplot2 object with clones, trees, and tissue image plotted together.
#'
img_plot <- function(raster_img,
                     raster_img_left,
                     raster_img_right,
                     raster_img_top,
                     raster_img_bottom,
                     coordinates_df_scaled,
                     coordinates_df_scaled_left,
                     coordinates_df_scaled_right,
                     coordinates_df_scaled_top,
                     coordinates_df_scaled_bottom,
                     newick_df,
                     newick_df_left,
                     newick_df_right,
                     newick_df_top,
                     newick_df_bottom,
                     from_to_df,
                     from_to_df_left,
                     from_to_df_right,
                     from_to_df_top,
                     from_to_df_bottom,
                     plot_points = TRUE,
                     plot_polygon = FALSE,
                     point_alpha = 0.8,
                     hull_alpha = 0,
                     hull_expansion = 0.005,
                     centroid_alpha = 0.9,
                     centroid_size = 8,
                     segment_alpha = 0.8,
                     segment_width = 2,
                     segment_colour = "grey",
                     fig_offset_x = 0,
                     fig_offset_y = 0,
                     palette,
                     multisample,
                     shared_clones,
                     plot_connections,
                     connections_coords,
                     connection_width,
                     connection_colour
) {

  # adjust plot coordinates if multiple plots
  plot_xlim <- c(0,1)
  plot_ylim <- c(0,1)

  if(!all(is.na(raster_img_left))){
    plot_xlim[1] <- -1
  }
  if(!all(is.na(raster_img_right))){
    plot_xlim[2] <- 2
  }
  if(!all(is.na(raster_img_top))){
    plot_ylim[2] <- 2
  }
  if(!all(is.na(raster_img_bottom))){
    plot_ylim[1] <- -1
  }

  # Main plot
  p <- ggplot()+
    coord_fixed(xlim = plot_xlim, ylim = plot_ylim)+
    annotation_raster(raster_img, xmin = 0,  xmax = 1+fig_offset_x,
                      ymin = 0+fig_offset_y, ymax = 1)

  # options for extra tissue images
  if(!all(is.na(raster_img_left))){
    p <- p + annotation_raster(raster_img_left, xmin = -(1+fig_offset_x),  xmax = 0,
                               ymin = 0+fig_offset_y, ymax = 1)
  }
  if(!all(is.na(raster_img_right))){
    p <- p + annotation_raster(raster_img_right, xmin = 1+fig_offset_x,  xmax = 2+(2*fig_offset_x),
                               ymin = 0+fig_offset_y, ymax = 1)
  }
  if(!all(is.na(raster_img_top))){
    p <- p + annotation_raster(raster_img_top, xmin = 0,  xmax = 1+fig_offset_x,
                               ymin = 1+fig_offset_y, ymax = 2)
  }
  if(!all(is.na(raster_img_bottom))){
    p <- p + annotation_raster(raster_img_bottom, xmin = 0,  xmax = 1+fig_offset_x,
                               ymin = -1+fig_offset_y, ymax = 0)
  }

  # plotting points
  if(plot_points){
    p <- p +
      geom_point(data = coordinates_df_scaled, aes(x = new_x_scaled, y = new_y_scaled, colour = Clone),
                 alpha = point_alpha, show.legend = NA)
    # plotting multisample points
    if(multisample & !shared_clones){
      if(!all(is.na(coordinates_df_scaled_left))){
        p <- p +
          geom_point(data = coordinates_df_scaled_left, aes(x = new_x_scaled, y = new_y_scaled, colour = Clone),
                     alpha = point_alpha, show.legend = NA)
      }
      if(!all(is.na(coordinates_df_scaled_right))){
        p <- p +
          geom_point(data = coordinates_df_scaled_right, aes(x = new_x_scaled, y = new_y_scaled, colour = Clone),
                     alpha = point_alpha, show.legend = NA)
      }
      if(!all(is.na(coordinates_df_scaled_top))){
        p <- p +
          geom_point(data = coordinates_df_scaled_top, aes(x = new_x_scaled, y = new_y_scaled, colour = Clone),
                     alpha = point_alpha, show.legend = NA)
      }
      if(!all(is.na(coordinates_df_scaled_bottom))){
        p <- p +
          geom_point(data = coordinates_df_scaled_bottom, aes(x = new_x_scaled, y = new_y_scaled, colour = Clone),
                     alpha = point_alpha, show.legend = NA)
      }
    }
  }

  # If plotting polygons
  if(plot_polygon){
    p <- p +
      geom_mark_hull(data = coordinates_df_scaled, aes(x = new_x_scaled, y = new_y_scaled, fill = Clone),
                     alpha = hull_alpha, expand=hull_expansion, show.legend = NA)

    ## If plotting multisample
    if(multisample & !shared_clones){
      if(!all(is.na(coordinates_df_scaled_left))){
        p <- p +
          geom_mark_hull(data = coordinates_df_scaled_left, aes(x = new_x_scaled, y = new_y_scaled, fill = Clone),
                         alpha = hull_alpha, expand=hull_expansion, show.legend = NA)
      }
      if(!all(is.na(coordinates_df_scaled_right))){
        p <- p +
          geom_mark_hull(data = coordinates_df_scaled_right, aes(x = new_x_scaled, y = new_y_scaled, fill = Clone),
                         alpha = hull_alpha, expand=hull_expansion, show.legend = NA)
      }
      if(!all(is.na(coordinates_df_scaled_top))){
        p <- p +
          geom_mark_hull(data = coordinates_df_scaled_top, aes(x = new_x_scaled, y = new_y_scaled, fill = Clone),
                         alpha = hull_alpha, expand=hull_expansion, show.legend = NA)
      }
      if(!all(is.na(coordinates_df_scaled_bottom))){
        p <- p +
          geom_mark_hull(data = coordinates_df_scaled_bottom, aes(x = new_x_scaled, y = new_y_scaled, fill = Clone),
                         alpha = hull_alpha, expand=hull_expansion, show.legend = NA)
      }
    }
  }

  # plotting segments
  p <- p +
    geom_segment(data = from_to_df, aes(x = from_x, y = from_y, xend = to_x, yend = to_y),
                 colour = segment_colour, linewidth = segment_width, alpha = segment_alpha)

  ## plotting multisample
  if(multisample & !shared_clones){
    if(!all(is.na(from_to_df_left))){
      p <- p +
        geom_segment(data = from_to_df_left, aes(x = from_x, y = from_y, xend = to_x, yend = to_y),
                     colour = segment_colour, linewidth = segment_width, alpha = segment_alpha)
    }
    if(!all(is.na(from_to_df_right))){
      p <- p +
        geom_segment(data = from_to_df_right, aes(x = from_x, y = from_y, xend = to_x, yend = to_y),
                     colour = segment_colour, linewidth = segment_width, alpha = segment_alpha)
    }
    if(!all(is.na(from_to_df_top))){
      p <- p +
        geom_segment(data = from_to_df_top, aes(x = from_x, y = from_y, xend = to_x, yend = to_y),
                     colour = segment_colour, linewidth = segment_width, alpha = segment_alpha)
    }
    if(!all(is.na(from_to_df_bottom))){
      p <- p +
        geom_segment(data = from_to_df_bottom, aes(x = from_x, y = from_y, xend = to_x, yend = to_y),
                     colour = segment_colour, linewidth = segment_width, alpha = segment_alpha)
    }
  }

  # Plotting connections between multisample plots
  if(multisample & !shared_clones){
    if(plot_connections){
      if(!all(is.na(newick_df_left))){
        p <- p +
          geom_segment(data = connections_coords, aes(x = centre_x, y = centre_y, xend = left_x, yend = left_y),
                       colour = segment_colour, linewidth = connection_width, alpha = segment_alpha, linetype = 2)
      }
      if(!all(is.na(newick_df_right))){
        p <- p +
          geom_segment(data = connections_coords, aes(x = centre_x, y = centre_y, xend = right_x, yend = right_y),
                       colour = segment_colour, linewidth = connection_width, alpha = segment_alpha, linetype = 2)
      }
      if(!all(is.na(newick_df_top))){
        p <- p +
          geom_segment(data = connections_coords, aes(x = centre_x, y = centre_y, xend = top_x, yend = top_y),
                       colour = segment_colour, linewidth = connection_width, alpha = segment_alpha, linetype = 2)
      }
      if(!all(is.na(newick_df_bottom))){
        p <- p +
          geom_segment(data = connections_coords, aes(x = centre_x, y = centre_y, xend = bottom_x, yend = bottom_y),
                       colour = segment_colour, linewidth = connection_width, alpha = segment_alpha, linetype = 2)
      }
    }
  }

  # plotting centroids

  p <- p +
    geom_point(data = subset(newick_df, !is.na(colour)), aes(x = centroid_x, y = centroid_y, colour = colour, fill = colour),
               size = centroid_size, alpha = centroid_alpha, show.legend = NA, shape = 21, colour = "black")

  # if multisample
  if(multisample & !shared_clones){
    if(!all(is.na(newick_df_left))){
      p <- p +
        geom_point(data = subset(newick_df_left, !is.na(colour)), aes(x = centroid_x, y = centroid_y, colour = colour, fill = colour),
                   size = centroid_size, alpha = centroid_alpha, show.legend = NA, shape = 21, colour = "black")
    }
    if(!all(is.na(newick_df_right))){
      p <- p +
        geom_point(data = subset(newick_df_right, !is.na(colour)), aes(x = centroid_x, y = centroid_y, colour = colour, fill = colour),
                   size = centroid_size, alpha = centroid_alpha, show.legend = NA, shape = 21, colour = "black")
    }
    if(!all(is.na(newick_df_top))){
      p <- p +
        geom_point(data = subset(newick_df_top, !is.na(colour)), aes(x = centroid_x, y = centroid_y, colour = colour, fill = colour),
                   size = centroid_size, alpha = centroid_alpha, show.legend = NA, shape = 21, colour = "black")
    }
    if(!all(is.na(newick_df_bottom))){
      p <- p +
        geom_point(data = subset(newick_df_bottom, !is.na(colour)), aes(x = centroid_x, y = centroid_y, colour = colour, fill = colour),
                   size = centroid_size, alpha = centroid_alpha, show.legend = NA, shape = 21, colour = "black")
    }
  }


  # Plot tidy and colour
  p <- p +
    xlab("")+
    ylab("")+
    theme_void()+
    # coord_cartesian(xlim = c(0,1),
    #                 ylim = c(0,1))+
    scale_color_manual(values = palette)+
    scale_fill_manual(values = palette)+
    guides(colour = "none", fill = guide_legend("Clone"))

  return(p)

}


# gif version
