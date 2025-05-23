# Function to bring together plotting and tree wrangling

#' Plot spatial phylogenetic trees
#'
#' A function for plotting spatial phylogenetic trees on Visium 10X data.
#'
#' @param visium_object A `Seurat` Visium 10X object for the tumour sample of interest.
#' @param newick_file A file path to the phylogenetic tree in the Newick format saved as a `.new` file.
#' @param image_file A file path to the hires Visium 10X png image to plot.
#' @param clone_df A `data.frame` listing which barcode belongs to which clone.
#' @param clone_group_column The name of the column in `clone_df` that indicates the clones.
#' @param clone_barcode_column The name of the column in `clone_df` that indicates the barcode.
#'
#' @inheritParams img_plot
#'
#' @param visium_object_left If running in multisample mode; a `Seurat` Visium 10X object for the tissue to be plotted to the left of the main sample.
#' @param visium_object_right If running in multisample mode; a `Seurat` Visium 10X object for the tissue to be plotted to the right of the main sample.
#' @param visium_object_top If running in multisample mode; a `Seurat` Visium 10X object for the tissue to be plotted above the main sample.
#' @param visium_object_bottom If running in multisample mode; a `Seurat` Visium 10X object for the tissue to be plotted below the main sample.
#' @param image_file_left If running in multisample mode; a file path to the hires Visium 10X png image to plot to the left of the main sample.
#' @param image_file_right If running in multisample mode; a file path to the hires Visium 10X png image to plot to the right of the main sample.
#' @param image_file_top If running in multisample mode; a file path to the hires Visium 10X png image to plot above the main sample.
#' @param image_file_bottom If running in multisample mode; a file path to the hires Visium 10X png image to plot below the main sample.
#' @param clone_df_left If running in multisample mode; a  `data.frame` listing which barcode belongs to which clone belonging to the sample plotted to the left of the main sample.
#' @param clone_df_right If running in multisample mode; a  `data.frame` listing which barcode belongs to which clone belonging to the sample plotted to the right of the main sample.
#' @param clone_df_top If running in multisample mode; a  `data.frame` listing which barcode belongs to which clone belonging to the sample plotted above the main sample.
#' @param clone_df_bottom If running in multisample mode; a  `data.frame` listing which barcode belongs to which clone belonging to the sample plotted below the main sample.
#'
#' @import ggplot2
#' @import ggforce
#'
#' @returns A ggplot2 object with the tissue, clones, and phylogenetic trees plotted together.
#' @export
#'
#' @examples
SpatialPhyloPlot <- function(visium_object,
                             newick_file,
                             image_file,
                             clone_df,
                             clone_group_column,
                             clone_barcode_column = "Barcode",
                             palette = "default",
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
                             fig_offset_x = 0.011,
                             fig_offset_y = 0.011,
                             # Multisample options
                             multisample = FALSE,
                             plot_connections = FALSE,
                             visium_object_left = NA,
                             visium_object_right = NA,
                             visium_object_top = NA,
                             visium_object_bottom = NA,
                             image_file_left = NA,
                             image_file_right = NA,
                             image_file_top = NA,
                             image_file_bottom = NA,
                             clone_df_left = NA,
                             clone_df_right = NA,
                             clone_df_top = NA,
                             clone_df_bottom = NA,
                             connection_width = 1,
                             connection_colour = "grey",
                             ...) {
  ############# Check inputs
  if(all(!("Seurat" %in% class(visium_object)))){
    stop("visium_object input must be of class 'Seurat'")
  }
  if(!(clone_group_column %in% colnames(clone_df))){
    stop("clone_group_column must refer to a column present in clone_df.")
  }
  if(class(clone_df) != "data.frame"){
    stop("clone_df must be an object of class 'data.frame'.")
  }
  if(multisample){
    if(all(is.na(c(visium_object_left, visium_object_right, visium_object_top, visium_object_bottom)))){
      stop("Must provide at least one of 'c(visium_object_left, visium_object_right, visium_object_top, visium_object_bottom)' if running in multisample mode.")
    }
    if(all(is.na(c(image_file_left, image_file_right, image_file_top, image_file_bottom)))){
      stop("Must provide at least one of 'c(image_file_left, image_file_right, image_file_top, image_file_bottom)' if running in multisample mode.")
    }
    if(all(is.na(c(clone_df_left, clone_df_right, clone_df_top, clone_df_bottom)))){
      stop("Must provide at least one of 'c(clone_df_left, clone_df_right, clone_df_top, clone_df_bottom)' if running in multisample mode.")
    }
  }
  if(plot_connections){
    if(!multisample){
      stop("Must set 'multisample = TRUE' and provide multiple samples to plot connections between clones in samples. ")
    }
  }
  if(!multisample){
    if(any(!is.na(c(visium_object_left, visium_object_right, visium_object_top, visium_object_bottom,
                    image_file_left, image_file_right, image_file_top, image_file_bottom,
                    clone_df_left, clone_df_right, clone_df_top, clone_df_bottom)))){
      stop("Must run 'multisample = TRUE' if providing more than one sample. ")
    }
  }
  # newick file is checked in newick function

  ############# Tree functions
  ############# Get coordinates and match - base/central sample
  # get coordinates
  coords <- as.data.frame(Seurat::GetTissueCoordinates(visium_object))
  # rotate coordinates
  coords <- rotate_coordinates(coords)
  # scale coordinates
  coords <- scale_coordinates(coords, visium = visium_object)
  # Match to list of clone barcodes
  coords <- match_clone_barcodes(coordinates_df_scaled = coords,
                                 clones_df = clone_df,
                                 clone_group_name = clone_group_column,
                                 coordinate_barcode_name = "cell",
                                 clone_barcode_name = clone_barcode_column)

  if(multisample){

    if(!is.na(visium_object_left)){
      # get coordinates
      coords_left <- as.data.frame(Seurat::GetTissueCoordinates(visium_object_left))
      # rotate coordinates
      coords_left <- rotate_coordinates(coords_left)
      # scale coordinates
      coords_left <- scale_coordinates(coords_left, visium = visium_object_left)

      # Adjust coordinates for plot location
      coords_left$new_x_scaled <- coords_left$new_x_scaled - (1+fig_offset_x)

      # Match to list of clone barcodes
      coords_left <- match_clone_barcodes(coordinates_df_scaled = coords_left,
                                          clones_df = clone_df_left,
                                          clone_group_name = clone_group_column,
                                          coordinate_barcode_name = "cell",
                                          clone_barcode_name = clone_barcode_column)
    }else{
      coords_left <- NA
    }
    if(!is.na(visium_object_right)){
      # get coordinates
      coords_right <- as.data.frame(Seurat::GetTissueCoordinates(visium_object_right))
      # rotate coordinates
      coords_right <- rotate_coordinates(coords_right)
      # scale coordinates
      coords_right <- scale_coordinates(coords_right, visium = visium_object_right)

      # Adjust coordinates for plot location
      coords_right$new_x_scaled <- coords_right$new_x_scaled + 1 + fig_offset_x

      # Match to list of clone barcodes
      coords_right <- match_clone_barcodes(coordinates_df_scaled = coords_right,
                                           clones_df = clone_df_right,
                                           clone_group_name = clone_group_column,
                                           coordinate_barcode_name = "cell",
                                           clone_barcode_name = clone_barcode_column)
    }else{
      coords_right <- NA
    }
    if(!is.na(visium_object_top)){
      # get coordinates
      coords_top <- as.data.frame(Seurat::GetTissueCoordinates(visium_object_top))
      # rotate coordinates
      coords_top <- rotate_coordinates(coords_top)
      # scale coordinates
      coords_top <- scale_coordinates(coords_top, visium = visium_object_top)

      # Adjust coordinates for plot location
      coords_top$new_y_scaled <- coords_top$new_y_scaled + 1

      # Match to list of clone barcodes
      coords_top <- match_clone_barcodes(coordinates_df_scaled = coords_top,
                                         clones_df = clone_df_top,
                                         clone_group_name = clone_group_column,
                                         coordinate_barcode_name = "cell",
                                         clone_barcode_name = clone_barcode_column)
    }else{
      coords_top <- NA
    }
    if(!is.na(visium_object_bottom)){
      # get coordinates
      coords_bottom <- as.data.frame(Seurat::GetTissueCoordinates(visium_object_bottom))
      # rotate coordinates
      coords_bottom <- rotate_coordinates(coords_bottom)
      # scale coordinates
      coords_bottom <- scale_coordinates(coords_bottom, visium = visium_object_bottom)

      # Adjust coordinates for plot location
      coords_bottom$new_y_scaled <- coords_bottom$new_y_scaled - 1

      # Match to list of clone barcodes
      coords_bottom <- match_clone_barcodes(coordinates_df_scaled = coords_bottom,
                                            clones_df = clone_df_bottom,
                                            clone_group_name = clone_group_column,
                                            coordinate_barcode_name = "cell",
                                            clone_barcode_name = clone_barcode_column)
    }else{
      coords_bottom <- NA
    }
  }else{
    coords_left<- NA
    coords_right <- NA
    coords_top <- NA
    coords_bottom <- NA
  }

  ############# Tree match
  # Newick tree to graph df
  newick_tree_df <- newick_to_graph_df(newick_file)

  # Calculate centroids
  ## multisample first as single sample overwritten as newick_tree_df
  if(multisample){
    if(!all(is.na(coords_left))){
      newick_df_left <- calculate_centroids(newick_df = newick_tree_df,
                                            coordinates_df_scaled = coords_left)
    }else{
      newick_df_left <- NA
    }
    if(!all(is.na(coords_right))){
      newick_df_right <- calculate_centroids(newick_df = newick_tree_df,
                                            coordinates_df_scaled = coords_right)
    }else{
      newick_df_right <- NA
    }
    if(!all(is.na(coords_top))){
      newick_df_top <- calculate_centroids(newick_df = newick_tree_df,
                                            coordinates_df_scaled = coords_top)
    }else{
      newick_df_top <- NA
    }
    if(!all(is.na(coords_bottom))){
      newick_df_bottom <- calculate_centroids(newick_df = newick_tree_df,
                                            coordinates_df_scaled = coords_bottom)
    }else{
      newick_df_bottom <- NA
    }
  }else{
    newick_df_left <- NA
    newick_df_right <- NA
    newick_df_top <- NA
    newick_df_bottom <- NA
  }

  newick_tree_df <- calculate_centroids(newick_df = newick_tree_df,
                                        coordinates_df_scaled = coords)

  # Create segments
  segments <- create_segments(newick_tree_df)

  if(multisample){
    if(!all(is.na(newick_df_left))){
      segments_left <- create_segments(newick_df_left)
    }else{
      segments_left <- NA
    }
    if(!all(is.na(newick_df_right))){
      segments_right <- create_segments(newick_df_right)
    }else{
      segments_right <- NA
    }
    if(!all(is.na(newick_df_top))){
      segments_top <- create_segments(newick_df_top)
    }else{
      segments_top <- NA
    }
    if(!all(is.na(newick_df_bottom))){
      segments_bottom <- create_segments(newick_df_bottom)
    }else{
      segments_bottom <- NA
    }
  }else{
    segments_left <- NA
    segments_right <- NA
    segments_top <- NA
    segments_bottom <- NA
  }

  ############# Plotting functions

  # get raster image
  raster_image <- img_get_raster(image_file)

  if(!is.na(image_file_left)){
    raster_image_left <- img_get_raster(image_file_left)
  }else{
    raster_image_left <- NA
  }
  if(!is.na(image_file_right)){
    raster_image_right <- img_get_raster(image_file_right)
  }else{
    raster_image_right <- NA
  }
  if(!is.na(image_file_top)){
    raster_image_top <- img_get_raster(image_file_top)
  }else{
    raster_image_top <- NA
  }
  if(!is.na(image_file_bottom)){
    raster_image_bottom <- img_get_raster(image_file_bottom)
  }else{
    raster_image_bottom <- NA
  }

  # combine newick objects for multisample
  if(multisample){
    multi_newick <- newick_tree_df
    multi_newick$Origin <- "Centre"

    if(!all(is.na(newick_df_left))){
      newick_df_left$Origin <- "Left"
      multi_newick <- rbind(multi_newick, newick_df_left)
    }
    if(!all(is.na(newick_df_right))){
      newick_df_right$Origin <- "Right"
      multi_newick <- rbind(multi_newick, newick_df_right)
    }
    if(!all(is.na(newick_df_top))){
      newick_df_top$Origin <- "Top"
      multi_newick <- rbind(multi_newick, newick_df_top)
    }
    if(!all(is.na(newick_df_bottom))){
      newick_df_bottom$Origin <- "Bottom"
      multi_newick <- rbind(multi_newick, newick_df_bottom)
    }
  }

  # get connections between plots for multisample
  if(multisample){
    connections_coords <- connect_multisample(multi_newick)
  }else{
    connections_coords <- NA
  }


  # get colour palette
  if(any(palette == "default")){
    my_palette <- img_name_colours(newick_tree_df)

    if(multisample){
      my_palette <- img_name_colours(multi_newick, ...)
    }
  }else{
    my_palette <- palette
  }


  plot <- img_plot(
    raster_img = raster_image,
    raster_img_left = raster_image_left,
    raster_img_right = raster_image_right,
    raster_img_top = raster_image_top,
    raster_img_bottom = raster_image_bottom,
    coordinates_df_scaled = coords,
    coordinates_df_scaled_left = coords_left,
    coordinates_df_scaled_right = coords_right,
    coordinates_df_scaled_top = coords_top,
    coordinates_df_scaled_bottom = coords_bottom,
    newick_df = newick_tree_df,
    newick_df_left = newick_df_left,
    newick_df_right = newick_df_right,
    newick_df_top = newick_df_top,
    newick_df_bottom = newick_df_bottom,
    from_to_df = segments,
    from_to_df_left = segments_left,
    from_to_df_right = segments_right,
    from_to_df_top = segments_top,
    from_to_df_bottom = segments_bottom,
    palette = my_palette,
    plot_points = plot_points,
    plot_polygon = plot_polygon,
    point_alpha = point_alpha,
    hull_alpha = hull_alpha,
    hull_expansion = hull_expansion,
    centroid_alpha = centroid_alpha,
    centroid_size = centroid_size,
    segment_alpha = segment_alpha,
    segment_width = segment_width,
    segment_colour = segment_colour,
    fig_offset_x = fig_offset_x,
    fig_offset_y = fig_offset_y,
    multisample = multisample,
    plot_connections = plot_connections,
    connections_coords = connections_coords,
    connection_colour = connection_colour,
    connection_width = connection_width
  )

  return(plot)
}
