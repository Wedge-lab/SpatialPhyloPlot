# Function to bring together plotting and tree wrangling

#' Title
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
  # newick file is checked in newick function

  ############# Tree functions
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
  # Newick tree to graph df
  newick_tree_df <- newick_to_graph_df(newick_file)
  # Calculate centroids
  newick_tree_df <- calculate_centroids(newick_df = newick_tree_df,
                                        coordinates_df_scaled = coords)
  # Create segments
  segments <- create_segments(newick_tree_df)

  ############# Plotting functions

  # get raster image
  raster_image <- img_get_raster(image_file)

  # get colour palette
  if(any(palette == "default")){
    my_palette <- img_name_colours(newick_tree_df, ...)
  }else{
    my_palette <- palette
  }


  plot <- img_plot(
    raster_img = raster_image,
    coordinates_df_scaled = coords,
    newick_df = newick_tree_df,
    from_to_df = segments,
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
    fig_offset_y = fig_offset_y
  )

  return(plot)
}
