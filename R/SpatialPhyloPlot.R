# Function to bring together plotting and tree wrangling

#' Plot spatial phylogenetic trees
#'
#' A function for plotting spatial phylogenetic trees on Visium 10X data.
#'
#' @param visium_object A `Seurat` Visium 10X object for the tumour sample of interest.
#' @param visium_version Which version of Visium was used. Should be either `"V1"` or `"V2"`
#' @param newick_file A file path to the phylogenetic tree in the Newick format saved as a `.new` file.
#' @param image_file A file path to the hires Visium 10X png image to plot.
#' @param tissue_positions_file A file path to the Visium 10X `tissue_positions_list.csv` file for V1 or `tissue_positions.csv` file for V2.
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
#' @param tissue_positions_file_left If running in multisample mode; a path to the Visium 10X `tissue_positions_list.csv` file for V1 or `tissue_positions.csv` file for V2, belonging to the sample plotted to the left of the main sample.
#' @param tissue_positions_file_right If running in multisample mode; a path to the Visium 10X `tissue_positions_list.csv` file for V1 or `tissue_positions.csv` file for V2, belonging to the sample plotted to the right of the main sample.
#' @param tissue_positions_file_top If running in multisample mode; a path to the Visium 10X `tissue_positions_list.csv` file for V1 or `tissue_positions.csv` file for V2, belonging to the sample plotted above the main sample.
#' @param tissue_positions_file_bottom If running in multisample mode; a path to the Visium 10X `tissue_positions_list.csv` file for V1 or `tissue_positions.csv` file for V2, belonging to the sample plotted below the main sample.
#'
#' @import ggplot2
#' @import ggforce
#' @import Seurat
#'
#' @returns A ggplot2 object with the tissue, clones, and phylogenetic trees plotted together.
#' @export
#'
#' @examples
SpatialPhyloPlot <- function(visium_object,
                             visium_version,
                             newick_file,
                             image_file,
                             tissue_positions_file,
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
                             fig_offset_x = 0,
                             fig_offset_y = 0,
                             # Multisample options
                             multisample = FALSE,
                             shared_clones = FALSE,
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
                             tissue_positions_file_left = NA,
                             tissue_positions_file_right = NA,
                             tissue_positions_file_top = NA,
                             tissue_positions_file_bottom = NA,
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
    if(all(suppressWarnings(is.na(c(visium_object_left, visium_object_right, visium_object_top, visium_object_bottom))))){
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
    if(any(suppressWarnings(!is.na(c(visium_object_left, visium_object_right, visium_object_top, visium_object_bottom,
                                     image_file_left, image_file_right, image_file_top, image_file_bottom,
                                     clone_df_left, clone_df_right, clone_df_top, clone_df_bottom))))){
      stop("Must run 'multisample = TRUE' if providing more than one sample. ")
    }
  }
  # newick file is checked in newick function

  ############# Tree functions
  ############# Get coordinates and match - base/central sample
  # read in tissue positions file
  tpos_colnames <- c("barcode","in_tissue","array_row","array_col",
                     "pxl_row_in_fullres", "pxl_col_in_fullres")
  if(visium_version == "V2"){

    tissue_positions <- read.csv(tissue_positions_file)

    if(multisample){
      if(!is.na(tissue_positions_file_left)){
        tissue_positions_left <- read.csv(tissue_positions_file_left)
      }
      if(!is.na(tissue_positions_file_right)){
        tissue_positions_right <- read.csv(tissue_positions_file_right)
      }
      if(!is.na(tissue_positions_file_top)){
        tissue_positions_top <- read.csv(tissue_positions_file_top)
      }
      if(!is.na(tissue_positions_file_bottom)){
        tissue_positions_bottom <- read.csv(tissue_positions_file_bottom)
      }
    }

  }else if(visium_version == "V1"){
    tissue_positions <- read.csv(tissue_positions_file, header = FALSE)
    colnames(tissue_positions) <- tpos_colnames

    if(multisample){
      if(!is.na(tissue_positions_file_left)){
        tissue_positions_left <- read.csv(tissue_positions_file_left, header = FALSE)
        colnames(tissue_positions_left) <- tpos_colnames
      }
      if(!is.na(tissue_positions_file_right)){
        tissue_positions_right <- read.csv(tissue_positions_file_right, header = FALSE)
        colnames(tissue_positions_right) <- tpos_colnames
      }
      if(!is.na(tissue_positions_file_top)){
        tissue_positions_top <- read.csv(tissue_positions_file_top, header = FALSE)
        colnames(tissue_positions_top) <- tpos_colnames
      }
      if(!is.na(tissue_positions_file_bottom)){
        tissue_positions_bottom <- read.csv(tissue_positions_file_bottom, header = FALSE)
        colnames(tissue_positions_bottom) <- tpos_colnames
      }
    }

  }else{
    stop("'visium_version' must be 'V1' or 'V2'")
  }

  # get scale factor
  hires_scale <- visium_object@images[[1]]@scale.factors$hires

  if(multisample){
    if(!is.na(image_file_left)){
      hires_scale_left <- visium_object_left@images[[1]]@scale.factors$hires
    }
    if(!is.na(image_file_right)){
      hires_scale_right <- visium_object_right@images[[1]]@scale.factors$hires
    }
    if(!is.na(image_file_top)){
      hires_scale_top <- visium_object_top@images[[1]]@scale.factors$hires
    }
    if(!is.na(image_file_bottom)){
      hires_scale_bottom <- visium_object_bottom@images[[1]]@scale.factors$hires
    }
  }

  # get coordinates
  coords <- as.data.frame(Seurat::GetTissueCoordinates(visium_object))
  # rotate coordinates
  coords <- rotate_coordinates(coords)
  # scale coordinates
  coords <- scale_coordinates(coords, tissue_positions)
  # Match to list of clone barcodes
  coords <- match_clone_barcodes(coordinates_df_scaled = coords,
                                 clones_df = clone_df,
                                 clone_group_name = clone_group_column,
                                 coordinate_barcode_name = "cell",
                                 clone_barcode_name = clone_barcode_column)

  if(multisample){

    if(suppressWarnings(!is.na(visium_object_left))){
      # get coordinates
      coords_left <- as.data.frame(Seurat::GetTissueCoordinates(visium_object_left))
      # rotate coordinates
      coords_left <- rotate_coordinates(coords_left)
      # scale coordinates
      coords_left <- scale_coordinates(coords_left, tissue_positions_left)

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
    if(suppressWarnings(!is.na(visium_object_right))){
      # get coordinates
      coords_right <- as.data.frame(Seurat::GetTissueCoordinates(visium_object_right))
      # rotate coordinates
      coords_right <- rotate_coordinates(coords_right)
      # scale coordinates
      coords_right <- scale_coordinates(coords_right, tissue_positions_right)

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
    if(suppressWarnings(!is.na(visium_object_top))){
      # get coordinates
      coords_top <- as.data.frame(Seurat::GetTissueCoordinates(visium_object_top))
      # rotate coordinates
      coords_top <- rotate_coordinates(coords_top)
      # scale coordinates
      coords_top <- scale_coordinates(coords_top, tissue_positions_top)

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
    if(suppressWarnings(!is.na(visium_object_bottom))){
      # get coordinates
      coords_bottom <- as.data.frame(Seurat::GetTissueCoordinates(visium_object_bottom))
      # rotate coordinates
      coords_bottom <- rotate_coordinates(coords_bottom)
      # scale coordinates
      coords_bottom <- scale_coordinates(coords_bottom, tissue_positions_bottom)

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

  if(shared_clones){

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

  }else{
    # TODO: add code here to calculate centroids together as one object
    # merge coords and calculate centroids
    coords_df <- coords

    if(multisample){
      if(!all(is.na(coords_left))){
        coords_df <- rbind(coords_df, coords_left)
      }
      if(!all(is.na(coords_right))){
        coords_df <- rbind(coords_df, coords_right)
      }
      if(!all(is.na(coords_top))){
        coords_df <- rbind(coords_df, coords_top)
      }
      if(!all(is.na(coords_bottom))){
        coords_df <- rbind(coords_df, coords_bottom)
      }

      # NA to other Newick data frames as we have a merged one
      newick_df_left <- NA
      newick_df_right <- NA
      newick_df_top <- NA
      newick_df_bottom <- NA

      newick_tree_df <- calculate_centroids(newick_df = newick_tree_df,
                                            coordinates_df_scaled = coords_df)

      # NA to side segments as we have a joint one
      segments_left <- NA
      segments_right <- NA
      segments_top <- NA
      segments_bottom <- NA

      segments <- create_segments(newick_tree_df)

      coords <- coords_df
    }
  }



  ############# Plotting functions

  # get raster image
  raster_image <- img_get_raster(image_file)

  # crop image to bounding box
  raster_image <- img_crop_raster(raster_image, hires_scale, tissue_positions)
  raster_image <- as.raster(raster_image)

  if(!is.na(image_file_left)){
    raster_image_left <- img_get_raster(image_file_left)
    raster_image_left <- img_crop_raster(raster_image_left, hires_scale_left, tissue_positions_left)
    raster_image_left <- as.raster(raster_image_left)
  }else{
    raster_image_left <- NA
  }
  if(!is.na(image_file_right)){
    raster_image_right <- img_get_raster(image_file_right)
    raster_image_right <- img_crop_raster(raster_image_right, hires_scale_right, tissue_positions_right)
    raster_image_right <- as.raster(raster_image_right)
  }else{
    raster_image_right <- NA
  }
  if(!is.na(image_file_top)){
    raster_image_top <- img_get_raster(image_file_top)
    raster_image_top <- img_crop_raster(raster_image_top, hires_scale_top, tissue_positions_top)
    raster_image_top <- as.raster(raster_image_top)
  }else{
    raster_image_top <- NA
  }
  if(!is.na(image_file_bottom)){
    raster_image_bottom <- img_get_raster(image_file_bottom)
    raster_image_bottom <- img_crop_raster(raster_image_bottom, hires_scale_bottom, tissue_positions_bottom)
    raster_image_bottom <- as.raster(raster_image_bottom)
  }else{
    raster_image_bottom <- NA
  }




  # get colour palette
  if(any(palette == "default")){
    my_palette <- img_name_colours(newick_tree_df)

    if(multisample & shared_clones){
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
    shared_clones = shared_clones,
    plot_connections = plot_connections,
    connections_coords = connections_coords,
    connection_colour = connection_colour,
    connection_width = connection_width
  )

  return(plot)
}
