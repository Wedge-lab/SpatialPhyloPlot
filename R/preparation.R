# Functions for wrangling phylogenetic trees and data

#' Rotate Image Coordinates
#'
#' @param coordinates_df A `data.frame` of spot coordinates extracted from the `Seurat` visium 10X data.
#'
#' @returns A `data.frame` with rotated Visium 10X spot coordinates.
#'
rotate_coordinates <- function(coordinates_df) {

  coordinates_df$new_x <- coordinates_df$y
  coordinates_df$new_y <- -as.numeric(coordinates_df$x)

  coordinates_df$new_x <- as.numeric(coordinates_df$new_x)
  coordinates_df$new_y <- as.numeric(coordinates_df$new_y)

  return(coordinates_df)
}

#' Scale Coordinates
#'
#' @param tissue_positions A `data.frame` with tissue positions for each spot.
#' @param coordinates_df A `data.frame` with rotated spot coordinates.
#'
#' @returns A `data.frame` with Visium 10X coordinates scaled from 0-1.
#'
scale_coordinates <- function(coordinates_df, tissue_positions){

  xmax_coord <- max(tissue_positions$pxl_col_in_fullres)
  ymax_coord <- max(tissue_positions$pxl_row_in_fullres)

  xmin_coord <- min(tissue_positions$pxl_col_in_fullres)
  ymin_coord <- min(tissue_positions$pxl_row_in_fullres)

  coordinates_df$new_x_scaled <- (coordinates_df$new_x - xmin_coord)/(xmax_coord - xmin_coord)
  coordinates_df$new_y_scaled <- (coordinates_df$new_y + ymax_coord)/(ymax_coord - ymin_coord)

  return(coordinates_df)

}

#' Match clone names and spot barcodes
#'
#' @param coordinates_df_scaled A `data.frame` with scaled spot coordinates.
#' @param clones_df A `data.frame` indicating which barcodes correspond to which clones.
#' @param clone_group_name The name of the column that contains the clone grouping information.
#' @param coordinate_barcode_name The name of the column that contains the barcodes in the coordinates `data.frame`.
#' @param clone_barcode_name The name of the column that contains the barcodes in the clones `data.frame`.
#'
#' @returns A `data.frame` with Visium 10X barcodes and their coordinates matched to clone labels.
#'
match_clone_barcodes <- function(coordinates_df_scaled,
                                 clones_df,
                                 clone_group_name,
                                 coordinate_barcode_name = "cell",
                                 clone_barcode_name = "Barcode") {

  if(!any(coordinates_df_scaled[,coordinate_barcode_name] %in% clones_df[,clone_barcode_name])){
    stop("Clone barcodes do not match Visium coordinate barcodes. ")
  }
  if(!all(coordinates_df_scaled[,coordinate_barcode_name] %in% clones_df[,clone_barcode_name]) |
     !all(clones_df[,clone_barcode_name] %in% coordinates_df_scaled[,coordinate_barcode_name])){
    warning("Some barcodes are not present in all inputs. Only barcodes present in both objects will be used.")
  }

  colnames(clones_df)[which(colnames(clones_df) == clone_group_name)] <- "Clone"

  coordinates_df <- merge(coordinates_df_scaled,
                          clones_df,
                          by.x = coordinate_barcode_name,
                          by.y = clone_barcode_name)

  return(coordinates_df)

}

#' Newick tree to graph data frame
#'
#' @param newick_file The path to a `.new` file containing the phylogenetic tree in Newick format.
#'
#' @returns A `data.frame` containing the phylogenetic tree information re-formatted as a `data.frame`.
#'
newick_to_graph_df <- function(newick_file) {

  if(!file.exists(newick_file)){
    stop("Newick tree file cannot be found.")
  }
  if(tools::file_ext(newick_file) != "new"){
    stop("Phylogenetic tree must be provided in the Newick format as a file with '.new' extension")
  }

  ## TODO: add support for tree provided as a data frame later

  newick_tree <- ape::read.tree(newick_file)
  newick_graph <- alakazam::phyloToGraph(newick_tree)
  newick_df <- igraph::as_data_frame(newick_graph)

  return(newick_df)
}

#' Calculate centroids
#'
#' @param newick_df A `data.frame` containing the reformatted Newick phylogenetic tree information.
#' @param coordinates_df_scaled A `data.frame` containing scaled spot coordinates and the barcodes that they correspond to.
#'
#' @returns A `data.frame` with centroids calculated for each clone in the phylogenetic tree.
#'
calculate_centroids <- function(newick_df, coordinates_df_scaled) {
  # Try new way of formatting tree
  newick_df$centroid_x <- NA
  newick_df$centroid_y <- NA

  for(i in 1:nrow(newick_df)){
    spot_sub <- coordinates_df_scaled[coordinates_df_scaled[,"Clone"] == newick_df$to[i],]#subset(coordinates_df_scaled, clone_group_column == newick_df$to[i])
    if(nrow(spot_sub) > 0){
      newick_df$centroid_x[i] <- mean(spot_sub$new_x_scaled)
      newick_df$centroid_y[i] <- mean(spot_sub$new_y_scaled)
    }
  }

  # calculate a centroid of centroids for in between points
  for(i in 1:nrow(newick_df)){

    if(is.na(newick_df$centroid_x[i])){
      parents <- subset(newick_df, to == newick_df$from[i])
      relations <- subset(newick_df, from == newick_df$to[i])

      family <- rbind(parents, relations)

      newick_df$centroid_x[i] <- mean(na.omit(family$centroid_x))
      newick_df$centroid_y[i] <- mean(na.omit(family$centroid_y))

    }
  }

  newick_df$colour <- ifelse(newick_df$to %in% coordinates_df_scaled[,"Clone"], yes = newick_df$to, no = NA)

  return(newick_df)

}

#' Create phylogenetic tree segments
#'
#' @param newick_df A `data.frame` containing reformatted Newick phylogenetic tree information.
#'
#' @returns A `data.frame` with the coordinates of each segment on the phylogenetic tree.
#'
create_segments <- function(newick_df) {
  from_to <- newick_df[,c("to","from","centroid_x","centroid_y")]
  colnames(from_to) <- c("to","from","to_x","to_y")
  from_to <- merge(from_to, newick_df[,c("to","centroid_x","centroid_y")], by.x = "from", by.y = "to", all = TRUE)
  colnames(from_to) <- c("to","from","to_x","to_y","from_x","from_y")

  return(from_to)

}

#' Create connections between clones in multiple samples
#'
#' @param combined_newick_df A concatenated Newick `data.frame` with the coordinates of each clone for each sample.
#'
#' @returns A `data.frame` with the start and end coordinates from the main plot to the others.
#'
connect_multisample <- function(combined_newick_df) {

  clones <- unique(na.omit(combined_newick_df$colour))

  connections_coordinates <- as.data.frame(matrix(nrow = length(clones), ncol = 11))
  colnames(connections_coordinates) <- c("Clone","centre_x","centre_y",
                                         "left_x","left_y",
                                         "right_x","right_y",
                                         "top_x","top_y",
                                         "bottom_x","bottom_y")
  connections_coordinates$Clone <- clones

  for(i in 1:nrow(connections_coordinates)){

    clone_info <- subset(combined_newick_df, colour == connections_coordinates$Clone[i])

    if(nrow(subset(clone_info, Origin == "Centre")) > 0){
      connections_coordinates$centre_x[i] <- subset(clone_info, Origin == "Centre")$centroid_x
      connections_coordinates$centre_y[i] <- subset(clone_info, Origin == "Centre")$centroid_y
    }

    if(nrow(subset(clone_info, Origin == "Left")) > 0){
      connections_coordinates$left_x[i] <- subset(clone_info, Origin == "Left")$centroid_x
      connections_coordinates$left_y[i] <- subset(clone_info, Origin == "Left")$centroid_y
    }

    if(nrow(subset(clone_info, Origin == "Right")) > 0){
      connections_coordinates$right_x[i] <- subset(clone_info, Origin == "Right")$centroid_x
      connections_coordinates$right_y[i] <- subset(clone_info, Origin == "Right")$centroid_y
    }

    if(nrow(subset(clone_info, Origin == "Top")) > 0){
      connections_coordinates$top_x[i] <- subset(clone_info, Origin == "Top")$centroid_x
      connections_coordinates$top_y[i] <- subset(clone_info, Origin == "Top")$centroid_y
    }

    if(nrow(subset(clone_info, Origin == "Bottom")) > 0){
      connections_coordinates$bottom_x[i] <- subset(clone_info, Origin == "Bottom")$centroid_x
      connections_coordinates$bottom_y[i] <- subset(clone_info, Origin == "Bottom")$centroid_y
    }

  }

  connections_coordinates <- subset(connections_coordinates, !is.na(centre_x))

  return(connections_coordinates)

}
