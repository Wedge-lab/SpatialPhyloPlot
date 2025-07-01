# testing data wrangling related functions

### Testing Inputs
test_that("Runs if input is Seurat object, otherwise fails", {
  expect_error(SpatialPhyloPlot(demo_visium,
                                newick_file_path,
                                clone_barcodes,
                                clone_barcode_column = "Barcodes1"))
})
test_that("Clone barcode input is data frame", {
  expect_error(SpatialPhyloPlot(demo_visium,
               newick_file_path,
               as.matrix(clone_barcodes),
               clone_barcode_column = "Barcodes1"))
})
test_that("Clone barcode input has >0 lines", {

})
###

# get coordinates
# get full image coordinates
test_that("All coordinates are in the range 0-1",{
  out_coords <- scale_coordinates(clone_barcodes_rotate, tissue_positions_data)

  expect_false(any(out_coords$xmin < 0))
  expect_false(any(out_coords$xmin > 1))
  expect_false(any(out_coords$xmax < 0))
  expect_false(any(out_coords$xmax > 1))
  expect_false(any(out_coords$ymin < 0))
  expect_false(any(out_coords$ymin > 1))
  expect_false(any(out_coords$ymax < 0))
  expect_false(any(out_coords$ymax > 1))

})

# rotate coordinates

# Match to list of clone barcodes
test_that("Runs if input has >0 coordinates in tissue and fails if coordinates do not match. Also tests custom barcode name", {
  expect_no_error(match_clone_barcodes(coordinates_df_scaled = clone_barcodes_scaled,
                                      clones_df = demo_clones,
                                      clone_group_name = "group",
                                      coordinate_barcode_name = "Barcodes1",
                                      clone_barcode_name = "Barcodes2"))
  expect_error(match_clone_barcodes(coordinates_df_scaled = clone_barcodes_scaled,
                                    clones_df = clone_barcodes_mismatched,
                                    clone_group_name = "group",
                                    coordinate_barcode_name = "Barcodes1",
                                    clone_barcode_name = "Barcodes2"))
})
test_that("There is a warning if there isn't complete overlap between clone barcodes?", {
  expect_warning(match_clone_barcodes(clone_barcodes_scaled,
                                      demo_clones[-c(1:10),],
                                      clone_group_name = "group",
                                      coordinate_barcode_name = "Barcodes1",
                                      clone_barcode_name = "Barcodes2"))
})

# Newick tree to graph df
test_that("Input tree file is .new", {
  expect_error(newick_to_graph_df(image_file_path))
  expect_error(newick_to_graph_df("text"))
  expect_no_error(newick_to_graph_df(newick_file_path))
})


# Calculate centroids

# Create segments
test_that("connect_multisample produces a data.frame with a single row and specific columns", {
  expect_no_error(connect_multisample(test_multinewick))

  output <- connect_multisample(test_multinewick)

  expect_equal(nrow(output),1)
  expect_equal(colnames(output),c("Clone","centre_x","centre_y",
                                  "left_x","left_y",
                                  "right_x","right_y",
                                  "top_x","top_y",
                                  "bottom_x","bottom_y"))
})
