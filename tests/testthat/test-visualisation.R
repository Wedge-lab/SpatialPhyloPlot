# Testing imaging related function

# Plotting functions

# get raster image
test_that("Test function fails if input data is inappropritate or missing but runs if correct.",{
  expect_error(img_get_raster("not_an_image.txt"))
  expect_no_error(img_get_raster(image_file_path))
})

# name colours (optional)

test_that("Default palette works if no option supplied.",{
  expect_no_error(SpatialPhyloPlot(newick_file = newick_file_path,
                                   image_file = image_file_path,
                                   scale_factors = scale_factor_path,
                                   clone_df = demo_clones,
                                   clone_group_column = "group",
                                   clone_barcode_column = "Barcodes2",
                                   visium_version = "V1",
                                   tissue_positions_file = tissue_positions_path
                                   ))
})
test_that("Custom palette works if supplied as named list.",{
  expect_no_error(SpatialPhyloPlot(newick_file = newick_file_path,
                                   image_file = image_file_path,
                                   scale_factors = scale_factor_path,
                                   clone_df = demo_clones,
                                   palette = test_palette,
                                   clone_group_column = "group",
                                   clone_barcode_column = "Barcodes2",
                                   visium_version = "V1",
                                   tissue_positions_file = tissue_positions_path))
})
test_that("Error if clone name not in list. ",{
  expect_error(SpatialPhyloPlot(newick_file = newick_file_path,
                                image_file = image_file_path,
                                scale_factors = scale_factor_path,
                                clone_df = demo_clones,
                                palette = wrong_palette,
                                clone_group_column = "group",
                                clone_barcode_column = "Barcodes1",
                                visium_version = "V1",
                                tissue_positions_file = tissue_positions_path))
})

# plot

# test plot coordinates

# test plot elements

test_that("Correct layers produced if you provide the right input. ", {
  output <- SpatialPhyloPlot(newick_file = newick_file_path,
                             image_file = image_file_path,
                             scale_factors = scale_factor_path,
                             clone_df = clone_barcodes,
                             clone_group_column = "group",
                             clone_barcode_column = "Barcodes1",
                             visium_version = "V1",
                             tissue_positions_file = tissue_positions_path)
  expect_equal(sum(get_geoms(output) == "rasterann"),1)
  expect_equal(sum(get_geoms(output) == "point"),2)

})

test_that("Correct layers produced if you provide the right input in multisample. ", {
  output <- SpatialPhyloPlot(newick_file = newick_file_path,
                             image_file = image_file_path,
                             scale_factors = scale_factor_path,
                             clone_df = clone_barcodes,
                             clone_group_column = "group",
                             clone_barcode_column = "Barcodes1",
                             multisample = TRUE,
                             image_file_left = image_file_path,
                             clone_df_left = clone_barcodes,
                             visium_version = "V1",
                             tissue_positions_file = tissue_positions_path,
                             tissue_positions_file_left = tissue_positions_path,
                             scale_factors_left = scale_factor_path,
                             shared_clones = TRUE)

  expect_equal(sum(get_geoms(output) == "rasterann"),2)
  expect_equal(sum(get_geoms(output) == "point"),4)

})

test_that("Correct layers produced if you provide the right input in multisample with hull. ", {
  output <- SpatialPhyloPlot(newick_file = newick_file_path,
                             image_file = image_file_path,
                             scale_factors = scale_factor_path,
                             clone_df = clone_barcodes,
                             clone_group_column = "group",
                             clone_barcode_column = "Barcodes1",
                             multisample = TRUE,
                             image_file_left = image_file_path,
                             clone_df_left = clone_barcodes,
                             plot_polygon = TRUE,
                             visium_version = "V1",
                             tissue_positions_file = tissue_positions_path,
                             tissue_positions_file_left = tissue_positions_path,
                             scale_factors_left = scale_factor_path,
                             shared_clones = TRUE)

  expect_equal(sum(get_geoms(output) == "rasterann"),2)
  expect_equal(sum(get_geoms(output) == "point"),4)
  expect_equal(sum(get_geoms(output) == "markhull"),2)

  output <- SpatialPhyloPlot(newick_file = newick_file_path,
                             image_file = image_file_path,
                             scale_factors = scale_factor_path,
                             clone_df = clone_barcodes,
                             clone_group_column = "group",
                             clone_barcode_column = "Barcodes1",
                             multisample = TRUE,
                             image_file_right = image_file_path,
                             clone_df_right = clone_barcodes,
                             plot_polygon = TRUE,
                             visium_version = "V1",
                             tissue_positions_file = tissue_positions_path,
                             tissue_positions_file_right = tissue_positions_path,
                             scale_factors_right = scale_factor_path,
                             shared_clones = TRUE)

  expect_equal(sum(get_geoms(output) == "rasterann"),2)
  expect_equal(sum(get_geoms(output) == "point"),4)
  expect_equal(sum(get_geoms(output) == "markhull"),2)

  output <- SpatialPhyloPlot(newick_file = newick_file_path,
                             image_file = image_file_path,
                             scale_factors = scale_factor_path,
                             clone_df = clone_barcodes,
                             clone_group_column = "group",
                             clone_barcode_column = "Barcodes1",
                             multisample = TRUE,
                             image_file_top = image_file_path,
                             clone_df_top = clone_barcodes,
                             plot_polygon = TRUE,
                             visium_version = "V1",
                             tissue_positions_file = tissue_positions_path,
                             tissue_positions_file_top = tissue_positions_path,
                             scale_factors_top = scale_factor_path,
                             shared_clones = TRUE)

  expect_equal(sum(get_geoms(output) == "rasterann"),2)
  expect_equal(sum(get_geoms(output) == "point"),4)
  expect_equal(sum(get_geoms(output) == "markhull"),2)

  output <- SpatialPhyloPlot(newick_file = newick_file_path,
                             image_file = image_file_path,
                             scale_factors = scale_factor_path,
                             clone_df = clone_barcodes,
                             clone_group_column = "group",
                             clone_barcode_column = "Barcodes1",
                             multisample = TRUE,
                             image_file_bottom = image_file_path,
                             clone_df_bottom = clone_barcodes,
                             plot_polygon = TRUE,
                             visium_version = "V1",
                             tissue_positions_file = tissue_positions_path,
                             tissue_positions_file_bottom = tissue_positions_path,
                             scale_factors_bottom = scale_factor_path,
                             shared_clones = TRUE)

  expect_equal(sum(get_geoms(output) == "rasterann"),2)
  expect_equal(sum(get_geoms(output) == "point"),4)
  expect_equal(sum(get_geoms(output) == "markhull"),2)

})
test_that("Correct layers produced if you provide the right input in multisample with connections ", {
  output <- SpatialPhyloPlot(newick_file = newick_file_path,
                             image_file = image_file_path,
                             scale_factors = scale_factor_path,
                             clone_df = clone_barcodes,
                             clone_group_column = "group",
                             clone_barcode_column = "Barcodes1",
                             multisample = TRUE,
                             image_file_left = image_file_path,
                             clone_df_left = clone_barcodes,
                             plot_connections = TRUE,
                             visium_version = "V1",
                             tissue_positions_file = tissue_positions_path,
                             tissue_positions_file_left = tissue_positions_path,
                             scale_factors_left = scale_factor_path,
                             shared_clones = TRUE)

  expect_equal(sum(get_geoms(output) == "rasterann"),2)
  expect_equal(sum(get_geoms(output) == "point"),4)
  expect_equal(sum(get_geoms(output) == "segment"),3)

  output <- SpatialPhyloPlot(newick_file = newick_file_path,
                             image_file = image_file_path,
                             scale_factors = scale_factor_path,
                             clone_df = clone_barcodes,
                             clone_group_column = "group",
                             clone_barcode_column = "Barcodes1",
                             multisample = TRUE,
                             image_file_right = image_file_path,
                             clone_df_right = clone_barcodes,
                             plot_connections = TRUE,
                             visium_version = "V1",
                             tissue_positions_file = tissue_positions_path,
                             tissue_positions_file_right = tissue_positions_path,
                             scale_factors_right = scale_factor_path,
                             shared_clones = TRUE)

  expect_equal(sum(get_geoms(output) == "rasterann"),2)
  expect_equal(sum(get_geoms(output) == "point"),4)
  expect_equal(sum(get_geoms(output) == "segment"),3)

  output <- SpatialPhyloPlot(newick_file = newick_file_path,
                             image_file = image_file_path,
                             scale_factors = scale_factor_path,
                             clone_df = clone_barcodes,
                             clone_group_column = "group",
                             clone_barcode_column = "Barcodes1",
                             multisample = TRUE,
                             image_file_top = image_file_path,
                             clone_df_top = clone_barcodes,
                             plot_connections = TRUE,
                             visium_version = "V1",
                             tissue_positions_file = tissue_positions_path,
                             tissue_positions_file_top = tissue_positions_path,
                             scale_factors_top = scale_factor_path,
                             shared_clones = TRUE)

  expect_equal(sum(get_geoms(output) == "rasterann"),2)
  expect_equal(sum(get_geoms(output) == "point"),4)
  expect_equal(sum(get_geoms(output) == "segment"),3)

  output <- SpatialPhyloPlot(newick_file = newick_file_path,
                             image_file = image_file_path,
                             scale_factors = scale_factor_path,
                             clone_df = clone_barcodes,
                             clone_group_column = "group",
                             clone_barcode_column = "Barcodes1",
                             multisample = TRUE,
                             image_file_bottom = image_file_path,
                             clone_df_bottom = clone_barcodes,
                             plot_connections = TRUE,
                             visium_version = "V1",
                             tissue_positions_file = tissue_positions_path,
                             tissue_positions_file_bottom = tissue_positions_path,
                             scale_factors_bottom = scale_factor_path,
                             shared_clones = TRUE)

  expect_equal(sum(get_geoms(output) == "rasterann"),2)
  expect_equal(sum(get_geoms(output) == "point"),4)
  expect_equal(sum(get_geoms(output) == "segment"),3)

})


# gif version
