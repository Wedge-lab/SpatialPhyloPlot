# Testing imaging related function

# Plotting functions

# get raster image
test_that("Test function fails if input data is inappropritate or missing but runs if correct.",{
  expect_error(img_get_raster("not_an_image.txt"))
  expect_no_error(img_get_raster(image_file_path))
})

# name colours (optional)

test_that("Default palette works if no option supplied.",{
  expect_no_error(SpatialPhyloPlot(visium_object = demo_visium,
                                   newick_file = newick_file_path,
                                   image_file = image_file_path,
                                   clone_df = demo_clones,
                                   clone_group_column = "group",
                                   clone_barcode_column = "Barcodes2"))
})
test_that("Custom palette works if supplied as named list.",{
  expect_no_error(SpatialPhyloPlot(visium_object = demo_visium,
                                   newick_file = newick_file_path,
                                   image_file = image_file_path,
                                   clone_df = demo_clones,
                                   palette = test_palette,
                                   clone_group_column = "group",
                                   clone_barcode_column = "Barcodes2"))
})
test_that("Error if clone name not in list. ",{
  expect_error(SpatialPhyloPlot(visium_object = demo_visium,
                                newick_file = newick_file_path,
                                image_file = image_file_path,
                                clone_df = demo_clones,
                                palette = wrong_palette,
                                clone_group_column = "group",
                                clone_barcode_column = "Barcodes1"))
})

# plot

# gif version
