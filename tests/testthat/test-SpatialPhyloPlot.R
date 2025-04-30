# Test inputs and ouptuts
test_that("Error if visium object is not a seurat object. ",{
  expect_error(SpatialPhyloPlot(visium_object = demo_clones,
                                newick_file = newick_file_path,
                                image_file = image_file_path,
                                clone_df = demo_clones,
                                clone_group_column = "group",
                                clone_barcode_column = "Barcodes1"))
})
test_that("Error if clone group name not right or not supplied. ",{
  expect_error(SpatialPhyloPlot(visium_object = demo_clones,
                                newick_file = newick_file_path,
                                image_file = image_file_path,
                                clone_df = demo_clones,
                                clone_group_column = "nonsense",
                                clone_barcode_column = "Barcodes1"))
  expect_error(SpatialPhyloPlot(visium_object = demo_clones,
                                newick_file = newick_file_path,
                                image_file = image_file_path,
                                clone_df = demo_clones,
                                clone_barcode_column = "Barcodes1"))
})
test_that("Error if clone_df is not a data.frame. ",{
  expect_error(SpatialPhyloPlot(visium_object = demo_clones,
                                newick_file = newick_file_path,
                                image_file = image_file_path,
                                clone_df = as.matrix(demo_clones),
                                clone_group_column = "group",
                                clone_barcode_column = "Barcodes1"))
})
