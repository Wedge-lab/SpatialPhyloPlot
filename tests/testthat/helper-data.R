### Demo data
# library(Seurat)
# library(SeuratData)
# remotes::install_github("rstudio/ggcheck")
library(ggcheck) # for checking plot layers and geoms

# withr::with_package("SeuratData", {demo_visium <- LoadData("stxBrain", type = "anterior1")})
demo_visium <- SeuratData::LoadData("stxBrain", type = "anterior1")
# demo_false_colours

# demo_newick_file
newick_file_path <- "./testdata/demo.new"

# demo clone barcode data frame
clone_barcodes <- as.data.frame(GetTissueCoordinates(demo_visium))
colnames(clone_barcodes) <- c("x","y","Barcodes1")
clone_barcodes$group <- c(rep("A",1000),rep("B",1000),rep("C",606), rep("D",90))
demo_clones <- clone_barcodes[,c("Barcodes1","group")]
colnames(demo_clones) <- c("Barcodes2","group")

# mismatched barcodes df
clone_barcodes_mismatched <- demo_clones
clone_barcodes_mismatched$Barcodes2 <- gsub("A","X",clone_barcodes_mismatched$Barcodes2)
clone_barcodes_mismatched$Barcodes2 <- gsub("C","Y",clone_barcodes_mismatched$Barcodes2)


# scaled barcode df
clone_barcodes_rotate <- clone_barcodes
clone_barcodes_rotate$new_x <- as.numeric(clone_barcodes_rotate$y)
clone_barcodes_rotate$new_y <- as.numeric(-clone_barcodes_rotate$x)

# scaling barcodes - this data isn't strictly speaking correct as the real approach uses
# the image dimensions to get the maximum coordinates, which is more correct, as this
# just uses the available barcodes/spots
clone_barcodes_scaled <- clone_barcodes_rotate
clone_barcodes_scaled$new_x_scaled <- clone_barcodes_rotate$new_x/(max(clone_barcodes_rotate$new_x))
clone_barcodes_scaled$new_y_scaled <- (clone_barcodes_rotate$new_y-min(clone_barcodes_rotate$new_y))/(-min(clone_barcodes_rotate$new_y))

# newick in wrong file format

# ? demo newickdf?

# imgae file
# img.url <- 'http://cf.10xgenomics.com/samples/spatial-exp/1.0.0/V1_Mouse_Brain_Sagittal_Anterior/V1_Mouse_Brain_Sagittal_Anterior_spatial.tar.gz'
# curl::curl_download(url = img.url, destfile = basename(path = img.url))
# untar(tarfile = basename(path = img.url))
# move hires png to test data
image_file_path <- "./testdata/tissue_hires_image.png"

# Make a colour palette
test_palette <- c("blue","grey","black","red")
names(test_palette) <- c("A","B","C","D")

wrong_palette <- c("blue","grey","black","red")
names(wrong_palette) <- c("A","B","C","J")

# Multinewick file (partial columns)
test_multinewick <- as.data.frame(matrix(nrow = 5, ncol = 4))
colnames(test_multinewick) <- c("colour","Origin","centroid_x","centroid_y")
test_multinewick$colour <- rep("Clone_X",5)
test_multinewick$Origin <- c("Centre","Left","Right","Top","Bottom")
test_multinewick$centroid_x <- rep(0,5)
test_multinewick$centroid_y <- rep(0,5)

