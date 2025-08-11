---
title: 'SpatialPhyloPlot: An R Package for Automated Plotting of Phylogenetic Trees
  on Spatial Transcriptomics Data'
tags:
- R
- Spatial transcriptomics
- Phylogenetic trees
- Data visualisation
- Tumour evolution
date: "11 August 2025"
output: pdf_document
authors:
  - name: G. Maria Jakobsdottir
    orcid: "0000-0001-8225-6541"
    corresponding: true
    affiliation: 1, 2
  - name: Arfa Mehmood
    affiliation: 1, 2
  - name: Miaomiao Gao
    affiliation: 1
  - name: David C. Wedge
    orcid: "0000-0002-7572-3196"
    corresponding: true
    affiliation: 1, 2
bibliography: paper.bib
affiliations:
  - name: Division of Cancer Sciences, The University of Manchester, 555 Wilmslow Road, Manchester, M20 4GJ, UK
    index: 1
  - name: Christie Hospital, The Christie NHS Foundation Trust, Manchester Academic Health Science Centre, Manchester, M20 4BX, UK
    index: 2
---

# Summary
Data visualisation and presentation are often key steps in data interpretation and dissemination of results.
The advent of spatially resolved transcriptomics has provided a novel visual aspect to traditional sequencing technologies, allowing behaviours such as tumour progression to be visualised in their combined spatial and omic context.
Recent advances, such as the inference of copy number profiles from spatial transcriptomics [@erickson2022spatially], and calculation of phylogenetic trees from copy number data [@kaufmann2022medicc2], have provided the opportunity to visualise the spread of clones throughout a tissue.
However, no tools are currently available to plot phylogenetic trees in their spatial context.

We present `SpatialPhyloPlot`, an `R` package which automates the process of plotting phylogenetic trees on top of their associated tissue regions.
`SpatialPhyloPlot` builds on popular `R` based tools for visualisation and analysis of data such as `ggplot2` [@wickham2016ggplot2], resulting in spatial phylogenetic plots that are readily customisable and easily integrated with existing user workflows.
`SpatialPhyloPlot` is currently compatible with Visium 10X spatial transcriptomics data and phylogenetic trees provided in the Newick file format, and it is freely available through github.

# Statement of Need

To our knowledge, this is the first package offering an automated method for plotting phylogenetic trees derived from spatial transcriptomic data in their spatial context.
This allows direct visualisation of the spread of clones through a tissue and the evolutionary relationship between these clones.

# Functionality

`SpatialPhyloPlot` is available to download and install through github at \url{https://github.com/Wedge-lab/SpatialPhyloPlot}.
The latest package release, or specific branches under development, can be installed from github, along with any dependencies, using the `devtools` `R` package [@wickham2022devtools]:

```
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("Wedge-lab/SpatialPhyloPlot", dependencies = TRUE)
```

`SpatialPhyloPlot` imports a number of common `R` packages to assist with plotting and data wrangling.
Firstly, the plots produced by `SpatialPhyloPlot` utilise the `ggplot2` [@wickham2016ggplot2] plotting framework, making them easily customisable, modifiable, and simple to integrate with other `ggplot2` based plot elements, including polygons or hulls created with the `ggplot2` extension `ggforce` [@pedersen2024ggforce].
Plotting and appropriate cropping of tissue images is facilitated by the `magick` `R` package [@ooms2025magick].
The `alakazam` [@gupta2015alakazam], `ape` [@paradis2019ape], and `igraph` [@csardi2006igraph; @csardi2025igraph] `R` packages are utilised internally to facilitate conversion of Newick format phylogenetic trees to `R` `data.frame` structures for easy manipulation and plotting.
Finally, distinct colour palettes for differentiating clones are provided through the `pals` `R` package [@wright2025pals].

`SpatialPhyloPlot` is currently compatible with Visium spatial transcriptomic data produced by 10X Genomics and phylogentic trees presented as `.new` files in the Newick format.
It has one main plotting function with 6 essential inputs, which are either produced by the 10X Genomics SpaceRanger pipeline, or software that produces phylogenetic trees, such as MEDICC2 [@kaufmann2022medicc2] (\autoref{fig:overview}A).

Essential parameters:

  - `newick_file`: Path to a Newick format phylogenetic tree file with `.new` file extension.
  - `image_file`: Path to the `tissue_hires_image.png` file produced by SpaceRanger.
  - `scale_factors`: Path to the `scalefactors_json.json` file produced by SpaceRanger.
  - `clone_df`: A `data.frame` indicating which clone in the phylogenetic tree corresponds to each barcode in the Visium assay.
  - `visium_version`: The version of the Visium assay or SpaceRanger used (currently "V1" or "V2").
  - `tissue_positions_file`: Path to the `tissue_positions.csv` or `tissue_positions_list.csv` file produced by SpaceRanger.

![A. A summary of the essential inputs required to run `SpatialPhyloPlot`. B. The main plot elements of `SpatialPhyloPlot`, and the main parameters used to customise them. C. Examples of how plot elements can be combined. D. The coordinate system and layout of a single central `SpatialPhyloPlot` plot. E. The coordinate system and layout of multiple samples in a `SpatialPhyloPlot` plot. \label{fig:overview}](paper_figs/Overview_figure.png)

Internally, `SpatialPhyloPlot` links the clones to the location of the points in the Visium assay and then calculates a central point, or "centroid", for each of the clones based on the average x and y location of the points within it.
The centroids of internal nodes in the tree are calculated as the average of the centroids to which it links.
This process is repeated until centroids have been calculated for all connecting nodes.
`SpatialPhyloPlot` can then plot the individual Visium points, centroids, and segments linking the centroids with the option of representing internal nodes, on top of the tissue image, cropped to the appropriate size.

The core plot elements of the figures produced by `SpatialPhyloPlot` are individual Visium points coloured by clone, polygons encapsulating clones, and centroids with connecting segments (\autoref{fig:overview}B).
These plot elements can be mixed and matched by the user depending on what is appropriate and aesthetically pleasing in each case (\autoref{fig:overview}C).
The coordinates of each plot and corresponding spots are scaled to x and y intervals of $[0,1]$, to create a consistent coordinate system to allow easy manipulation and annotation by the user (\autoref{fig:overview}D).
If multiple samples are presented together, for example a primary and metastatic sample from the same patient or multiple primary samples taken from the same patient, this coordinate system is extended by one unit in the x or y direction depending on whether the second sample is plotted to the left, right, above, or below the primary plot (\autoref{fig:overview}E).

![A. Code used to produce the plot in B. B. A phylogenetic tree plotted on top of the 10X Genomics Human Breast Cancer  (Block A Section 1) Visium V1 spatial transcriptomic data set. C. The phylogenetic tree produced by MEDICC2 and used as input to the code in A to produce the plot in B. \label{fig:breast_cancer}](paper_figs/breast_cancer.png)

We provide `SpatialPhyloPlot` users with an [extensive guide](https://wedge-lab.github.io/SpatialPhyloPlot/articles/exploring_plotting_options_in_SpatialPhyloPlot.html) that works through the different ways to represent the data.
Users can then try different approaches and select the representation that is the most appropriate for their data.
We also illustrate `SpatialPhyloPlot`'s performance here using a real example.
The Human Breast Cancer Visium V1 spatial transcriptomic data set processed using Space Ranger [@spaceranger] was downloaded from 10X Genomics [@visium2020breast].
We then ran InferCNV [@trinityinfercnv} using the SpatialInferCNV [@erickson2022spatially] workflow to identify tumour clones and created a data frame matching them to each Visium barcode.
Copy number based phylogenetic trees were built using MEDICC2 [@kaufmann2022medicc2].
We then passed the Visium image (`image_primary`), Visium tissue positions (`tissue_pos_file`), Visium scale factors (`scale_file`) barcode clone annotation data frame (`clones_primary`), and Newick format phylogenetic tree file (`newfile`) to the `SpatialPhyloPlot` function to produce \autoref{fig:breast_cancer}.


# Use Cases

`SpatialPhyloPlot` is a new `R` package and thus has not yet been used in any scholarly publications.
However, we currently have at least two manuscripts in preparation featuring figures produced using `SpatialPhyloPlot`.
In the first publication, `SpatialPhyloPlot` is used to map the relationship between primary colorectal cancers and their associated peritoneal metastases [@pmet2026gao].
In the latter, `SpatialPhyloPlot` is being used to map clones across a set of 12 heterogeneous prostate cancer samples [@prostate2026mehmood].

# Acknowledgements

This research was supported by the NIHR Manchester Biomedical Research Centre (NIHR203308).
The views expressed are those of the author(s) and not necessarily those of the NIHR or the Department of Health and Social Care.

# References
