% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preparation.R
\name{connect_multisample}
\alias{connect_multisample}
\title{Create connections between clones in multiple samples}
\usage{
connect_multisample(combined_newick_df)
}
\arguments{
\item{combined_newick_df}{A concatenated Newick `data.frame` with the coordinates of each clone for each sample.}
}
\value{
A `data.frame` with the start and end coordinates from the main plot to the others.
}
\description{
Create connections between clones in multiple samples
}
