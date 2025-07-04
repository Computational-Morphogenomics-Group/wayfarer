% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multivariate.R
\name{plotMoranPCs}
\alias{plotMoranPCs}
\title{Plot Moran's I of bin projection in PCA space}
\usage{
plotMoranPCs(sfes, npcs = 20, reduction = "PCA")
}
\arguments{
\item{sfes}{A list of SFE objects, with Moran's I computed for bin projections
in the dimension reduction of interest.}

\item{npcs}{Number of components whose variance explained or eigenvalues are
to be plotted.}

\item{reduction}{Name of the dimension reduction of interest}
}
\value{
A \code{ggplot} object
}
\description{
Basically see how spatially structured each principal component is and how it
may relate to bin size.
}
