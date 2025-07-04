% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{plotClusterMedians}
\alias{plotClusterMedians}
\title{Plot cluster median curves}
\usage{
plotClusterMedians(df, cluster_col)
}
\arguments{
\item{df}{Data frame output from either \code{\link{clusterLeeCurves}} or
\code{\link{clusterMoranCurves}}. There must be a column named either
"moran" or "lee" and a column called "side" in addition to the column
specified in \code{cluster_col}.}

\item{cluster_col}{Name of a categorical column in \code{df} to use as
cluster labels.}
}
\value{
A \code{ggplot2} object.
}
\description{
To better juxtapose different cluster patterns in the same plot.
}
