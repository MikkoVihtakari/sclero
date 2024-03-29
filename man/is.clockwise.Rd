% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is.clockwise.R
\name{is.clockwise}
\alias{is.clockwise}
\title{Check whether points for an \code{\link[spatstat]{owin}} are clockwise}
\usage{
is.clockwise(x)
}
\arguments{
\item{x}{a data frame with x coordinates in the first column and y coordinates in the second.}
}
\value{
A logical telling whether the polygon is arranged clockwise.
}
\description{
Check whether points for an \code{\link[spatstat]{owin}} are clockwise
}
\details{
Similarly to \code{\link[spatstat]{owin}}, the polygon should not be closed
}
\author{
Mikko Vihtakari, but the idea has been scavenged from https://stackoverflow.com/a/1165943/1082004
}
