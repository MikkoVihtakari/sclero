% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auto.detect.R
\name{auto.detect}
\alias{auto.detect}
\title{Automatically detect the most likely action from a set of ImageJ ROIs}
\usage{
auto.detect(roi.types)
}
\arguments{
\item{roi.types}{character vector containing the ROI types in an ImageJ zip file.}
}
\value{
Returns a list of the required parameter options
}
\description{
An internal function which tries to automatically detect the spots, gbs, and main parameters for the \code{\link{read.ijdata}} function.
}
\author{
Mikko Vihtakari
}
\keyword{internal}
