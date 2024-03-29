% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dataset-docs.R
\docType{data}
\name{shellspots}
\alias{shellspots}
\title{Coordinates of sample spots and growth lines along a bivalve mollusk shell section}
\format{
\code{IJDATA} object. A list containing ImageJ ROI information about sample spots and growth lines along a bivalve mollusk shell section. The information is acquired by running "shellspots.zip" dataset through the \code{\link{read.ijdata}} function (see Examples). The list contains following data frames / vectors:
\itemize{
  \item spots.x. A data frame containing x-coordinates for sample spots with each sample spot sequence in its own column.
  \item spots.y. A data frame containing y-coordinates for sample spots.
  \item gbs.x. A data frame containing x-coordinates for growth lines with each growth line in its own column.
  \item gbs.y. A data frame containing y-coordinates for growth lines
  \item main.x. A data frame containing x-coordinates for the measurement axis 
  \item main.y. A data frame containing y-coordinates for the measurement axis 
  \item sample.name. A vector containing the sample name
  \item scaling.factor. A numeric value defining the scale of photograph in pixels / \code{unit}.
  \item unit. A charater vector defining the unit of measurements.
}
}
\usage{
data(shellspots)
}
\description{
Coordinates of sample spots and growth lines along a bivalve mollusk shell section
}
\keyword{datasets}
