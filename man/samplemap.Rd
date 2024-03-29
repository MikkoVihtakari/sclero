% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/samplemap.R
\name{samplemap}
\alias{samplemap}
\title{Common plotting function for plot.rawDist and plot.spotDist functions}
\usage{
samplemap(
  x,
  ...,
  sname,
  sptype,
  size,
  scol,
  mtype,
  gaxis,
  colpalette,
  hlight,
  hlcol,
  values = NULL
)
}
\arguments{
\item{x}{\link[=convert.ijdata]{rawDist}} or \code{\link[=spot.dist]{spotDist} object.}

\item{...}{Arguments to be passed to other methods, such as \link[=par]{graphical parameters}.}

\item{sname}{equals to \code{\link[=plot.rawDist]{sample.name}}}

\item{sptype}{equals to \code{\link[=plot.rawDist]{sample.type}}}

\item{size}{equals to \code{\link[=plot.rawDist]{spot.size}}}

\item{scol}{equals to \code{\link[=plot.rawDist]{spot.color}}}

\item{mtype}{equals to \code{\link[=plot.rawDist]{main.type}}}

\item{gaxis}{equals to \code{\link[=plot.rawDist]{growth.axis}}}

\item{colpalette}{equals to \code{\link[=plot.rawDist]{color.palette}}}

\item{hlight}{equals to \code{\link[=plot.rawDist]{highlight.gbs}}}

\item{hlcol}{equals to \code{\link[=plot.rawDist]{highlight.col}}}

\item{values}{the location of \code{\link{assign.value}} data.}
}
\description{
Used as an internal function to make the code more comprehensible. See \code{\link{plot.rawDist}} and \code{\link{plot.spotDist}} for plotting \code{\link[=convert.ijdata]{rawDist}} and \code{\link[=spot.dist]{spotDist}} objects, respectively.
}
\seealso{
\code{\link{plot.rawDist}}, \code{\link{plot.spotDist}}
}
\author{
Mikko Vihtakari
}
\keyword{internal}
