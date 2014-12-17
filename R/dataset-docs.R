#' Coordinates of sample spots and growth lines along a bivalve mollusk shell section
#'
#' @docType data
#' @keywords datasets
#' @name shellspots
#' @usage data(shellspots)
#' @format \code{IJDATA} object. A list containing ImageJ ROI information about sample spots and growth lines along a bivalve mollusk shell section. The information is acquired by running \link{ijzip} dataset through the \code{\link{read.ijdata}} function. The list contains following data frames / vectors:
#' \itemize{
#'   \item spots.x. A data frame containing x-coordinates for sample spots with each sample spot sequence in its own column.
#'   \item spots.y. A data frame containing y-coordinates for sample spots.
#'   \item gbs.x. A data frame containing x-coordinates for growth lines with each growth line in its own column.
#'   \item gbs.y. A data frame containing y-coordinates for growth lines
#'   \item main.x. A data frame containing x-coordinates for the measurement axis 
#'   \item main.y. A data frame containing y-coordinates for the measurement axis 
#'   \item sample.name. A vector containing the sample name
#'   \item scaling.factor. A numeric value defining the scale of photograph in pixels / \code{unit}.
#'   \item unit. A charater vector defining the unit of measurements.
#' }
NULL

#' \code{ijzip} object from \code{\link[RImageJROI]{read.ijzip}} function (RImageJROI package)
#'
#' @docType data
#' @keywords datasets
#' @name ijzip
#' @usage data(ijzip)
#' @format \code{ijzip} object obtained using \code{\link[RImageJROI]{read.ijzip}} function from RImageJROI package. A list containing ImageJ ROI information about sample spots and growth lines along a bivalve mollusk shell section.
NULL