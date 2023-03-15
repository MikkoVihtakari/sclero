#' @title Check whether points for an \code{\link[spatstat]{owin}} are clockwise
#' @param x a data frame with x coordinates in the first column and y coordinates in the second. 
#' @details Similarly to \code{\link[spatstat]{owin}}, the polygon should not be closed
#' @return A logical telling whether the polygon is arranged clockwise.
#' @author Mikko Vihtakari, but the idea has been scavenged from https://stackoverflow.com/a/1165943/1082004
#' @export

is.clockwise <- function(x) {
  # x.ppp <- spatstat::as.ppp(x)
  # x.coords <- c(x.ppp$x
  # y.coords <- x.ppp$y
  
  x.coords <- c(x[[1]], x[[1]][1])
  y.coords <- c(x[[2]], x[[2]][1])
  
  double.area <- sum(sapply(2:length(x.coords), function(i) {
    (x.coords[i] - x.coords[i-1])*(y.coords[i] + y.coords[i-1])
  }))
  
  double.area > 0
} 
