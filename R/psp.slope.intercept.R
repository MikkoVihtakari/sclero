## Slope and intercept functions needed in rotation of objects

#' @title Calculate slopes for a spatstat psp object
#' @description Calculates a slope for each line segment in a psp object
#' @param x a \code{\link[spatsta]{psp}} object
#' @param orthogonal logical indicating whether a slope for orthogonal line to \code{x} should be returned instead of slope for \code{x}
#' @return Returns a numeric vector containing slopes of line segments in the same order than \code{x}
#' @details The function assumes that each line segment has an infinite domain and follows the line equation $y = slope * x + intercept$.
#' @keywords internal
#' @author Mikko Vihtakari
#' @export

slope.psp <- function(x, orthogonal = FALSE) {
  out <- tan(angles.psp(x))
  
  if(orthogonal) {
    -1/out
  } else {
    out
  }
}

#' @title Calculate intercept for a spatstat psp object
#' @description Calculates an intercept for each line segment in a psp object
#' @param x a \code{\link[spatsta]{psp}} object
#' @return Returns a numeric vector containing y-intercepts of line segments in the same order than \code{x}
#' @details The function assumes that each line segment has an infinite domain and follows the line equation $y = slope * x + intercept$.
#' @keywords internal
#' @author Mikko Vihtakari
#' @export

intercept.psp <- function(x) {
  x$ends$y0 - tan(angles.psp(x)) * x$ends$x0
}

#' @title Find points that lie along line based on x- or y coordinates
#' @description Given a line segment pattern, find a series of points that lie along it based on x or y coordinates
#' @param X A line segment pattern (object of class "psp").
#' @param x,y Known x or y coordinates given as a numeric vector or a ppp object. One of the parameters has to be \code{NULL}. 
#' @return A point pattern (object of class "ppp") in the same window as X.
#' @author Mikko Vihtakari
#' @export

pointsAlongLines <- function(X = NULL, x = NULL, y = NULL) {
  
  ## Tests
  if(is.null(X)) stop("X has to be provided. Please add a psp object")
  if(!"psp" %in% class(X)) stop("X has to be a psp object.")
  if(is.null(x) & is.null(y)) stop("Either x or y coordinates have to be provided.")
  if(!is.null(x) & !is.null(y)) stop("Both x and y cannot be provided.")
  
  ## If y is provided
  if(!is.null(y)) {
    if(!(is.numeric(y) | "ppp" %in% class(y))) stop("y has to be a numeric vector or a ppp object")
    if("ppp" %in% class(y)) y <- y$y
      
    out <- lapply(y, function(k) {
      segs <- X[X$ends$y0 <= k & X$ends$y1 >= k]
      
      if(segs$n == 0) {
        spatstat::ppp(x = NULL, y = NULL, window = spatstat::Window(X))
      } else {
        spatstat::ppp(x = (k - intercept.psp(segs))/tan(spatstat::angles.psp(segs)), y = k, window = spatstat::Window(X))
      }
    })
  }
  
  ## If x is provided
  
  if(!is.null(x)) {
    if(!(is.numeric(x) | "ppp" %in% class(x))) stop("x has to be a numeric vector or a ppp object")
    if("ppp" %in% class(x)) x <- x$x
    
    if(!is.numeric(x))  stop("x has to be a numeric value or vector")
    
    out <- lapply(x, function(k) {
      segs <- X[X$ends$x0 <= k & X$ends$x1 >= k]
      
      if(segs$n == 0) {
        ppp(x = NULL, y = NULL, window = Window(X))
      } else {
        ppp(x = k, y = tan(spatstat::angles.psp(segs))*k + intercept.psp(segs), window = Window(X))
      }
    })
  }
  
  return(do.call(spatstat::superimpose, out))
}
## ####
