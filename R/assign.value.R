#' @title Assign values to 'rawDist' objects for plotting.
#'
#' @description Assigns values to \code{\link[=convert.ijdata]{rawDist}} objects for \link[=plot.rawDist]{spatial density plotting}.
#' 
#' @param rawDist \code{\link[=convert.ijdata]{rawDist}} object to which the values should be assigned.
#' @param value A list of \code{\link{data.frame}}s of the same length than spot spequences (\code{\link[=read.ijdata]{spots}}). Each \code{data.frame} has to have identical length to spots. First column specifies the order. Second column containing the values.
#' @return Returns a list of class 'rawDist' containing spot value information.
#' @author Mikko Vihtakari
#' @export
#' 

assign.value <- function(rawDist, value){

TMP <- lapply(seq_along(rawDist$spots), function(k){
  tmp.s <- rawDist$spots[[k]]
  tmp.v <- value[[names(rawDist$spots[k])]]
  if(is.null(tmp.v)) {
    tmp.ret <- data.frame(spot.marks = marks(tmp.s), value = rep(NA, npoints(tmp.s)))} else {
    if(class(tmp.v) == "data.frame") {
      if(ncol(tmp.v) != 2) 
        stop("Each data.frame has to have 2 columns. First column defines the spot number, second the value") else {
          tmp.ret <- merge(data.frame(spot.marks = marks(tmp.s)), tmp.v, by = 1, sort = FALSE, all = TRUE)
        }
      }
    }
  return(tmp.ret)})
  
  names(TMP) <- names(rawDist$spots)
  X <- rawDist
  X$values <- TMP
  return(X)
}



