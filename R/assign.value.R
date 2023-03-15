#' @title Assign values to 'rawDist' objects for plotting.
#'
#' @description Assigns values to \code{\link[=convert.ijdata]{rawDist}} objects for \link[=plot.rawDist]{spatial density plotting}.
#' 
#' @param rawDist \code{\link[=convert.ijdata]{rawDist}} object to which the values should be assigned.
#' @param value a \code{\link{data.frame}} or list of data.frames of the same length than spot spequences (\code{\link[=read.ijdata]{spots}}). Each \code{data.frame} has to have identical length to spots. First column specifies the order. Second column containing the values (see Details).
#' @details This function can be used to plot values as color-densities on \link[=plot.rawDist]{sample maps}. The function is useful e.g. for examining the spatial distribution of geochemical data, such as element ratios or isotope ratios, along sample materials. If the \code{\link[=convert.ijdata]{rawDist}} object contains only one sample spot sequence, the \code{value} parameter should be expressed as a data.frame with two columns. If the \code{rawDist} object consists of several sample spot sequences, the \code{value} parameter should be a list of data.frames with length equivalent to number of spot sequences. The first column in all \code{value} data.frames represents spot number and should be equivalent to \code{$spots} marks in the \code{rawDist} object. The second column represents the values to be assigned. Column names are ignored.
#' @param value.name a character or function defining the name for the assigned value. Can be \code{\link[base]{expression}}
#' @param coord.type A character or numeric argument defining the coordinate type to for assignment. Used only to acquire spot sequence names. Alternatives:
#' \itemize{
#'   \item \strong{"scaled"} or \strong{1}. Only scaled coordinates.
#'   \item \strong{"original"} or \strong{2}. Only original ImageJ ROI coordinates.
#'   }
#' @return Returns a list of class \code{\link[=convert.ijdata]{rawDist}} containing spot value information.
#' @author Mikko Vihtakari
#' @seealso \code{\link{read.ijdata}} \code{\link{spot.dist}} \code{\link{assign.size}} \code{\link{plot.rawDist}}
#' @examples data(barium)
#' data(shellsizes)
#' 
#' ## rawDist
#' shellvalues <- assign.value(shellsizes, barium, value.name = "Ba/Ca")
#' plot(shellvalues, spot.size = "actual", spot.type = "value", main.type = "none")
#' 
#' ## spotDist
#' shellvalues.aligned <- spot.dist(shellvalues)
#' plot(shellvalues.aligned, spot.size = "actual", spot.type = "idvalue", 
#' spot.color = "darkgrey", highlight.gbs = c("WG_start", "WG_end"))
#' @export
#' 

# value.name = NULL; coord.type = "scaled"
assign.value <- function(rawDist, value, value.name = NULL, coord.type = "scaled"){
  
  if(all(c("scaled", "original") %in% names(rawDist))) {
    if(coord.type == "scaled" | coord.type == 1) {
      X <- rawDist$scaled
    } else if(coord.type == "original" | coord.type == 2) {
      X <- rawDist$original
    } else {
      stop("coord.type has to be either 'scaled' or 'original'")  
    }
  } else {
    X <- rawDist
  }
  
  if(class(value) == "data.frame" & length(X$spots) == 1) {
    value <- list(value)
    names(value) <- names(X$spots)
  }  
  
  # i = 1
  TMP <- lapply(seq_along(X$spots), function(i){
    tmp.s <- X$spots[[i]]
    tmp.v <- value[[names(X$spots[i])]]
    if(is.null(tmp.v)) {
      tmp.ret <- data.frame(spot.marks = marks(tmp.s), value = rep(NA, npoints(tmp.s)))}
    else {
      if(class(tmp.v) == "data.frame") {
        if(ncol(tmp.v) != 2) 
          stop("Each data.frame has to have 2 columns. First column defines the spot number, second the value ot be added") else {
            tmp.ret <- merge(data.frame(spot.marks = marks(tmp.s)), tmp.v, 
                             by = 1, sort = FALSE, all = TRUE)
          }
      }
    }
    return(tmp.ret)
  })
  
  names(TMP) <- names(X$spots)
  out <- rawDist
  out$values <- TMP
  out$value.name <- value.name
  
  class(out) <- "rawDist"
  
  return(out)
}



