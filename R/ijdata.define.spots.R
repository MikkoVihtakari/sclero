#' @title Define the spot sequences from an ImageJ zip file
#' @description An internal function which defines spot sequences
#' @param spots a character specifying the type of ROI that should be considered as sampling spot sequences. See the argument with the same name in \code{\link{read.ijdata}}.
#' @param types a data.frame of ROI types. Calculated internally.
#' @param rois a numeric value defining the scale of photograph in pixels / \code{unit}.
#' @param names a character argument. See the argument with the same name in \code{\link{read.ijdata}}.
#' @param spot.names a character vector specifying the manually added spot names. See the argument with the same name in \code{\link{read.ijdata}}.
#' @keywords internal
#' @return Returns a list of x and y coordinates for samling sequences.
#' @author Mikko Vihtakari
#' @export

ijdata.define.spots <- function(spots, types, scale, names, spot.names) {

  ## Add oval
  if(!is.null(spots)) {

    if(class(spots) == "character" & length(spots) == 1 & any(spots %in% types[,1])) hole.seqs <- which(types[,1] %in% spots)
    if(class(spots) == "character" & !any(spots %in% types[,1])) hole.seqs <- which(rownames(types) %in% spots)
    if(class(spots) == "integer" | class(spots) == "numeric") hole.seqs <- spots

    ## Find coordinates for spot sequences, x-axis

    tmp <-  lapply(hole.seqs, function(i) dat[[i]]$coords[,1]/scale) # Use scale argument to scale coordinates
    n <- max(unlist(lapply(tmp, length)))
    spots.x <- as.data.frame(do.call("cbind", lapply(tmp, function(x) {length(x) <- n; x})))
    spot.nms <- ifelse(is.null(spot.names) | ncol(spots.x) != length(spot.names) & names == "manual", "generate.invalid", names)

    ## X names
    colnames(spots.x) <- ijdata.define.names(names = spot.nms, types = types, rois = hole.seqs, data = spots.x, type.names = spot.names)

    ## Find coordinates for spot sequences, y-axis

    tmp <-  lapply(hole.seqs, function(i) dat[[i]]$coords[,2]/scale) # Use scale argument to scale coordinates
    n <- max(unlist(lapply(tmp, length)))
    spots.y <- as.data.frame(do.call("cbind", lapply(tmp, function(x) {length(x) <- n; x})))

    ## Y names
    colnames(spots.y) <- ijdata.define.names(names = spot.nms, types = types, rois = hole.seqs, data = spots.y, type.names = spot.names)

    ## Return

    list(spots.x = spots.x, spots.y = spots.y)

  } else {
    list(spots.x = NULL, spots.y = NULL)
  }
}

