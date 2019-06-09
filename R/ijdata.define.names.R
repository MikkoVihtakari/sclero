#' @title Define names for ROI data
#' @description An internal function which defines names for various ROI objects
#' @param names see the \code{\link[=read.ijdata]{names}} argument in \code{read.ijdata}
#' @param types a data.frame of ROI types. Calculated internally
#' @param rois integer vector defining the order of ROIs in \code{types} which will be affected by the function
#' @param data a data.frame containing the coordinate data. Calculated internally
#' @param type.names a character vector defining the manual names. See the \code{\link[=read.ijdata]{spot.names,gbs.names}} argument in \code{read.ijdata}.
#' @return Returns a character vector of names
#' @keywords internal
#' @author Mikko Vihtakari
#' @export

# rois = hole.seqs; data = spots.x; type.names = spot.names
ijdata.define.names <- function(names, types, rois, data, type.names) {

  switch(names,
         force.keep = rownames(types)[rois],
         keep = make.names(row.names(types))[rois],
         generate.invalid = {
           change <- make.names(row.names(types))[rois] != row.names(types)[rois]
           ifelse(change, paste0("s", 1:ncol(data)), row.names(types)[rois])
         },
         generate = paste0("s", 1:ncol(data)),
         manual = {
           if(is.null(type.names)) stop("type.names not provided")
           if(ncol(data) != length(type.names)) stop("number of data and type.names differ")
           type.names
         },
         stop("Unrecognized names argument")
  )
}
