#' @title Automatically detect the most likely action from a set of ImageJ ROIs
#' @description An internal function which tries to automatically detect the spots, gbs, and main parameters for the \code{\link{read.ijdata}} function.
#' @param roi.types character vector containing the ROI types in an ImageJ zip file.
#' @return Returns a list of the required parameter options
#' @keywords internal
#' @author Mikko Vihtakari
#' @export

# roi.types = types[[1]]
auto.detect <- function(roi.types) {

  ## Checks
  if(is.factor(roi.types)) roi.types <- as.character(roi.types)

  if(sum(roi.types == "line") > 1) stop("Only one *Straight* line object is allowed per an ImageJ zip file as this line defines the main axis")
  if(!any(roi.types %in% "line")) stop("One *Straight* line object has to be provided in an ImageJ zip file.
    This line will be used as the main axis to which all other measurements are scaled. See ?read.ijdata")

  ## Definitions

  unique.types <- sort(unique(roi.types))
  detect.types <- unique.types[unique.types %in% c("point", "polyline")]
  avg.error <- any(unique.types %in% "oval")

  ## Detect the required action, see the tutorial for action schematic



  if (all(c("point", "polyline") %in% detect.types)) {
    list(type = "spotDist", spot.rois = "point", gbs.rois = "polyline", main.roi = "line", avg.error = ifelse(avg.error, "oval", NA))

    } else if (detect.types == "point") {
      list(type = "growthDist", spot.rois = NA, gbs.rois = "point", main.roi = "line", avg.error = NA)

    } else if (detect.types == "polyline") {
      list(type = "growthDist", spot.rois = NA, gbs.rois = "polyline", main.roi = "line", avg.error = NA)

    } else {
    stop("The auto.detect function does not seem to work with your ImageJ zip file. The file contains unforseen ROIs. Please adjust the spots, gbs, main and avg.error parameters manually (see ?read.ijdata) or contact the maintainer.")
  }

}
