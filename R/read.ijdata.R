#' @title Read ImageJ zip file containing several ROI files and extract coordinate information.
#'
#' @description A wrapper function, which reads an ImageJ zip file containing a collection of ROI files and outputs a list of data frames ready for \code{\link{convert.ijdata}} function.
#'
#' @param X a character string defining the name (including extension) or file path of an ImageJ zip file. Alternatively an \code{\link[=read.ijzip]{ijzip}} object.
#' @param spots,gbs optional. Sample spots and growth lines, respectively. One of the following options:
#' \itemize{
#'   \item \strong{\code{NULL}} (default). Type of the ROI object will be chosen automatically. Leads to selection of "point" and "polyline" objects as sample spots and growth lines.
#'   \item \strong{numeric or character vector} specifying the order or the names of elements of ROI objects that should be assigned as the corresponding type. Might be useful in debugging.
#'   \item \strong{a character string} specifying the type of ROI objects that should be considered as the corresponding type. This option is included for package development and backwards compatibility.
#' }
#' @param main optional. Main (measurement) axis. Same options than for \code{spots,gbs}. Only one measurement axis per ImageJ .zip file is allowed. Defaults to \code{"line"}.
#' @param growth optional. A character string specifying the name of the "polyline" ROI that was used as a growth axis. Only if one was used, otherwise use \code{NULL}.
#' @param names optional. A character argument specifying how the names of \code{spots} and \code{gbs} should be generated. These names will be used in further functions (\code{\link{convert.ijdata}}, \code{\link{spot.dist}}). In general, it is adviced to use simple ROI names without special characters (for example \code{-} is not allowed in a ROI name; see 'Details'). Possible \code{names} options are:
#' \itemize{
#' \item \code{"generate.invalid"} (default). Uses the ROI object names, except when they are not valid \code{\link{data.frame}} column names. In the latter case sequential names will be generated.
#' \item \code{"generate"}. Generates sequential names for all elements.
#' \item \code{"keep"}. Uses the ROI object names, except when they are not valid \code{\link{data.frame}} column names. In the latter case \code{\link{make.names}} function will be used to generate \code{data.frame} combatible column names.
#' \item \code{"force.keep"}. Uses the ROI object names as they are in the .zip file. Using this option might cause problems in consequent functions and is not recommended.
#' \item \code{"manual"}. Names for \code{spots} and \code{gbs} are searched from \code{spot.names} and \code{gbs.names} arguments, respectively. If these arguments are not \code{NULL} and have the same length than \code{spots} or \code{gbs}, the given names are used for the elements. If the condition is not true, the names are generated following \code{"generate.invalid"}.
#' }
#' @param spot.names,gbs.names optional. A character vector of equal length to \code{spots}/\code{gbs} defining the names of sample spot sequences/growth bands.
#' @param main.name optional. A character string defining the name of the measurement axis (\code{main}). If \code{main.name = "keep"}, the ROI object name will be used (not recommended, see "Details"). Otherwise the name will be taken from the argument. Defaults to \code{"main"}.
#' @param sample.name optional. A character string defining the name of the sample. File name without the extension or alternatively \code{ijzip} object name is used when \code{NULL} (default).
#' @param scale optional. A numeric value defining the scale of photograph in pixels / \code{unit}. Defaults to 1.
#' @param unit optional. A charater string defining the unit of measurements. See \code{scale}.
#'
#' @details In order to minimize the amount of text to be typed by a user, ROI objects of type "point" (this includes the "Multi-point Tool" points) are considered as sample spot sequences (\code{spots}) by default. Further, all "polyline" objects are assumed as growth bands (\code{gbs}) and "line" objects as the measurement axis (\code{main}) resulting to that only one "line" object is allowed per .zip file using the default settings. Alternatively, the user can specify the \code{spots}, \code{gbs}, and \code{main} objects manually using the order of the ImageJ .zip file with the exception that \bold{only one measurement axis is allowed} per \code{\link[=convert.ijdata]{rawDist}} or \code{\link[=spot.dist]{spotDist}} object.
#'
#'Punctuation characters other than \code{_} or \code{.} should not be used as names of \code{spots} or \code{gbs}, because they tend to confuse the internal \code{\link[base]{grep}} functions in \code{\link{spot.dist}} function. Hence it is adviced to use one of the options renaming invalid names of \code{spots} and \code{gbs} (\code{"generate.invalid"}, \code{"generate"}, \code{"keep"}).
#'
#' Custom growth axis, if different from \code{main}, can be defined by using one named "polyline" in the ROI file and typing the name to \code{growth} argument. In such case, \code{main} is used to define the direction of growth, while the \code{growth} axis is used in growth measurements and sample spot alignment (distances are projected first to the \code{growth} axis and then transfered to the \code{main} axis).
#'
#' @return Returns an "IJDATA" object, which is a list of data frames containing the x and y coordinates for sampling spot sequences (\code{spots.x} and \code{spots.y}), growth bands (\code{gbs.x} and \code{gbs.y}), measurement axis (\code{main.x} and \code{main.y}) and growth axis (\code{growth.x} and \code{growth.y}) together with sample name, scaling factor and unit of measurement.
#' @author Mikko Vihtakari
#' @seealso \code{\link{order.ijdata}} for ordering and subsetting \code{read.ijdata} output.
#'
#' \code{\link{convert.ijdata}} for converting the coordinate information to \link[spatstat]{spatstat} point patterns.
#'
#' \code{\link{spot.dist}} for aligning sample spot sequences.
#'
#' \code{\link[RImageJROI]{read.ijroi}} and \code{\link[RImageJROI]{read.ijzip}} for reading ImageJ ROI and .zip files.
#'
#' @examples
#' # Locate the example zip file
#' path <- file.path(system.file("extdata", package = "sclero"), "shellspots.zip")
#'
#' # You can replace 'path' by 'Your_file_name.zip'
#' dat <- read.ijdata(path)
#' summary(dat)
#'
#' ## Works also for IJZIP objects
#' dat2 <- read.ijzip(path)
#' dat2 <- read.ijdata(dat2)
#' dat[!(dat %in% dat2)] # Only the sample name differs
#' @import RImageJROI
#' @export

## Debugging parameters
# X = file; spots = NULL; gbs = NULL; main = NULL; growth = NULL; names = "generate.invalid"; spot.names = NULL; gbs.names = NULL; main.name = "main"; sample.name = NULL; scale = 174; unit = "mm"
# growth = "growth"

read.ijdata <- function(X, spots = NULL, gbs = NULL, main = NULL, growth = NULL, names = "generate.invalid", spot.names = NULL, gbs.names = NULL, main.name = "main", sample.name = NULL, scale = 1, unit = NULL){

## Read zip file or load IJZIP object ----
if(class(X) == "character") {
  tmp <- grep(".zip", X, value = T)
  if(length(tmp) != 1) stop("Something wrong with the file name")
  if(length(tmp) == 1 & tmp==X){
  dat <- RImageJROI::read.ijzip(X, names = TRUE)
  filepath <- X
    }
  } else {
    if(class(X) == "ijzip") {
      dat <- X
      filepath <- NULL
      } else {
        stop("X is not path to a ImageJ .zip file nor a ijzip object")
  }}

## Define parameters ---

types <- as.data.frame(do.call("rbind", lapply(dat, function(x) x$strType)))
params <- auto.detect(roi.types = types[[1]])

### Redefine parameters

spots <- ifelse(is.null(spots), params$spot.rois, spots)
gbs <- ifelse(is.null(gbs), params$gbs.rois, gbs)
main <- ifelse(is.null(main), params$main.roi, main)

################
## 1. Spots ####

## Add oval
if(is.na(spots)) {
  spots.x = NULL; spots.y = NULL
} else if (spots == "point") {

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

 } else {
   spots.x = NULL; spots.y = NULL
}


#######################
## 2. Growth bands ####

if(class(gbs) == "character" & length(gbs) == 1 & any(gbs %in% types[,1])) {

  if(is.null(growth)) {

    gbss <- which(types[,1] %in% gbs)

    } else {

      if(length(growth) != 1) stop("the growth argument must have a length of 1")
      if(!growth %in% rownames(types)) stop(paste("name", growth, "not found from the ROI file. Check the growth argument"))
      if(unname(types[growth,]) != "polyline") stop("growth axis must be a polyline (i.e. segmented line ROI)")

      gbss <- which(types[,1] %in% gbs)
      growth.id <- which(gbss == which(rownames(types) == growth))

    }
  } else if(class(gbs) == "integer" | class(gbs) == "numeric") {

      gbss <- gbs

    } else {
      stop("did not find growth lines")
}

## Find coordinates, x-axis ####

tmp <-  lapply(gbss, function(i) dat[[i]]$coords[,1]/scale) # Use scale argument to scale coordinates
n <- max(unlist(lapply(tmp, length)))
gbs.x <- as.data.frame(do.call("cbind", lapply(tmp, function(x) {length(x) <- n; x})))
gbs.nms <- ifelse(is.null(gbs.names) | ncol(gbs.x) != length(gbs.names) & names == "manual", "generate.invalid", names)

## X names
colnames(gbs.x) <- ijdata.define.names(names = gbs.nms, types = types, rois = gbss, data = gbs.x, type.names = gbs.names)

## Find coordinates, y-axis

tmp <-  lapply(gbss, function(i) dat[[i]]$coords[,2]/scale) # Use scale argument to scale coordinates
n <- max(unlist(lapply(tmp, length)))
gbs.y <- as.data.frame(do.call("cbind", lapply(tmp, function(x) {length(x) <- n; x})))

## Y names
colnames(gbs.y) <- ijdata.define.names(names = gbs.nms, types = types, rois = gbss, data = gbs.y, type.names = gbs.names)

######################
## 3. Growth axis ####

if(is.null(growth)) {

  growth.x <- NULL
  growth.y <- NULL

  } else {

  growth.x <- gbs.x[growth.id]
  growth.y <- gbs.y[growth.id]

  gbs.x <- gbs.x[-growth.id]
  gbs.y <- gbs.y[-growth.id]
}

####################
## 4. Main axis ####

if(length(main) != 1) stop("Only one main axis allowed")
if(class(main) == "character" & any(main %in% types[,1])) main.l <- which(types[,1] %in% main)
if(class(main) == "character" & !any(main %in% types[,1])) main.l <- which(rownames(types) %in% main)
if(class(main) == "integer" | class(main) == "numeric") main.l <- main

## Find coordinates

main.x <- data.frame(main = dat[[main.l]]$coords[,1]/scale) # Use scale argument to scale coordinates
main.y <- data.frame(main = dat[[main.l]]$coords[,2]/scale)

## Names

if(length(main.name) != 1) {
  stop("The required length of main.name is 1. Use either 'keep' or a custom name.")} else {
    if(main.name != "keep") {
      colnames(main.x) <- main.name
      colnames(main.y) <- main.name
    }
  }


######################
## 4. Sample name ####

if (is.null(sample.name)) {
  if(class(X) == "ijzip") {
    name <- deparse(substitute(X))
  } else {
    if(grepl(".zip", X)) {
      if(grepl('/', X)) {
        tmp <- unlist(strsplit(X, '/'))
        name <- sub(".zip", "", tmp[length(tmp)])
      } else {
        name <- sub(".zip", "", X)}
    }
  }} else {
    if(length(sample.name) > 1) stop("Length of sample name > 1")
    name <- sample.name
}

############################
## 5. Compile to a list ####

dat <- list(spots.x = spots.x, spots.y = spots.y, gbs.x = gbs.x, gbs.y = gbs.y, main.x = main.x, main.y = main.y, growth.x = growth.x, growth.y = growth.y, sample.name = name, scaling.factor = scale, unit = unit, parameters = params, file.path = filepath)

class(dat) <- "IJDATA"

dat

## END ####
}

