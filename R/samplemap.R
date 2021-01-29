#' @title Common plotting function for plot.rawDist and plot.spotDist functions
#' @description Used as an internal function to make the code more comprehensible. See \code{\link{plot.rawDist}} and \code{\link{plot.spotDist}} for plotting \code{\link[=convert.ijdata]{rawDist}} and \code{\link[=spot.dist]{spotDist}} objects, respectively.
#' @param x \link[=convert.ijdata]{rawDist}} or \code{\link[=spot.dist]{spotDist} object.
#' @param ... Arguments to be passed to other methods, such as \link[=par]{graphical parameters}.
#' @param sname equals to \code{\link[=plot.rawDist]{sample.name}}
#' @param sptype equals to \code{\link[=plot.rawDist]{sample.type}}
#' @param size equals to \code{\link[=plot.rawDist]{spot.size}}
#' @param scol equals to \code{\link[=plot.rawDist]{spot.color}}
#' @param mtype equals to \code{\link[=plot.rawDist]{main.type}}
#' @param gaxis equals to \code{\link[=plot.rawDist]{growth.axis}}
#' @param colpalette equals to \code{\link[=plot.rawDist]{color.palette}}
#' @param hlight equals to \code{\link[=plot.rawDist]{highlight.gbs}}
#' @param hlcol equals to \code{\link[=plot.rawDist]{highlight.col}}
#' @author Mikko Vihtakari
#' @keywords internal
#' @seealso \code{\link{plot.rawDist}}, \code{\link{plot.spotDist}}
#' @import spatstat
#' @export

# x = X$scaled; sname = sample.name; sptype = spot.type; size = spot.size; scol = spot.color; mtype = main.type; gaxis = growth.axis; colpalette = color.palette; hlight = highlight.gbs; hlcol = highlight.col

samplemap <- function(x, ..., sname, sptype, size, scol, mtype, gaxis, colpalette, hlight, hlcol){
  
  ## Define the unit for axes ####
  unit <- x$unit
  
  ## Standard plot arguments (allows edits using the ...) ###
  
  arguments <- list(x = x$window$xrange, 
    y = x$window$yrange,
    type = "n",
    axes = F, 
    ..., 
    asp = 1, 
    xlab = paste0("x ", "(", unit, ")"), 
    ylab = paste0("y ", "(", unit, ")")
  )
  
  arguments <- arguments[!duplicated(names(arguments))]
  
  ## Plot header (main)
  if(sname == "keep") {
    samplename <- x$sample.name
  } else {
    samplename <- sname
  }
  
  ## Establish colors for hole sequences
  if(is.null(scol)) {
    hole.cols <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999")
  } else {
    hole.cols <- scol
  }
  
  ## Spots
  
  if(!is.null(x$mid.point.type)) {
    if(x$mid.point.type == "spots") {
      spots <- x$spots
    } else {
      if(x$mid.point.type == "centroids") {
        spots <- x$spot.area$centroids
      }
    }
  } else {
    if(is.null(x$spot.area)) {
      spots <- x$spots
    } else {
      spots <- x$spot.area$centroids
    }
  }
  
  ################
  ## Plotting ####
  
  ## Make layout if sptype == "value" | sptype == "idvalue"
  if(sptype == "value" | sptype == "idvalue") layout(matrix(c(1,2), nrow = 1), widths = c(9,1))
  
  par(mar = c(4, 4, 2.5, 1) + 0.1)
  
  do.call("plot", arguments) # To account for the modifiable arguments
  
  ## Title ####
  
  title(samplename, line = -0.2)
  
  ## Plot the axes
  
  axis(1)
  axis(2, las = 2)
  
  ## Plot the growth axis
  
  if(!is.null(x$growth) & !is.null(gaxis)) {
    plot(x$growth, add = T, col = gaxis)
  }
  
  ## Plot growth bands
  
  plot(x$gbs, col = "darkgrey", add = TRUE)
  
  ## Highlighted growth bands
  
  if(!is.null(hlight)) {
    plot(x$gbs[marks(x$gbs) %in% hlight], col = hlcol, add = TRUE)
  }
  
  ## Colors (if sptype == "value" or "idvalue")
  
  if(sptype %in% c("value", "idvalue")) {
    
    if(all(unlist(lapply(x$values, is.null)))) stop("'values' for the color mapping must be specified. See ?assign.value")
    color.spec = TRUE
    colmap <- colourmap(colpalette, range = range(unlist(lapply(x$values, function(k) range(k[,2], na.rm = TRUE))), na.rm = TRUE))
    COLOR <- lapply(seq_along(x$spots), function(i){
      ifelse(is.na(x$values[[i]][,2]), "white", colmap(x$values[[i]][,2]))
    })
    
  } else {
    
    color.spec = FALSE
  }
  
  ## Should sample number be printed?
  
  if(sptype %in% c("id", "idvalue")) id.spec = TRUE else id.spec = FALSE
  
  ## Plot sample spots
  
  if(!is.null(x$spots)) {
    
    if(size == "actual") {
      
      ## Plot actual sample spot size (if size == TRUE)
      if(is.null(x$spot.area$spot.dat)) stop("x does not contain spot size information. See ?assign.size")
      
      if(color.spec) {
        lapply(seq_along(x$spot.area$spot.dat), function(i) {
          lapply(seq_along(x$spot.area$spot.dat[[i]]$spot.owins), function(j) plot(x$spot.area$spot.dat[[i]]$spot.owins[[j]], add = TRUE, col = COLOR[[i]][[j]]))})
      } else {
        lapply(x$spot.area$spot.dat, function(k) {
          lapply(k$spot.owins, function(j) plot(j, add = TRUE))})
      }
      if(id.spec) lapply(seq_along(spots), function(i){
        text(spatstat::coords(spots[[i]])[,1], spatstat::coords(spots[[i]])[,2], spatstat::marks(spots[[i]]), col = hole.cols[i], cex = 0.6, font = 2)
      })
      
    } else {
      
      ## Other sample spot plotting operations
      if(color.spec) {
        lapply(seq_along(spots), function(i) {
          points(x = spatstat::coords(spots[[i]])[,1], y = spatstat::coords(spots[[i]])[,2], bg = COLOR[[i]], cex = size, pch = 21)
        })
      } else {
        lapply(seq_along(spots), function(i) {
          points(x = spatstat::coords(spots[[i]])[,1], y = spatstat::coords(spots[[i]])[,2], cex = size, pch = 21)
        })
      }
      if(id.spec) {
        lapply(seq_along(spots), function(i) {
          text(spatstat::coords(spots[[i]])[,1], spatstat::coords(spots[[i]])[,2], spatstat::marks(spots[[i]]), col = hole.cols[i], cex = 0.6, font = 2)
        })
      }
    }
  }
  
  ## Plot the main axis
  
  if(mtype %in% c("all", "axis")) plot(x$main, add = T, col = "darkslateblue")
  
  ## Plot the end and start positions for the main axis
  
  if(mtype %in% c("all", "ends")) {
    text(x$start.main$x, x$start.main$y, x$start.main$marks, col = "dark green", font = 2)
    text(coords(x$end.main)[,1], coords(x$end.main)[,2], x$end.main$marks, col = "dark green", font = 2)
  }
  
  if(color.spec) output <- list(spots = spots, hole.cols = hole.cols, colmap = colmap) else output <- list(spots = spots, hole.cols = hole.cols)
  return(output)
}
