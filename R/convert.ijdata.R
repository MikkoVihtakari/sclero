#' @title Convert IJDATA object to a list of spatstat objects
#' @description Converts an \code{\link[=read.ijdata]{IJDATA}} to a list of \link[spatstat]{spatstat} patterns ready for \link[=plot.rawDist]{plotting} or \link[=spot.dist]{sample spot alignment}.
#' @param X an \code{\link[=read.ijdata]{IJDATA}} object to be converted.
#' @param Accuracy numeric value defining the number to round to. See \code{\link{round_any}}.
#' @return Returns a list of class \code{rawDist}, which contains \link[spatstat]{spatstat} point patterns. The returned \code{rawDist} can be plotted using the generic plotting command.
#' @author Mikko Vihtakari
#' @seealso \code{\link{read.ijdata}} for generating \code{IJDATA} objects.
#'
#' \code{\link{plot.rawDist}} for plotting.
#'
#' \code{\link{spot.dist}} for aligning sample spots.
#'
#' @examples data(shellspots)
#' shell_map <- convert.ijdata(shellspots)
#' plot(shell_map)
#' @import spatstat
#' @export

# X = dat; Accuracy = 0.1
convert.ijdata <- function(X, Accuracy = 0.1) {
  
  ## 1. Define observation window ####
  
  tmp <- lapply(c("main.x", "spots.x", "gbs.x", "growth.x"), function(i) {
    if(is.null(X[[i]])) {
      NULL
    } else {
      range(X[[i]], na.rm = TRUE)
    }
  })
  
  range.x <- c(round_any(do.call("min", tmp), Accuracy, floor), round_any(do.call("max", tmp), Accuracy, ceiling))
  
  tmp <- lapply(c("main.y", "spots.y", "gbs.y", "growth.y"), function(i) {
    if(is.null(X[[i]])) {
      NULL
    } else {
      range(X[[i]], na.rm = TRUE)
    }
  })
  
  range.y <- c(round_any(do.call("min", tmp), Accuracy, floor), round_any(do.call("max", tmp), Accuracy, ceiling))
  
  owin <- spatstat::owin(xrange = range.x, yrange =  range.y, unitname = X$unit)
  
  ## 2. Main axis ###
  
  x <- X$main.x
  x <- x[!is.na(x)]
  y <- X$main.y
  y <- y[!is.na(y)]
  
  main <- spatstat::psp(x0 = x[1], x1 = x[2], y0 = y[1], y1 = y[2], marks = "main", window = owin)
  
  start <- spatstat::ppp(x[1], y[1], marks = "start", window = owin)
  end <- spatstat::ppp(x[2], y[2], marks = "end", window = owin)
  
  main.sl <- slope.psp(main)
  
  main.rot <- rotate.psp(main, angle = pi - atan(main.sl))  ### Rotated main axis
  main.shift.x <- -main.rot$ends$x0
  main.shift.y <- -main.rot$ends$y0
  
  ## 3. Spots ###
  
  if(is.na(X$parameters$spot.rois)) {
    
    spots <- NULL
    spots.rot <- NULL
    
  } else {
    spots <- lapply(1:ncol(X$spots.x), function(i) {
      x <- na.omit(X$spots.x[,i])
      y <- na.omit(X$spots.y[,i])
      spatstat::ppp(x, y, marks = as.factor(seq(1, length(x))), window = owin)
    })
    
    names(spots) <- colnames(X$spots.x)
    
    ### Rotated and shifted spots
    
    spots.rot <- lapply(spots, function(k) {
      rot <- rotate.ppp(k, angle = pi - atan(main.sl))
      shift.ppp(rot, vec = c(main.shift.x, main.shift.y))
    })
    
  }
  
  ## 4. Growth lines ###
  
  x <- X$gbs.x
  y <- X$gbs.y
  
  if(ncol(x) != ncol(y)) stop("x and y coordinate data.frames for growth lines do not match")
  
  gbs.list <- lapply(1:ncol(x), function(j) {
    xx <- x[,j][!is.na(x[,j])]
    yy <- y[,j][!is.na(y[,j])]
    spatstat::psp(x0 = xx[1:(length(xx)-1)], x1 = xx[2:length(xx)], y0 = yy[1:(length(yy)-1)], y1 = yy[2:length(yy)], marks = rep(colnames(x)[j], (length(xx)-1)), window = owin)
  })
  
  gbs <- do.call(spatstat::superimpose, gbs.list)
  
  names(gbs.list) <- colnames(X$gbs.x)
  
  start.gbs <- marks(gbs[spatstat::nearestsegment(start, gbs)])
  
  # k <- gbs.list[[1]]
  tmp <- lapply(gbs.list, function(k) {
    data.frame(gbs.name = unique(spatstat::marks(k)), min = min(spatstat::crossdist(gbs.list[[start.gbs]], k)), max = max(spatstat::crossdist(gbs.list[[start.gbs]], k)), stringsAsFactors = FALSE)
  })
  
  gbs.order <- do.call(rbind, tmp)
  gbs.list <- gbs.list[order(gbs.order$min)]
  
  gbs <- do.call(spatstat::superimpose, gbs.list)
  
  marks(gbs) <- factor(marks(gbs)[[1]], unique(marks(gbs)[[1]]))
  
  gbs.rot <- rotate.psp(gbs, angle = pi - atan(main.sl)) ### Rotated growth lines
  gbs.rot <- shift.psp(gbs.rot, vec = c(main.shift.x, main.shift.y)) ### Shift growth lines
  
  ## 5. Growth axis ###
  
  if(is.null(X$growth.x)) {
    growth <- NULL
    growth.rot <- NULL
    
  } else {
    
    x <- X$growth.x[!is.na(X$growth.x)]
    y <- X$growth.y[!is.na(X$growth.y)]
    
    if(length(x) != length(y)) stop("x and y coordinates for the growth axis do not match")
    
    growth <- spatstat::psp(x0 = x[1:(length(x)-1)], x1 = x[2:length(x)], y0 = y[1:(length(y)-1)], y1 = y[2:length(y)], marks = rep("growth", (length(x)-1)), window = owin)
    
    growth.rot <- rotate.psp(growth, angle = pi - atan(main.sl)) ### Rotated growth axis
    growth.rot <- shift.psp(growth.rot, vec = c(main.shift.x, main.shift.y)) ### Shift growth axis
  }
  
  ## 6. Clean up the rotated coordinates, which will be used as main coordinates from now on ####
  
  range.fun <- function(k) {
    if(is.null(k)) {
      NULL
    } else if("psp" %in% class(k)) {
      xrange <- range(k$ends[grepl("x", names(k$ends))], na.rm = TRUE)
      yrange <- range(k$ends[grepl("y", names(k$ends))], na.rm = TRUE)
      
      list(x = xrange, y = yrange)
    } else if("ppp" %in% class(k)) {
      xrange <- range(k$x, na.rm = TRUE)
      yrange <- range(k$y, na.rm = TRUE)
      
      list(x = xrange, y = yrange)
    } else {
      stop("Something went wrong in step 6")
    }
  }
  
  main.rot <- shift.psp(main.rot, vec = c(main.shift.x, main.shift.y)) ## Shift to the main axis begin from start
  
  tmp <- lapply(list(main.rot, spots.rot, gbs.rot, growth.rot), function(k) {
    
    if(all(class(k) %in% "list")) {
      out <- lapply(k, function(g) range.fun(g))
      list(x = unname(unlist(lapply(out, function(g) g$x))), y = unname(unlist(lapply(out, function(g) g$y))))
    } else {
      range.fun(k)
    }
  })
  
  range.x <- c(round_any(do.call("min", lapply(tmp, function(k) k$x)), Accuracy, floor), round_any(do.call("max", lapply(tmp, function(k) k$x)), Accuracy, ceiling))
  range.y <- c(round_any(do.call("min", lapply(tmp, function(k) k$y)), Accuracy, floor), round_any(do.call("max", lapply(tmp, function(k) k$y)), Accuracy, ceiling))
  
  owin.rot <- spatstat::owin(xrange = range.x, yrange = range.y, unitname = X$unit)
  
  start.rot <- spatstat::ppp(main.rot$ends$x0, main.rot$ends$y0, marks = "start", window = owin.rot)
  end.rot <- spatstat::ppp(main.rot$ends$x1, main.rot$ends$y1, marks = "end", window = owin.rot)
  
  Window(main.rot) <- owin.rot
  if(!is.null(spots.rot)) spots.rot <- lapply(spots.rot, function(k) {Window(k) <- owin.rot; k})
  Window(gbs.rot) <- owin.rot
  if(!is.null(growth.rot)) Window(growth.rot) <- owin.rot
  
  ## Parameters
  
  params <- X$parameters
  params <- c(params, main.sl = main.sl, main.shift.x = main.shift.x, main.shift.y = main.shift.y)
  
  ## Return ----
  
  df <- 
    list(
      scaled = list(main = main.rot, gbs = gbs.rot, spots = spots.rot, growth = growth.rot,
                    window = owin.rot, start.main = start.rot, end.main = end.rot, sample.name = X$sample.name,
                    scaling.factor = X$scaling.factor, unit = X$unit, parameters = params),
      original = list(main = main, gbs = gbs, spots = spots, growth = growth, 
                      window = owin, start.main = start, end.main = end, sample.name = X$sample.name,
                      scaling.factor = X$scaling.factor, unit = X$unit, parameters = params)
    )
  
  class(df) <- "rawDist"
  
  return(df)
}
