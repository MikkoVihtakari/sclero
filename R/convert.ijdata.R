#' @title Convert IJDATA object to a list of spatstat.geom objects
#' @description Converts an \code{\link[=read.ijdata]{IJDATA}} to a list of \link[spatstat.geom]{spatstat.geom} patterns ready for \link[=plot.rawDist]{plotting} or \link[=spot.dist]{sample spot alignment}.
#' @param X an \code{\link[=read.ijdata]{IJDATA}} object to be converted.
#' @param Accuracy numeric value defining the number to round to. See \code{\link{round_any}}.
#' @return Returns a list of class \code{rawDist}, which contains \link[spatstat.geom]{spatstat.geom} point patterns. The returned \code{rawDist} can be plotted using the generic plotting command.
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
#' @import spatstat.geom
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
  
  range.x <- c(plyr::round_any(do.call("min", tmp), Accuracy, floor), plyr::round_any(do.call("max", tmp), Accuracy, ceiling))
  
  tmp <- lapply(c("main.y", "spots.y", "gbs.y", "growth.y"), function(i) {
    if(is.null(X[[i]])) {
      NULL
    } else {
      range(X[[i]], na.rm = TRUE)
    }
  })
  
  range.y <- c(round_any(do.call("min", tmp), Accuracy, floor), round_any(do.call("max", tmp), Accuracy, ceiling))
  
  owin <- spatstat.geom::owin(xrange = range.x, yrange =  range.y, unitname = X$unit)
  
  ## 2. Main axis ####
  
  x <- X$main.x
  x <- x[!is.na(x)]
  y <- X$main.y
  y <- y[!is.na(y)]
  
  main <- spatstat.geom::psp(x0 = x[1], x1 = x[2], y0 = y[1], y1 = y[2], marks = "main", window = owin)
  
  start <- spatstat.geom::ppp(x[1], y[1], marks = "start", window = owin)
  end <- spatstat.geom::ppp(x[2], y[2], marks = "end", window = owin)
  
  main.sl <- tan(spatstat.geom::angles.psp(main))
  
  main.rot <- spatstat.geom::rotate.psp(main, angle = pi - atan(main.sl))  ### Rotated main axis
  main.shift.x <- -main.rot$ends$x0
  main.shift.y <- -main.rot$ends$y0
  
  ## 3. Growth lines ####
  
  x <- X$gbs.x
  y <- X$gbs.y
  
  dups <- data.frame(x = duplicated(t(round(x, 0))), y = duplicated(t(y)))
  
  if(any(rowSums(dups) > 1)) {
    message("Growth bands ", paste(rownames(dups)[rowSums(dups) > 1], collapse = ", "), " are duplicated. Removed.")
    x <- x[!names(x) %in% rownames(dups)[rowSums(dups) > 1]]
    y <- y[!names(y) %in% rownames(dups)[rowSums(dups) > 1]]
  }
  
  if(ncol(x) != ncol(y)) stop("x and y coordinate data.frames for growth lines do not match")
  
  gbs.list <- lapply(1:ncol(x), function(j) {
    xx <- x[,j][!is.na(x[,j])]
    yy <- y[,j][!is.na(y[,j])]
    spatstat.geom::psp(x0 = xx[1:(length(xx)-1)], x1 = xx[2:length(xx)], y0 = yy[1:(length(yy)-1)], y1 = yy[2:length(yy)], marks = rep(colnames(x)[j], (length(xx)-1)), window = owin)
  })
  
  gbs <- do.call(spatstat.geom::superimpose, gbs.list)
  
  names(gbs.list) <- colnames(x)
  
  start.gbs <- marks(gbs[spatstat.geom::nearestsegment(start, gbs)])
  
  ## Part 3.2 Test whether the object is 'along' or 'cross' type (see the tutorial) ###
  
  test <- spatstat.geom::crossing.psp(main, spatstat.geom::superimpose(gbs))
  
  if(test$n == 0) {
    sec.type <- "along"
  } else {
    if(test$n == length(unique(gbs$marks))) {
      sec.type <- "cross"
    } else {
      stop(paste0("Number of growth line and main axis crossing points is neither 0 or ", length(unique(gbs$marks)), ". Cannot define the aligment type. See ?spot.dist"))
    }
  }
  
  ## 3.3. Closest distances
  
  # k <- gbs.list[[1]]
  tmp <- lapply(gbs.list, function(k) {
    
    if(sec.type == "along") {
      data.frame(gbs.name = unique(spatstat.geom::marks(k)), min = min(spatstat.geom::crossdist(gbs.list[[start.gbs]], k)), max = max(spatstat.geom::crossdist(gbs.list[[start.gbs]], k)), stringsAsFactors = FALSE)
    } else {
      crs.point <- spatstat.geom::crossing.psp(main, k)
      crs.point.start <- spatstat.geom::crossing.psp(main, gbs.list[[start.gbs]])
      data.frame(gbs.name = unique(spatstat.geom::marks(k)), min = spatstat.geom::crossdist(crs.point.start, crs.point))
    }
  })
  
  gbs.order <- do.call(rbind, tmp)
  gbs.list <- gbs.list[order(gbs.order$min)]
  
  gbs <- do.call(spatstat.geom::superimpose, gbs.list)
  
  marks(gbs) <- factor(marks(gbs)[[1]], unique(marks(gbs)[[1]]))
  
  gbs.rot <- spatstat.geom::rotate.psp(gbs, angle = pi - atan(main.sl)) ### Rotated growth lines
  gbs.rot <- spatstat.geom::shift.psp(gbs.rot, vec = c(main.shift.x, main.shift.y)) ### Shift growth lines
  
  ### Reflect around axes 
  
  refl.x <- sum(gbs.rot$ends$x1 < 0) > sum(gbs.rot$ends$x1 > 0) # Reflect around y-axis if more negative than positive x values 
  refl.y <- sum(gbs.rot$ends$y1 < 0) > sum(gbs.rot$ends$y1 > 0) # Reflect around x-axis if more negative than positive y values 
  
  if(refl.x) { 
    gbs.rot$ends[c("x0", "x1")] <- -1*gbs.rot$ends[c("x0", "x1")]
  }
  
  if(refl.y) { 
    gbs.rot$ends[c("y0", "y1")] <- -1*gbs.rot$ends[c("y0", "y1")]
  }
  
  ## 4. Sample spots ###
  
  if(is.na(X$parameters$spot.rois)) {
    
    spots <- NULL
    spots.rot <- NULL
    
  } else {
    spots <- lapply(1:ncol(X$spots.x), function(i) {
      x <- na.omit(X$spots.x[,i])
      y <- na.omit(X$spots.y[,i])
      spatstat.geom::ppp(x, y, marks = as.factor(seq(1, length(x))), window = owin)
    })
    
    names(spots) <- colnames(X$spots.x)
    
    ### Rotated and shifted spots
    
    spots.rot <- lapply(spots, function(k) {
      rot <- rotate.ppp(k, angle = pi - atan(main.sl))
      rot <- shift.ppp(rot, vec = c(main.shift.x, main.shift.y))
      
      if(refl.x) { # Reflect around y-axis if more negative than positive x values 
        rot$x <- -1*rot$x
      }
      
      if(refl.y) { # Reflect around x-axis if more negative than positive y values 
        rot$y <- -1*rot$y
      }
      
      rot  
    })
    
  }
  
  ## 5. Growth axis ###
  
  if(is.null(X$growth.x)) {
    growth <- NULL
    growth.rot <- NULL
    
  } else {
    
    x <- X$growth.x[!is.na(X$growth.x)]
    y <- X$growth.y[!is.na(X$growth.y)]
    
    if(length(x) != length(y)) stop("x and y coordinates for the growth axis do not match")
    
    growth <- spatstat.geom::psp(x0 = x[1:(length(x)-1)], x1 = x[2:length(x)], y0 = y[1:(length(y)-1)], y1 = y[2:length(y)], marks = rep("growth", (length(x)-1)), window = owin)
    
    growth.rot <- rotate.psp(growth, angle = pi - atan(main.sl)) ### Rotated growth axis
    growth.rot <- shift.psp(growth.rot, vec = c(main.shift.x, main.shift.y)) ### Shift growth axis
    
    if(refl.x) { 
      growth.rot$ends[c("x0", "x1")] <- -1*growth.rot$ends[c("x0", "x1")]
    }
    
    if(refl.y) { 
      growth.rot$ends[c("y0", "y1")] <- -1*growth.rot$ends[c("y0", "y1")]
    }
    
  }
  
  ## 6. Clean up the rotated coordinates, which will be used as main coordinates from now on ###
  
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
  
  main.rot <- spatstat.geom::shift.psp(main.rot, vec = c(main.shift.x, main.shift.y)) ## Shift to the main axis begin from start
  
  if(refl.x) { 
    main.rot$ends[c("x0", "x1")] <- -1*main.rot$ends[c("x0", "x1")]
  }
  
  if(refl.y) { 
    main.rot$ends[c("y0", "y1")] <- -1*main.rot$ends[c("y0", "y1")]
  }
  
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
  
  owin.rot <- spatstat.geom::owin(xrange = range.x, yrange = range.y, unitname = X$unit)
  
  start.rot <- spatstat.geom::ppp(main.rot$ends$x0, main.rot$ends$y0, marks = "start", window = owin.rot)
  end.rot <- spatstat.geom::ppp(main.rot$ends$x1, main.rot$ends$y1, marks = "end", window = owin.rot)
  
  Window(main.rot) <- owin.rot
  if(!is.null(spots.rot)) spots.rot <- lapply(spots.rot, function(k) {Window(k) <- owin.rot; k})
  Window(gbs.rot) <- owin.rot
  if(!is.null(growth.rot)) Window(growth.rot) <- owin.rot
  
  ## Parameters
  
  params <- X$parameters
  params <- c(params, main.sl = main.sl, main.shift.x = main.shift.x, main.shift.y = main.shift.y, refl.x = refl.x, refl.y = refl.y, sec.type = sec.type)
  
  ## Return ####
  
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
