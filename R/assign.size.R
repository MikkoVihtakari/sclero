#' @title Assign spot sizes to 'rawDist' objects for estimating spatial extent of sample averaging error.
#' 
#' @description Assigns spot sizes to \code{\link[=convert.ijdata]{rawDist}} objects for \link[=spot.dist]{estimating averaging error}.
#' 
#' @param rawDist \code{\link[=convert.ijdata]{rawDist}} object to which the values should be assigned.
#' @param file optional. ImageJ .zip file containing the spot size information. If \code{NULL} (default), the file name is assumed to be the same than from where \code{rawDist} \link[=read.ijdata]{data was read from}.
#' @param path optional. A character argument specifying the location of the \code{file}. If \code{NULL} (default), the \code{file} is assumed to be located in the \link[=getwd]{working directory}. See \code{\link[base]{dir}} for further information.
#' @param names optional. A character argument specifying how the names of \code{spots} should be generated. See \code{\link{read.ijdata}} for details. Defaults to "generate.invalid".
#' @param types optional. A character vector specifying the \code{strType} of ROI objects to be considered as sample spots (see \code{\link{plot.ijroi}} for possible pattern types). Defaults to \code{c("oval", "freehand", "rect")} meaning that oval and freehand selections, as well as rectangle tool selections will be used to calculate the spatial extent of sample spots.
#' @param scaled logical indicating whether to use scaled (\code{TRUE}) or original (\code{FALSE}) coordinates. Does not work well for scaled coordinates yet.
#' @return Returns a list of class 'rawDist' with a list of \code{\link{ppp}} objects containing locations of sample spot centroids and a list of \code{\link{hyperframe}}s containing spot size information.
#' @details If the .zip file containing spot size information is the same than from which the \code{rawDist} object was derived from and located in your working directory, assignment of spot sizes is simply specified by \code{assign.size(rawDist)}. Otherwise, use the \code{path} argument to specify the folder where the \code{file} is located.
#' @author Mikko Vihtakari
#' @seealso \code{\link{read.ijdata}} \code{\link{spot.dist}} \code{\link{assign.value}} \code{\link{plot.spotDist}}  
#' @examples data(shellspots)
#' shell <- convert.ijdata(shellspots)
#' path <- file.path(system.file("extdata", package = "sclero"))
#' sizes <- assign.size(shell, path = path)
#' sizes$spot.area
#' @import spatstat spatstat.geom RImageJROI
#' @export
#' 

# file = NULL; path = NULL; names = "generate.invalid"; types = c("oval", "freehand", "rect"); scaled = FALSE
assign.size <- function(rawDist, file = NULL, path = NULL, names = "generate.invalid", types = c("oval", "freehand", "rect"), scaled = FALSE) {
  
  ## ####
  
  if(scaled) {
    Z <- rawDist$scaled
  } else {
    Z <- rawDist$original
  }
  
  ## 1. Read in the .zip file and convert to spatstat
  
  if(is.null(file)) file <- paste0(Z$sample.name, ".zip")
  if(is.null(path)) path <- getwd()
  if(!file %in% dir(path)) {
    stop(paste0("Cannot find ", file, ". Check getwd() or define path."))
  }
  
  X <- RImageJROI::read.ijzip(file.path(path,file), names = TRUE)
  
  spot.owins <- RImageJROI::ij2spatstat(X, 
                                        convert.only = types, 
                                        scale = Z$scaling.factor, 
                                        unitname = Z$unit
  )
  
  if(length(spot.owins) == 0) {
    warning(paste0("Could not find ", 
                   paste(types, collapse = ", "), " objects from ", 
                   file.path(path,file), ". Returning ", 
                   ifelse(scaled, "scaled", "original"), " rawDist object")
    )
    class(Z) <- "rawDist"
    return(Z)
  }
  
  # k <- spot.owins[[1]]
  ## Scale coordinates
  
  if(scaled) {
    spot.owins <- lapply(spot.owins, function(k) {
      
      k <- rotate.owin(k, angle = pi - atan(Z$parameters$main.sl))
      k <- shift.owin(k, vec = c(Z$parameters$main.shift.x, Z$parameters$main.shift.y))
      
      if(Z$parameters$refl.x & Z$parameters$refl.y) {
        
        k$xrange <- sort(-1*k$xrange)
        k$bdry[[1]]$x <- -1*k$bdry[[1]]$x
        k$yrange <- sort(-1*k$yrange)
        k$bdry[[1]]$y <- -1*k$bdry[[1]]$y
        
      } else if(Z$parameters$refl.x) {
        
        k$xrange <- sort(-1*k$xrange)
        k$bdry[[1]]$x <- rev(-1*k$bdry[[1]]$x)
        
      } else if(Z$parameters$refl.y) {
        
        k$yrange <- sort(-1*k$yrange)
        k$bdry[[1]]$y <- rev(-1*k$bdry[[1]]$y)
        
      }
      
      k
      
      # ## Debug plotting for k, ignore
      # tmp <- data.frame(x = k$bdry[[1]]$x, y = k$bdry[[1]]$y, id = seq_along(k$bdry[[1]]$x))
      # ggplot2::ggplot(tmp, ggplot2::aes(x = x, y = y, label = id,color = as.numeric(id))) + 
      #   ggplot2::geom_path() + 
      #   ggplot2::geom_text(size = 2) + 
      #   ggplot2::scale_color_distiller(palette = "Spectral")
      # 
    })
  }
  
  ## 2. Relate the spot areas with sample spots  ####
  ## Relate spot areas with spots
  ## Calculate centroids and find closest spot to each centroid. This will be related to spot area
  
  centroids <- lapply(spot.owins, function(k) {
    suppressWarnings(spatstat::centroid.owin(k, as.ppp = TRUE)) # tmp <- 
    # spatstat::ppp(x = tmp$x, y = tmp$y, window = k)
  })
  
  spot.equivalent <- lapply(centroids, function(k) {
    dist2spots <- lapply(Z$spots, function(s) spatstat::nncross(k, s))
    min.dist <- dist2spots[which.min(do.call(rbind, dist2spots)$dist)]
    min.dist[[1]]$type <- names(min.dist)
    min.dist[[1]]$spot.mark <- as.integer(as.character(marks(Z$spots[[names(min.dist)]])[min.dist[[1]]$which]))
    min.dist[[1]]
  })
  
  G <- do.call(rbind, spot.equivalent)
  G$owin.name <- rownames(G)
  G <- G[with(G, order(type, spot.mark)),]
  
  out <- lapply(names(Z$spots), function(k) {
    tmp <- G[G$type %in% k,]$owin.name
    tmp.dt <- G[G$owin.name %in% tmp,]
    areas <- unlist(lapply(spot.owins[tmp], function(k) spatstat::area.owin(k)))
    diams <- unlist(lapply(spot.owins[tmp], function(k) spatstat::diameter(k)))
    spatstat::hyperframe(spot = tmp.dt$spot.mark, 
                         dist2spot = tmp.dt$dist, 
                         spot.owin.name = tmp.dt$owin.name, 
                         spot.owins = spot.owins[tmp], 
                         spot.area = areas, 
                         spot.diameter = diams)
  })
  
  names(out) <- names(Z$spots)
  
  centroids <- lapply(names(Z$spots), function(k) {
    tmp <- G[G$type %in% k,]$owin.name
    centroid.name <- G[G$owin.name %in% tmp,]$spot.mark
    tmp2 <- c()
    for(i in seq_along(tmp)){
      tmp3 <- centroids[[tmp[i]]]
      marks(tmp3) <- factor(centroid.name[i])
      tmp2 <- spatstat::superimpose(tmp2, tmp3, W = Z$window)
    }
    tmp2
  })
  
  names(centroids) <- names(Z$spots)  
  
  Z$spot.area <- list(centroids = centroids, spot.dat = out)
  
  class(Z) <- "rawDist"
  return(Z)
}