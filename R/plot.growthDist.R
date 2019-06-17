#' @title Plot a growthDist object
#' @description Plots a map of \code{\link[=growth.dist]{growthDist}} object.
#' @param X \code{\link[=growth.dist]{growthDist}} object
#' @param ... Arguments to be passed to other methods, such as \link[=par]{graphical parameters}.
#' @param growth.types A character vector specifying the names of growth axes to be plotted. If \code{NULL} (default), all growth axes within \code{X} will be plotted.
#' @param sample.name A character argument specifying the sample name to be plotted as an overall title for the plot (\code{\link[=plot.default]{main}}). Defaults to \code{"keep"} meaning that the sample name will be extracted from the \code{\link[=convert.ijdata]{rawDist}} object. The plot title can be omitted by specifying \code{sample.name = NULL}.
#' @param main.type A character argument with four possible levels (\code{"all"}, \code{"axis"}, \code{"ends"}, and \code{"none"}) indicating how the distance / main axis should be plotted. Defaults to \code{"all"} indicating that both the main axis and end points should be plotted. If \code{"axis"} only the main axis will be plotted. If \code{"ends"} only the end points will be plotted, and if \code{"none"} the main axis intormation will not be plotted.
#' @param highlight.gbs A character vector specifying the names of growth bands to be highlighted (i.e. colored with a different color than "darkgrey"). If \code{NULL} (default) all growth bands will be drawn using the standard color.
#' @param highlight.col A character argument specifying the color to be used in growth band highlighting (\code{highlight.gbs}).
#' @param color.palette Not implemented.
#' @details The \pkg{sclero} package currently uses the \pkg{graphics} package distributed with R for plotting (see \code{\link[graphics]{plot}}). Plotting sample maps is carried out by the \code{sclero:::samplemap} function, which works as an internal function and therefore has not been exported. Users willing to modify \pkg{sclero} plots beyond the flexibility allowed by the \code{\link{plot.rawDist}} function are instructed to modify the \code{\link{samplemap}} function, which consists of standard R graphics syntax. 
#'
#' Because the function plots a sample map, the \strong{aspect ratio} is forced to 1 and cannot be changed. If this causes troubles when trying to set the axis limits (\code{ylim} and \code{xlim}), try resizing the graphics window.
#' @author Mikko Vihtakari
#' @seealso \code{\link{growth.dist}} for calculation of growth distances.
#'
#' \code{\link[=plot]{plot.default}} and other methods; \code{\link[=points]{points}}, \code{\link[=lines]{lines}}, \code{\link[=par]{par}}.
#'
#' @import spatstat
#' @importFrom grDevices colorRampPalette
#' @importFrom stats na.omit
#' @importFrom graphics axis layout par plot points text
#' @export

# X = growthDist
# growth.types = NULL; sample.name = "keep"; main.type = "all"; color.palette = NULL; highlight.gbs = NULL; highlight.col = "red"
plot.growthDist <- function(X, ..., growth.types = NULL, sample.name = "keep", main.type = "all", color.palette = NULL, highlight.gbs = NULL, highlight.col = "red"){
  
  ## Plot the sample map
  
  output <- samplemap(x = X$rawDist, ..., sname = sample.name, sptype = "id", size = 2, scol = NULL, mtype = main.type, gaxis = NULL, colpalette = NULL, hlight = highlight.gbs, hlcol = highlight.col)
  
  ## Growth axes
  gaxes <- Filter(Negate(is.null), X$growthDist$growth.axes)
  
  if(!is.null(growth.types)) {
    gaxes <- gaxes[growth.types]
  }
  
  lapply(seq_along(gaxes), function(i) {
    plot(gaxes[[i]], add = TRUE, col = output$hole.cols[i])
    
    ps <- lapply(unique(marks(gaxes[[i]])), function(k) {
      tmp <- spatstat::endpoints.psp(gaxes[[i]][marks(gaxes[[i]]) == k])
      tmp[tmp$n]
    })
    
    ps <- do.call(spatstat::superimpose, ps)
    ps <- superimpose(endpoints.psp(gaxes[[i]])[1], ps)
    plot(ps, add = TRUE, col = output$hole.cols[i])
  })
  
  legend("topright", legend = names(gaxes), fill = output$hole.cols[1:length(gaxes)], title = "Growth axis", bty = "n")
  
}
