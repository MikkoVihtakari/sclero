#' @title Plot rawDist object
#' @description Plots a map of \code{\link[=convert.ijdata]{rawDist}} object.
#' 
#' @param x rawDist object
#' @param sample.name A character vector of length 1 specifying the sample name to be plotted as an overall title for the plot (\code{\link[=plot.default]{main}}). Defaults to \code{"keep"} meaning that the sample name will be extracted from the \code{\link[=convert.ijdata]{rawDist}} object.
#' @param spot.type levels "id", "value", "idvalue".
#' @param spot.size cex parameter at the moment
#' @param color.palette color palette for "value" and "idvalue". Passed to \code{\link{colorRampPalette}}.
#' @param ... Arguments to be passed to other methods, such as \link[=par]{graphical parameters}. Not implemented yet.
#' @details Because the function plots a sample map, the aspect ratio is forced to 1 and cannot be changed. This might cause troubles when trying to set the axis limits. Try resizing the graphics window.
#' @author Mikko Vihtakari
#' @seealso \code{\link{convert.ijdata}} for converting the coordinate information to \link[spatstat]{spatstat} point patterns. 
#' 
#' \code{\link{read.ijdata}} for reading zip files containing ImageJ ROIs.
#' 
#' \code{\link{plot.spotDist}} for plotting \code{\link[=spot.dist]{spotDist}} objects.
#' 
#' \code{\link[=plot]{plot.default}} and other methods; \code{\link[=points]{points}}, \code{\link[=lines]{lines}}, \code{\link[=par]{par}}.
#' 
#' @examples data(shellspots)
#' shell_map <- convert.ijdata(shellspots)
#' plot(shell_map)
#' 
#' @import spatstat
#' @export

plot.rawDist <- function(x, ..., sample.name = "keep", spot.type = "id", spot.size = 2, color.palette = colorRampPalette(c("blue", "cyan", "yellow", "red"), bias=1)(100)){

## Debugging parameters, remove when ready
#X = PT; sample.name = "keep"; spot.type = "idvalue"; spot.size = 2; legend = TRUE; color.palette = colorRampPalette(c("blue", "cyan", "yellow", "red"), bias=1)(100)

## Define the unit for axes
X <- x
unit <- strsplit(as.character(X$window$units), "/")[[1]]

## Drawing parameters
#if(class(xlim) == "function") xlim = X$window$xrange
#if(class(ylim) == "function") ylim = X$window$yrange
if(sample.name == "keep") sample.name = X$sample.name
  
## Establish colors for hole sequences
hole.cols <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999")

## Make layout if spot.type == "value" | spot.type == "idvalue"
  
if(spot.type == "value" | spot.type == "idvalue") layout(matrix(c(1,2), nrow = 1), widths = c(9,1))
  
## Plot base plot
par(mar = c(4, 4, 2.5, 1) + 0.1)
plot(x = X$window$xrange, y = X$window$yrange, type = "n", main = sample.name, asp = 1, axes = F, xlab = paste0("x ", "(", unit, ")"), ylab = paste0("y ", "(", unit, ")"))#, ...), xlim = xlim, ylim = ylim

## Plot growth bands
plot(X$gbs, col = "darkgrey", add = T)

## Plot the hole sequences
if(spot.type == "id") {
lapply(seq_along(X$spots), function(i){
  points(coords(X$spots[[i]])[,1], coords(X$spots[[i]])[,2], cex = 2)
  text(coords(X$spots[[i]])[,1], coords(X$spots[[i]])[,2], marks(X$spots[[i]]), col = hole.cols[i], cex = 0.6, font = 2)})
  } else {
    if(spot.type == "value") {
      if(all(unlist(lapply(X$values, is.null)))) stop("rawDist values must be specified")
      colmap <- colourmap(color.palette, range = range(unlist(lapply(X$value, function(k) range(k[,2]))), na.rm = TRUE))
      lapply(seq_along(X$spots), function(i){
      color <- ifelse(is.na(X$values[[i]][,2]), "white", colmap(X$values[[i]][,2]))
      plot(X$spots[[i]], add = TRUE, use.marks = FALSE, cex = spot.size, bg = color, pch = 21)})
    } else {
     if(spot.type == "idvalue") {
      if(all(unlist(lapply(X$values, is.null)))) stop("rawDist values must be specified")
      colmap <- colourmap(color.palette, range = range(unlist(lapply(X$value, function(k) range(k[,2]))), na.rm = TRUE))
      lapply(seq_along(X$spots), function(i){
      color <- ifelse(is.na(X$values[[i]][,2]), "white", colmap(X$values[[i]][,2]))
      plot(X$spots[[i]], add = TRUE, use.marks = FALSE, cex = spot.size, bg = color, pch = 21)
      text(coords(X$spots[[i]])[,1], coords(X$spots[[i]])[,2], marks(X$spots[[i]]), col = "black", cex = 0.4, font = 2
      )}) } else {
      stop("Invalid spot.type. Use 'identity' or 'value'.")
    }}}
  
## Plot the main axis
plot(X$main, add = T, col = "darkslateblue")

## Plot the end and start positions for the main axis
text(coords(X$start.main)[,1], coords(X$start.main)[,2], X$start.main$marks, col = "dark green", font = 2)
text(coords(X$end.main)[,1], coords(X$end.main)[,2], X$end.main$marks, col = "dark green", font = 2)

## Plot the axes
  
axis(1)
axis(2, las = 2)

## Add legend
  
if(spot.type == "value" | spot.type == "idvalue") {
par(mar=c(1,0.5,3,2))
   plot(colmap, vertical = TRUE, las = 2, main = expression(paste(delta^{18},"O", sep = "")), cex = 0.1)}
}

