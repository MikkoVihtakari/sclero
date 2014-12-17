#' @title Plot spotDist object
#' @description Plots a map of \code{\link[=spot.dist]{spotDist}} object.
#' 
#' @param x \code{\link[=spot.dist]{spotDist}} object 
#' @param sample.name A character vector of length 1 specifying the sample name to be plotted as an overall title for the plot (\code{\link[=plot.default]{main}}). Defaults to \code{"keep"} meaning that the sample name will be extracted from the \code{\link[=spot.dist]{spotDist}} object.
#' @param crossing.points Indicates whether the crossing points between sampling spot sequence traverses and growth lines should be plotted. Defaults to FALSE.
#' @param size optional. A logical vector of length 1 indicating whether actual spot size should be used in plotting. Defaults to FALSE.
#' @param ... Arguments to be passed to other methods, such as \link[=par]{graphical parameters}. Not implemented yet.
#' @details Because the function plots a sample map, the aspect ratio is forced to 1 and cannot be changed. This might cause troubles when trying to set the axis limits. Try resizing the graphics window.
#' @author Mikko Vihtakari
#' @seealso \code{\link{spot.dist}} for aligning sample spots.
#' 
#' \code{\link{convert.ijdata}} for converting the coordinate information to \link[spatstat]{spatstat} point patterns. 
#' 
#' \code{\link{read.ijdata}} for reading zip files containing ImageJ ROIs.
#' 
#' \code{\link{plot.rawDist}} for plotting \code{\link[=convert.ijdata]{rawDist}} objects.
#' 
#' \code{\link[=plot]{plot.default}} and other methods; \code{\link[=points]{points}}, \code{\link[=lines]{lines}}, \code{\link[=par]{par}}.
#' 
#' @examples data(shellspots)
#' shell_map <- convert.ijdata(shellspots)
#' x <- spot.dist(shell_map)
#' plot(x) 
#' @import spatstat
#' @export

plot.spotDist <- function(x, ..., sample.name = 'keep', crossing.points = FALSE, size = FALSE){

## Debugging parameters, remove when ready
#X = x; main = "sample.name"

X <- x
  
## Define the unit for axes

unit <- strsplit(as.character(X$window$units), "/")[[1]]

## Spots

if(X$mid.point.type == "spots") {
  spots <- X$spots
} else {
if(X$mid.point.type == "centroids") {
  spots <- X$spot.area$centroids
}
}
  
## Drawing parameters
#if(class(xlim) == "function") xlim = X$window$xrange
#if(class(ylim) == "function") ylim = X$window$yrange
if(sample.name == "keep") sample.name = X$sample.name
  
## Establish colors for hole sequences
hole.cols <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999")

## Plot base plot
plot(x = X$window$xrange, y = X$window$yrange, type = "n", main = sample.name, asp = 1, axes = F, xlab = paste0("x ", "(", unit, ")"), ylab = paste0("y ", "(", unit, ")"))

## Plot growth bands
plot(X$gbs, col = "darkgrey", add = T)

## Plot the hole sequences
if(size) {
  lapply(X$spot.area$spot.dat, function(k) {
  lapply(k$spot.owins, function(j) plot(j, add = TRUE))})
lapply(seq_along(spots), function(i){
  text(coords(spots[[i]])[,1], coords(spots[[i]])[,2], marks(spots[[i]]), col = hole.cols[i], cex = 0.6, font = 2)})
} else {  
  lapply(seq_along(spots), function(i){
  points(coords(spots[[i]])[,1], coords(spots[[i]])[,2], cex = 2)
  text(coords(spots[[i]])[,1], coords(spots[[i]])[,2], marks(spots[[i]]), col = hole.cols[i], cex = 0.6, font = 2)})}

## Plot the main axis
plot(X$main, add = T, col = "darkslateblue")

## Plot the end and start positions for the main axis
text(coords(X$start.main)[,1], coords(X$start.main)[,2], X$start.main$marks, col = "dark green", font = 2)
text(coords(X$end.main)[,1], coords(X$end.main)[,2], X$end.main$marks, col = "dark green", font = 2)

## Plot the growth band projections on the main axis
text(coords(X$gb.projections)[,1], coords(X$gb.projections)[,2], X$gb.projections$marks, col = "blue", cex = 0.6, srt = 90, font = 2)

## Plot the hole sequence projections next to the main axis

print <- lapply(X$output, function(i) data.frame(hole = i$spot, dist = i$dist))
end <- X$gb.end
start <- X$gb.start
x.new <- lapply(print, function(i) (i[,2]/crossdist(start, end))*coords(end)$x + (1-(i[,2]/crossdist(start, end)))*coords(start)$x)
y.new <- lapply(print, function(i) (i[,2]/crossdist(start, end))*coords(end)$y + (1-(i[,2]/crossdist(start, end)))*coords(start)$y)

adjust <- lapply(seq_along(spots), function(i) c(0,-i+1))  

lapply(seq_along(x.new), function(i) text(x.new[[i]], y.new[[i]], print[[i]]$hole, adj=adjust[[i]], font = 2, col = hole.cols[i], cex = 0.6, offset = 0))

## Optional parameters
if(crossing.points == TRUE){
line.cross <- X$line.cross
lapply(seq_along(line.cross), function(i) points(coords(line.cross[[i]])[,1], coords(line.cross[[i]])[,2], cex = 1, pch = 19, col = hole.cols[i]))
  
cross.segs <- X$cross.segs
lapply(seq_along(cross.segs), function(x) plot(cross.segs[[x]], add = T, col = hole.cols[x]))}
  
axis(1)
axis(2, las = 2)}
