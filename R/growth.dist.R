#' @title Calculate growth using a variety of measurement techniques
#' @description Calculates growth from rawDist \code{\link[=convert.ijdata]{rawDist}} objects and provides a range of results depending on the type of the sample.
#' @param rawDist a \code{\link[=convert.ijdata]{rawDist}} object for which the alignment should be done.
#' @param coord.type A character or numeric argument defining the coordinate type to do the calculations on. Alternatives:
#' \itemize{
#'   \item \strong{"scaled"} or \strong{1}. Only scaled coordinates.
#'   \item \strong{"original"} or any other character/number. Only original ImageJ ROI coordinates.
#' }
#' @details This function is still under developemt and may not work as intended.
#' @return A list of class \code{"growthDist"}. The first element (\code{$data}) contains a data.frame of growth measurements. The \code{$id} column gives the order of growth increments from the start point defined by the \code{main} axis. The \code{$gap} column gives the name of growth lines between which the growth is measured. Following columns give the growth measurements in units defined by the user in \code{\link{read.ijdata}} function (see the \code{scale} and \code{unit} arguments). The columns denoted as \code{.pr} give the percentage growth as compared to estimated total growth in the sample. The second list element contains \link[spatstat]{spatstat} segmented line (\code{psp}) patterns of the various growth measurement methods. List elements with \code{NULL} values were not applicable for the sample type. Abbreviations for measurement types are given in parenthesis. The types are:
#' \itemize{
#'   \item \code{\strong{$main}} (main).
#'   \item \code{\strong{$manual}} (man). Only applicable for \code{cross} section types.
#'   \item \code{\strong{$caliber}} (cal). Only applicable for \code{along} section types.
#'   \item \code{\strong{$maximum}} (max).
#'   \item \code{\strong{$maximum.along.x1}} (maxx1). Only applicable for \code{cross} section types.
#'   \item \code{\strong{$maximum.along.x2}} (maxx2). Only applicable for \code{cross} section types.
#'   \item \code{\strong{$guided.maximum}} (maxg). Only applicable for \code{cross} section types.
#'   \item \code{\strong{$minimum}} (min). Only applicable for \code{cross} section types.
#'   \item \code{\strong{$guided.minimum}} (ming). Only applicable for \code{cross} section types.
#'   \item \code{\strong{$direct}} (dire). Only applicable for \code{along} section types.
#' }
#' @seealso \code{\link{read.ijdata}} for reading zip files containing ImageJ ROIs.
#' \code{\link{convert.ijdata}} for converting the coordinate information to \link[spatstat]{spatstat} point patterns.
#' @import spatstat plyr
#' @export

# sample.name = NULL; coord.type = "scaled"
growth.dist <- function(rawDist, coord.type = "scaled") {

  ## Read in the data ####

  if(coord.type == "scaled" | coord.type == 1) {
    x <- rawDist$scaled
  } else {
    x <- rawDist$original
  }

  ## Part 1. Define parameters

  gbs <- x$gbs
  main <- x$main
  start <- x$start.main
  end <- x$end.main
  growth <- x$growth
  lines <- unique(gbs$marks)
  nlines <- length(lines)
  window <- x$window
  gbs.all <- spatstat::superimpose(gbs)

  ## Part 1.2 Test whether the object is 'along' or 'cross' type (see the tutorial) ###

  test <- crossing.psp(main, superimpose(gbs))

  if(test$n == 0) {

    sec.type <- "along"

  } else {

    if(test$n == length(unique(gbs$marks))) {

      sec.type <- "cross"
      if(!(coord.type == "scaled" | coord.type == 1)) message("Using the original coordinates disables maximum distance based methods.")

    } else {

      stop(paste0("Number of growth line and main axis crossing points is neither 0 or ", length(unique(gbs$marks)), ". Cannot define the aligment type. See ?spot.dist"))

    }
  }

  ## Part 1.3 Define start and endpoint for measurements

  if(sec.type == "cross") {

    startp <- crossing.psp(main, gbs.all)[1]
    marks(startp) <- "start point"

    endp <- crossing.psp(main, gbs.all)[crossing.psp(main, gbs.all)$n]
    marks(endp) <- "end point"

  } else {

    tmp <- gbs.all[marks(gbs.all) == lines[1]][1]$ends
    startp <- ppp(x = tmp$x0, y = tmp$y0, window = window, marks = "start point")
    startp.main <- project2segment(startp, main)$Xproj

    tmp <- gbs.all[marks(gbs.all) == lines[nlines]]
    tmp <- tmp[tmp$n]$ends
    endp <- ppp(x = tmp$x1, y = tmp$y1, window = window, marks = "end point")
    endp.main <- project2segment(endp, main)$Xproj

  }


  ## Part 2 Manual growth axis ####

  if(sec.type == "cross") {

    tmp <- spatstat::crossing.psp(growth, gbs.all, details = TRUE)
    tmp.marks <- spatstat::marks(gbs.all)[spatstat::marks(tmp)$jB]
    tmp.marks <- paste(tmp.marks[1:(length(tmp.marks)-1)], tmp.marks[2:length(tmp.marks)], sep = "-")
    tmp.marks <- factor(tmp.marks, tmp.marks)

    growth.man <- spatstat::psp(x0 = tmp$x[1:(length(tmp$x)-1)],
      x1 = tmp$x[2:length(tmp$x)],
      y0 = tmp$y[1:(length(tmp$y)-1)],
      y1 = tmp$y[2:length(tmp$y)],
      marks = tmp.marks,
      window = window)

    dat.man <- make.growth.data(growth.man)

  }

  ## Part 3 Main axis & caliber for along types

  if(sec.type == "cross") {

    tmp <- spatstat::crossing.psp(main, gbs.all, details = TRUE)
    tmp.marks <- unique(spatstat::marks(gbs.all))
    tmp.marks <- paste(tmp.marks[1:(length(tmp.marks)-1)], tmp.marks[2:length(tmp.marks)], sep = "-")
    tmp.marks <- factor(tmp.marks, tmp.marks)

    growth.main <- spatstat::psp(x0 = tmp$x[1:(length(tmp$x)-1)],
      x1 = tmp$x[2:length(tmp$x)],
      y0 = tmp$y[1:(length(tmp$y)-1)],
      y1 = tmp$y[2:length(tmp$y)],
      marks = tmp.marks,
      window = window)

    dat.main <- make.growth.data(growth.main)

  } else if (sec.type == "along") { # along types

    #  k <- lines[6]
    meas.points <- lapply(lines, function(k) {
      tp <- gbs.all[marks(gbs.all) == k]
      st <- ppp(x = tp[1]$ends$x0, y = tp[1]$ends$y0, window = window, marks = k)
      ed <- ppp(x = tp[tp$n]$ends$x1, y = tp[tp$n]$ends$y1, window = window, marks = k)
      cal <- psp(x0 = startp$x, y0 = startp$y, x1 = ed$x, y1 = ed$y, window = window, marks = k)
      dire <- psp(x0 = st$x, y0 = st$y, x1 = ed$x, y1 = ed$y, window = window, marks = k)
      list(end.point = ed, cal = cal, dire = dire)
    })

    ps <- do.call(spatstat::superimpose, lapply(meas.points, function(k) k$end.point))
    tmp <- project2segment(ps, main)$Xproj

    if(any(tmp$x < startp.main$x)) {
      warning("Projected growth line end points have smaller values than the starting point for the 'main' (projected) method. The growth.main results are biased.")
    }

    l1 <- psp(x0 = startp.main$x, y0 = startp.main$y, x1 = tmp$x[1], y1 = tmp$y[1], window = window, marks = lines[1])

    l2 <- spatstat::psp(x0 = tmp$x[1:(length(tmp$x)-1)],
      x1 = tmp$x[2:length(tmp$x)],
      y0 = tmp$y[1:(length(tmp$y)-1)],
      y1 = tmp$y[2:length(tmp$y)],
      marks = lines[-1],
      window = window)

    growth.main <- spatstat::superimpose(l1, l2)
    growth.dire <- do.call(spatstat::superimpose, lapply(meas.points, function(k) k$dire))
    growth.cal <- do.call(spatstat::superimpose, lapply(meas.points, function(k) k$cal))

    dat.main <- make.growth.data(growth.main)
    dat.dire <- make.growth.data(growth.dire)

    dat.cal <- make.growth.data(growth.cal)
    dat.cal$cal <- c(dat.cal$cal[1], diff(dat.cal$cal))
    dat.cal$cal.pr <- 100*dat.cal$cal/sum(dat.cal$cal)
  }

  ## Part 4 Minimum distance ####

  if(sec.type == "cross") {

    tmp <- lapply(1:(nlevels(spatstat::marks(gbs.all))-1), function(i) {

      gbs1 <- gbs.all[spatstat::marks(gbs.all) == levels(spatstat::marks(gbs.all))[i]]
      gbs2 <- gbs.all[spatstat::marks(gbs.all) == levels(spatstat::marks(gbs.all))[i+1]]

      gbs1.mid <- spatstat::midpoints.psp(gbs1)
      spatstat::marks(gbs1.mid) <- unique(spatstat::marks(gbs1))

      gbs1.ppp <- suppressWarnings(spatstat::superimpose(gbs1.mid, spatstat::as.ppp(gbs1)))
      gbs1.ppp <- gbs1.ppp[!spatstat::duplicated.ppp(gbs1.ppp)]

      proj1 <- spatstat::project2segment(gbs1.ppp, gbs2)

      p1 <- gbs1.ppp[which.min(proj1$d)]
      p2 <- proj1$Xproj[which.min(proj1$d)]

      tmp.marks <- paste(unique(spatstat::marks(gbs1)), unique(spatstat::marks(gbs2)), sep = "-")

      growth.out <- spatstat::psp(x0 = p1$x,
        x1 = p2$x,
        y0 = p1$y,
        y1 = p2$y,
        marks = tmp.marks,
        window = window)

      list(line = growth.out, dist = min(proj1$d))

    })

    growth.min <- do.call(spatstat::superimpose, lapply(tmp, function(k) k$line))

    dat.min <- make.growth.data(growth.min)

  }

  ## Part 5 Guided minimum distance ####

  if(sec.type == "cross" & (coord.type == "scaled" | coord.type == 1)) {
    for(i in 2:(nlevels(spatstat::marks(gbs.all)))) {

      gbs1 <- gbs.all[spatstat::marks(gbs.all) == levels(spatstat::marks(gbs.all))[i-1]]
      gbs2 <- gbs.all[spatstat::marks(gbs.all) == levels(spatstat::marks(gbs.all))[i]]

      p1min <- spatstat::crossing.psp(spatstat::clip.infline(spatstat::infline(h = min(gbs2$ends[c("y0", "y1")])), win = window), gbs1)
      p1max <- spatstat::crossing.psp(spatstat::clip.infline(spatstat::infline(h = max(gbs2$ends[c("y0", "y1")])), win = window), gbs1)

      p2min <- spatstat::crossing.psp(spatstat::clip.infline(spatstat::infline(h = min(gbs1$ends[c("y0", "y1")])), win = window), gbs2)
      p2max <- spatstat::crossing.psp(spatstat::clip.infline(spatstat::infline(h = max(gbs1$ends[c("y0", "y1")])), win = window), gbs2)

      gbs1.mid <- spatstat::midpoints.psp(gbs1)
      gbs1.ppp <- spatstat::as.ppp(gbs1)
      marks(gbs1.ppp) <- NULL
      gbs1.ppp <- suppressWarnings(spatstat::superimpose(gbs1.mid, gbs1.ppp))
      gbs1.ppp <- gbs1.ppp[!spatstat::duplicated.ppp(gbs1.ppp)]

      gbs2.mid <- spatstat::midpoints.psp(gbs2)
      gbs2.ppp <- spatstat::as.ppp(gbs2)
      marks(gbs2.ppp) <- NULL
      gbs2.ppp <- suppressWarnings(spatstat::superimpose(gbs2.mid, gbs2.ppp))
      gbs2.ppp <- gbs2.ppp[!spatstat::duplicated.ppp(gbs2.ppp)]

      ## Cut the sequences ###

      if(p1min$n == 1 & p1max$n == 1) {
        gbs1.ppp <- gbs1.ppp[gbs1.ppp$y >= p1min$y & gbs1.ppp$y <= p1max$y]
        gbs1.ppp <- spatstat::superimpose(p1min, gbs1.ppp, p1max)
      } else if(p1min$n == 1) {
        gbs1.ppp <- gbs1.ppp[gbs1.ppp$y >= p1min$y]
        gbs1.ppp <- spatstat::superimpose(p1min, gbs1.ppp)
      } else if(p1max$n == 1) {
        gbs1.ppp <- gbs1.ppp[gbs1.ppp$y <= p1max$y]
        gbs1.ppp <- spatstat::superimpose(gbs1.ppp, p1max)
      }

      if(p2min$n == 1 & p2max$n == 1) {
        gbs2.ppp <- gbs2.ppp[gbs2.ppp$y >= p2min$y & gbs2.ppp$y <= p2max$y]
        gbs2.ppp <- spatstat::superimpose(p2min, gbs2.ppp, p2max)
      } else if(p2min$n == 1) {
        gbs2.ppp <- gbs2.ppp[gbs2.ppp$y >= p2min$y]
        gbs2.ppp <- spatstat::superimpose(p2min, gbs2.ppp)
      } else if(p2max$n == 1) {
        gbs2.ppp <- gbs2.ppp[gbs2.ppp$y <= p2max$y]
        gbs2.ppp <- spatstat::superimpose(gbs2.ppp, p2max)
      }

      ## Remake gbs1 and gbs2 equally long ###

      tmp <- data.frame(x = gbs1.ppp$x, y = gbs1.ppp$y)
      tmp <- tmp[order(tmp$y),]

      gbs1 <- spatstat::psp(x0 = tmp$x[1:(length(tmp$x)-1)],
        x1 = tmp$x[2:length(tmp$x)],
        y0 = tmp$y[1:(length(tmp$y)-1)],
        y1 = tmp$y[2:length(tmp$y)],
        marks = rep(unique(marks(gbs1)), nrow(tmp)-1),
        window = window)

      tmp <- data.frame(x = gbs2.ppp$x, y = gbs2.ppp$y)
      tmp <- tmp[order(tmp$y),]

      gbs2 <- spatstat::psp(x0 = tmp$x[1:(length(tmp$x)-1)],
        x1 = tmp$x[2:length(tmp$x)],
        y0 = tmp$y[1:(length(tmp$y)-1)],
        y1 = tmp$y[2:length(tmp$y)],
        marks = rep(unique(marks(gbs2)), nrow(tmp)-1),
        window = window)

      ## Point calculation

      if(i == 2) {
        p1 <- startp
        p2 <- project2segment(p1, gbs2)$Xproj

        line <- spatstat::psp(x0 = p1$x,
          x1 = p2$x,
          y0 = p1$y,
          y1 = p2$y,
          marks = paste(unique(spatstat::marks(gbs1)), unique(spatstat::marks(gbs2)), sep = "-"),
          window = window)

      } else if(i == nlevels(spatstat::marks(gbs.all))) {
        p1 <- p2
        p2 <- endp

        addline <- spatstat::psp(x0 = p1$x,
          x1 = p2$x,
          y0 = p1$y,
          y1 = p2$y,
          marks = paste(unique(spatstat::marks(gbs1)), unique(spatstat::marks(gbs2)), sep = "-"),
          window = window)

        line <- append.psp(line, addline)

      } else {
        p1 <- p2
        p2 <- project2segment(p1, gbs2)$Xproj

        addline <- spatstat::psp(x0 = p1$x,
          x1 = p2$x,
          y0 = p1$y,
          y1 = p2$y,
          marks = paste(unique(spatstat::marks(gbs1)), unique(spatstat::marks(gbs2)), sep = "-"),
          window = window)

        line <- append.psp(line, addline)
      }
    }

    growth.ming <- line

    dat.ming <- make.growth.data(growth.ming)


  }

  ## Part 6 Maximum distances ####

  if(sec.type == "cross" & (coord.type == "scaled" | coord.type == 1)) {

    # i = 1
    tmp <- lapply(1:(nlevels(spatstat::marks(gbs.all))-1), function(i) {
      ## ####
      gbs1 <- gbs.all[spatstat::marks(gbs.all) == levels(spatstat::marks(gbs.all))[i]]
      gbs2 <- gbs.all[spatstat::marks(gbs.all) == levels(spatstat::marks(gbs.all))[i+1]]

      p1min <- spatstat::crossing.psp(spatstat::clip.infline(spatstat::infline(h = min(gbs2$ends[c("y0", "y1")])), win = window), gbs1)
      p1max <- spatstat::crossing.psp(spatstat::clip.infline(spatstat::infline(h = max(gbs2$ends[c("y0", "y1")])), win = window), gbs1)

      p2min <- spatstat::crossing.psp(spatstat::clip.infline(spatstat::infline(h = min(gbs1$ends[c("y0", "y1")])), win = window), gbs2)
      p2max <- spatstat::crossing.psp(spatstat::clip.infline(spatstat::infline(h = max(gbs1$ends[c("y0", "y1")])), win = window), gbs2)

      gbs1.mid <- spatstat::midpoints.psp(gbs1)
      gbs1.ppp <- spatstat::as.ppp(gbs1)
      marks(gbs1.ppp) <- NULL
      gbs1.ppp <- suppressWarnings(spatstat::superimpose(gbs1.mid, gbs1.ppp))
      gbs1.ppp <- gbs1.ppp[!spatstat::duplicated.ppp(gbs1.ppp)]

      gbs2.mid <- spatstat::midpoints.psp(gbs2)
      gbs2.ppp <- spatstat::as.ppp(gbs2)
      marks(gbs2.ppp) <- NULL
      gbs2.ppp <- suppressWarnings(spatstat::superimpose(gbs2.mid, gbs2.ppp))
      gbs2.ppp <- gbs2.ppp[!spatstat::duplicated.ppp(gbs2.ppp)]

      ## Cut the sequences ###

      if(p1min$n == 1 & p1max$n == 1) {
        gbs1.ppp <- gbs1.ppp[gbs1.ppp$y >= p1min$y & gbs1.ppp$y <= p1max$y]
        gbs1.ppp <- spatstat::superimpose(p1min, gbs1.ppp, p1max)
      } else if(p1min$n == 1) {
        gbs1.ppp <- gbs1.ppp[gbs1.ppp$y >= p1min$y]
        gbs1.ppp <- spatstat::superimpose(p1min, gbs1.ppp)
      } else if(p1max$n == 1) {
        gbs1.ppp <- gbs1.ppp[gbs1.ppp$y <= p1max$y]
        gbs1.ppp <- spatstat::superimpose(gbs1.ppp, p1max)
      }

      if(p2min$n == 1 & p2max$n == 1) {
        gbs2.ppp <- gbs2.ppp[gbs2.ppp$y >= p2min$y & gbs2.ppp$y <= p2max$y]
        gbs2.ppp <- spatstat::superimpose(p2min, gbs2.ppp, p2max)
      } else if(p2min$n == 1) {
        gbs2.ppp <- gbs2.ppp[gbs2.ppp$y >= p2min$y]
        gbs2.ppp <- spatstat::superimpose(p2min, gbs2.ppp)
      } else if(p2max$n == 1) {
        gbs2.ppp <- gbs2.ppp[gbs2.ppp$y <= p2max$y]
        gbs2.ppp <- spatstat::superimpose(gbs2.ppp, p2max)
      }

      ## Remake gbs1 and gbs2 equally long ###

      tmp <- data.frame(x = gbs1.ppp$x, y = gbs1.ppp$y)
      tmp <- tmp[order(tmp$y),]

      gbs1 <- spatstat::psp(x0 = tmp$x[1:(length(tmp$x)-1)],
        x1 = tmp$x[2:length(tmp$x)],
        y0 = tmp$y[1:(length(tmp$y)-1)],
        y1 = tmp$y[2:length(tmp$y)],
        marks = rep(unique(marks(gbs1)), nrow(tmp)-1),
        window = window)

      tmp <- data.frame(x = gbs2.ppp$x, y = gbs2.ppp$y)
      tmp <- tmp[order(tmp$y),]

      gbs2 <- spatstat::psp(x0 = tmp$x[1:(length(tmp$x)-1)],
        x1 = tmp$x[2:length(tmp$x)],
        y0 = tmp$y[1:(length(tmp$y)-1)],
        y1 = tmp$y[2:length(tmp$y)],
        marks = rep(unique(marks(gbs2)), nrow(tmp)-1),
        window = window)

      gbs1.ppp <- pointsOnLines(gbs1)
      gbs2.ppp <- pointsOnLines(gbs2)

      ## Distance calculations ###

      ## Maximum along x-axis alt 1 ###

      p1 <- gbs1.ppp
      p2 <- pointsAlongLines(gbs2, y = p1)

      line.alts <- spatstat::psp(x0 = p1$x,
        x1 = p2$x,
        y0 = p1$y,
        y1 = p2$y,
        window = window)

      maxx1 <- line.alts[which.max(spatstat::lengths.psp(line.alts))]
      marks(maxx1) <- paste(unique(spatstat::marks(gbs1)), unique(spatstat::marks(gbs2)), sep = "-")

      ## Maximum along x-axis alt 2 ###

      dist2gbs1 <- distfun(gbs1.ppp)

      maxp.gbs2 <- gbs2.ppp[which.max(dist2gbs1(x = gbs2.ppp$x, y = gbs2.ppp$y))]
      p1 <- pointsAlongLines(gbs1, y = maxp.gbs2)

      maxx2 <- spatstat::psp(x0 = p1$x,
        x1 = maxp.gbs2$x,
        y0 = p1$y,
        y1 = maxp.gbs2$y,
        marks = paste(unique(spatstat::marks(gbs1)), unique(spatstat::marks(gbs2)), sep = "-"),
        window = window)

      ## Maximum distance with closest point along gbs1

      maxp.gbs1 <- project2segment(maxp.gbs2, gbs1)$Xproj

      maxcp <- spatstat::psp(x0 = maxp.gbs1$x,
        x1 = maxp.gbs2$x,
        y0 = maxp.gbs1$y,
        y1 = maxp.gbs2$y,
        marks = paste(unique(spatstat::marks(gbs1)), unique(spatstat::marks(gbs2)), sep = "-"),
        window = window)

      ## Guided maximum distance from start point

      if(i == 1) {
        dist2sp <- distfun(startp)
        maxp.st.gbs2 <- gbs2.ppp[which.max(dist2sp(x = gbs2.ppp$x, y = gbs2.ppp$y))]

        maxg <- spatstat::psp(x0 = startp$x,
          x1 = maxp.st.gbs2$x,
          y0 = startp$y,
          y1 = maxp.st.gbs2$y,
          marks = paste(unique(spatstat::marks(gbs1)), unique(spatstat::marks(gbs2)), sep = "-"),
          window = window)
      } else if(i == nlevels(spatstat::marks(gbs.all))-1) {
        dist2ep <- distfun(endp)
        maxp.st.gbs1 <- gbs1.ppp[which.max(dist2ep(x = gbs1.ppp$x, y = gbs1.ppp$y))]

        maxg <- spatstat::psp(x0 = maxp.st.gbs1$x,
          x1 = endp$x,
          y0 = maxp.st.gbs1$y,
          y1 = endp$y,
          marks = paste(unique(spatstat::marks(gbs1)), unique(spatstat::marks(gbs2)), sep = "-"),
          window = window)
      } else {
        maxp.gbs1 <- project2segment(maxp.gbs2, gbs1)$Xproj

        maxg <- spatstat::psp(x0 = maxp.gbs1$x,
          x1 = maxp.gbs2$x,
          y0 = maxp.gbs1$y,
          y1 = maxp.gbs2$y,
          marks = paste(unique(spatstat::marks(gbs1)), unique(spatstat::marks(gbs2)), sep = "-"),
          window = window)
      }

      ## ####

      list(maxx1 = maxx1, maxx2 = maxx2, maxcp = maxcp, maxg = maxg)
    })

    growth.maxx1 <- do.call(spatstat::superimpose, lapply(tmp, function(k) k$maxx1))
    growth.maxx2 <- do.call(spatstat::superimpose, lapply(tmp, function(k) k$maxx2))
    growth.max <- do.call(spatstat::superimpose, lapply(tmp, function(k) k$maxcp))
    growth.maxg <- do.call(spatstat::superimpose, lapply(tmp, function(k) k$maxg))

    dat.maxx1 <- make.growth.data(growth.maxx1)
    dat.maxx2 <- make.growth.data(growth.maxx2)
    dat.max <- make.growth.data(growth.max)
    dat.maxg <- make.growth.data(growth.maxg)


  } else if(sec.type == "along" & x$parameters$gbs.rois == "polyline") { ## ALONG TYPE

  tmp <- data.frame(gap = marks(gbs.all), dist = lengths.psp(gbs.all))
  tmp <- ddply(tmp, .(gap), summarise, max = sum(dist))
  tmp <- cbind(data.frame(id = 1:nrow(tmp)), tmp)
  tmp$max.pr <- 100*tmp$max/sum(tmp$max)

  dat.max <- tmp
  growth.max <- gbs.all

  }

  ######################
  ## Part 7. Return ####


  possible.data <- c("dat.main", "dat.man", "dat.cal", "dat.max", "dat.maxg", "dat.maxx1", "dat.maxx2", "dat.min", "dat.ming", "dat.dire")
  available.data <- possible.data[possible.data %in% ls(envir = environment())]
  all.data <- lapply(available.data, function(k) {get(k, envir = environment())})
  dat <- Reduce(function(...) merge(..., all=TRUE), all.data)

  dat <- dat[c(names(dat)[!grepl(".pr", names(dat))], names(dat)[grepl(".pr", names(dat))])]
  dat <- dat[order(dat$id),]
  rownames(dat) <- 1:nrow(dat)

  out <- list(data = dat, growth.axes = list(
    main = growth.main,
    manual = get0("growth.man"),
    caliber = get0("growth.cal"),
    maximum = get0("growth.max"),
    guided.maximum = get0("growth.maxg"),
    maximum.along.x1 = get0("growth.maxx1"),
    maximum.along.x2 = get0("growth.maxx2"),
    minimum = get0("growth.min"),
    guided.minimum = get0("growth.ming"),
    direct = get0("growth.dire")
    ))

  class(out) <- "growthDist"
  return(out)
}

############################
## Return growth data.frame
# X <- growth.min
make.growth.data <- function(X) {
  var <- deparse(substitute(X))

  out <- data.frame(id = 1:X$n,
    gap = spatstat::marks(X),
    sub1 = spatstat::lengths.psp(X),
    sub2 = 100*spatstat::lengths.psp(X)/sum(spatstat::lengths.psp(X)))

  names(out)[names(out) == "sub1"] <- select(strsplit(var, "\\."), 2)
  names(out)[names(out) == "sub2"] <- paste(select(strsplit(var, "\\."), 2), "pr", sep = ".")
  return(out)
}
