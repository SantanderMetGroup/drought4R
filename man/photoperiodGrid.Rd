% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/photoperiodGrid.R
\name{photoperiodGrid}
\alias{photoperiodGrid}
\title{Photoperiod}
\usage{
photoperiodGrid(c4r.obj)
}
\arguments{
\item{c4r.obj}{climate4R object. Note that this is independent of the variable stored, as it will use only the calendar
an geographical information of the object to compute the photoperiod}
}
\description{
Compute photoperiod (daylength) as a funtion of latitude and date
}
\details{
This is just a wrapper from the \code{\link[geosphere]{daylength}} function of package \pkg{geosphere}.
}
\examples{
require(transformeR)
data("tasmin.eobs.iberia.daily")
daylength <- photoperiodGrid(tasmin.eobs.iberia.daily)
lat <- range(getCoordinates(daylength)[["y"]])
lat.ranges <- seq(from = lat[1], to = tail(lat, 1), length.out = 5)
lat.list <- lapply(1:(length(lat.ranges) - 1), function(x) {
    subsetGrid(daylength, latLim = c(lat.ranges[x], lat.ranges[x + 1]), outside = TRUE)
})
require(visualizeR)
temporalPlot(lat.list, xyplot.custom = list(main = "Photoperiod as a function of latitude",
                                            ylab = "Daylight hours",
                                            xlab = "Date",
                                            key = list(corner = c(1,.5),
                                                       lines = list(col = 1:4),
                                                       horizontal = TRUE,
                                                       text = list(c("37ºN","39ºN","41ºN","43ºN")))))
}
\author{
J. Bedia
}
