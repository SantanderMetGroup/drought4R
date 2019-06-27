##     photoperiod.R Compute photoperiod (daylength) as a funtion of latitude and date
##
##     Copyright (C) 2019 Santander Meteorology Group (http://www.meteo.unican.es)
##
##     This program is free software: you can redistribute it and/or modify
##     it under the terms of the GNU General Public License as published by
##     the Free Software Foundation, either version 3 of the License, or
##     (at your option) any later version.
## 
##     This program is distributed in the hope that it will be useful,
##     but WITHOUT ANY WARRANTY; without even the implied warranty of
##     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##     GNU General Public License for more details.
## 
##     You should have received a copy of the GNU General Public License
##     along with this program.  If not, see <http://www.gnu.org/licenses/>.

#' @title Photoperiod
#' @description Compute photoperiod (daylength) as a funtion of latitude and date
#' @param c4r.obj climate4R object. Note that this is independent of the variable stored, as it will use only the calendar
#' an geographical information of the object to compute the photoperiod 
#' @importFrom geosphere daylength
#' @importFrom transformeR getRefDates getCoordinates getShape redim subsetGrid mat2Dto3Darray bindGrid
#' @importFrom magrittr %>% %<>% 
#' @details This is just a wrapper from the \code{\link[geosphere]{daylength}} function of package \pkg{geosphere}.
#' @export
#' @author J. Bedia
#' @examples 
#' require(transformeR)
#' data("tasmin.eobs.iberia.daily")
#' daylength <- photoperiodGrid(tasmin.eobs.iberia.daily)
#' lat <- range(getCoordinates(daylength)[["y"]])
#' lat.ranges <- seq(from = lat[1], to = tail(lat, 1), length.out = 5)
#' lat.list <- lapply(1:(length(lat.ranges) - 1), function(x) {
#'     subsetGrid(daylength, latLim = c(lat.ranges[x], lat.ranges[x + 1]), outside = TRUE)
#' })
#' require(visualizeR)
#' temporalPlot(lat.list, xyplot.custom = list(main = "Photoperiod as a function of latitude",
#'                                             ylab = "Daylight hours",
#'                                             xlab = "Date",
#'                                             key = list(corner = c(1,.5),
#'                                                        lines = list(col = 1:4),
#'                                                        horizontal = TRUE,
#'                                                        text = list(c("37ºN","39ºN","41ºN","43ºN")))))

photoperiodGrid <- function(c4r.obj) {
    message("[", Sys.time(), "] Calculating photoperiod...")
    j <- getRefDates(c4r.obj) %>% as.POSIXlt() %>% format(format = "%j") %>% as.integer()
    coords <- getCoordinates(c4r.obj)
    lats <- expand.grid(coords$y, coords$x)[2:1][,2]
    c4r.obj %<>% redim(member = TRUE)
    nmem <- getShape(c4r.obj, "member")
    aux.list <- lapply(1:nmem, function(x) {
        ref <- subsetGrid(c4r.obj, members = x) %>% redim(member = FALSE) 
        ref$Data <- vapply(lats, FUN = geosphere::daylength, j, FUN.VALUE = numeric(length(j))) %>% mat2Dto3Darray(x = coords$x, y = coords$y)
        return(ref)
    })
    c4r.obj <- suppressWarnings(bindGrid(aux.list, dimension = "member"))
    aux.list <- NULL
    c4r.obj$Variable$varName <- "photoperiod"
    attr(c4r.obj$Variable, "description") <- "daylength"
    attr(c4r.obj$Variable, "longname") <- "photoperiod_length"
    attr(c4r.obj$Variable, "units") <- "hours.day-1"
    attr(c4r.obj, "R_package_desc") <- paste0("drought4R-v", packageVersion("drought4R"))
    attr(c4r.obj, "R_package_URL") <- "https://github.com/SantanderMetGroup/drought4R"
    attr(c4r.obj, "R_package_ref") <- "http://dx.doi.org/10.1016/j.envsoft.2018.09.009"
    return(c4r.obj)
}
