##     effectiveTemp.R Calculation of the effective temperature for improving the Thornthwaite's ET0 estimation method
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

#' @title Effective temperature
#' @description Calculation of the \dQuote{effective temperature} for improving the Thornthwaite's ET0 estimation method
#' @param tasmax Grid of maximum monthly temperature (degC)
#' @param tasmin Grid of minimum monthly temperature (degC)
#' @param k Calibration coefficient. Default to 0.69, as proposed by Pereira and Pruitt (2004)
#' @references \itemize{
#' \item de Camargo, A.P., Marin, F.R., Sentelhas, P. C., Picini, A.G., 1999. Adjust of the Thornthwaite’s method to estimate the potential evapotranspiration for arid and superhumid climates, based on daily temperature amplitude. Revista Brasileira de Agrometeorologia 7.
#' \item Pereira, A.R., Pruitt, W.O., 2004. Adaptation of the Thornthwaite scheme for estimating daily reference evapotranspiration. Agricultural Water Management 66, 251–257. https://doi.org/10.1016/j.agwat.2003.11.003
#' }
#' @details The function is internally used by \code{\link{petGrid}} when the Thornthwaite's method with calibration coefficient k is chosen.
#' @seealso \code{\link{petGrid}}
#' @author J. Bedia
#' @export

effectiveTempGrid <- function(tasmin = NULL, tasmax = NULL, k = 0.69) {
    if (is.null(tasmin) | is.null(tasmax)) {
        stop("Both \'tasmin\' and \'tasmax\' arguments are required for effective temperature computation", call. = FALSE)
    }
    if (is.null(k)) stop("A calibration factor \'k\' is required", call. = FALSE)
    checkDim(tasmin, tasmax, dimensions = c("time", "lat", "lon"))
    checkTemporalConsistency(tasmin, tasmax)
    if ((getTimeResolution(tasmax) != "DD" | getTimeResolution(tasmin) != "DD")) stop("Daily data are required for effective temperature computation", call. = FALSE)
    tasmin <- gridArithmetics(tasmax, 3, tasmin, operator = c("*", "-")) %>% gridArithmetics(., 0.5 * k, operator = "*") 
    tasmax <- NULL
    tasmin$Variable$varName <- "Effective_temperature"
    attr(tasmin$Variable, "description") <- "effective temperature for Thorthwaite's method calibration"
    attr(tasmin$Variable, "longname") <- "effective_temperature"
    attr(tasmin, "R_package_desc") <- paste0("drought4R-v", packageVersion("drought4R"))
    attr(tasmin, "R_package_URL") <- "https://github.com/SantanderMetGroup/drought4R"
    attr(tasmin, "R_package_ref") <- "http://dx.doi.org/10.1016/j.envsoft.2018.09.009"
    return(tasmin)
}
