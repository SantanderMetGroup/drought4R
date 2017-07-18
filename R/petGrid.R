##     petGrid.R Compute Potential Evapotranspiration Grids
##
##     Copyright (C) 2017 Santander Meteorology Group (http://www.meteo.unican.es)
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

#' @title Calculation of Potential Evapotranspiration
#' @description Implementation of several PET methods for the climate4R bundle
#' @param tas Grid of mean monthly temperature (degC)
#' @param tasmax Grid of maximum monthly temperature (degC)
#' @param tasmin Grid of minimum monthly temperature (degC)
#' @param pr Grid of total monthly precipitation amount (mm/month)
#' @param method Potential evapotranspiration method. Currently \code{"thornthwaite"} and 
#' \code{"hargreaves"} methods are available. See details.
#' @param ... Further arguments passed to the PET internals
#' @details
#' This function is a wrapper of the functions with the same name as given in \code{method}
#' from the \pkg{SPEI} package (Begueria and Vicente-Serrano). Monthly input data are thus required. 
#' The latitude of the sites/grid points is internally used.
#' In case of multimember grids (e.g. seasonal forecast data), the PET is calculated for each member sepparately. 
#' @importFrom transformeR array3Dto2Dmat mat2Dto3Darray checkDim getCoordinates getDim
#' @importFrom utils packageVersion
#' @importFrom abind abind
#' @importFrom magrittr %>% 
#' @author J Bedia
#' @export
#' @examples 
#' # Thorthwaite requires monthly mean temperature data as input:
#' data("tas.cru.iberia")
#' thpet <- petGrid(tas = tas.cru.iberia, method = "thornthwaite")
#' require(transformeR)
#' require(magrittr)
#' # This is the climatology (by seasons)
#' seasons <- list("DJF" = c(12,1,2), "MAM" = 3:5, "JJA" = 6:8, "SON" = 9:11)
#' season.list <- lapply(1:length(seasons), function(x) {
#'     subsetGrid(thpet, season = seasons[[x]]) %>% 
#'                aggregateGrid(aggr.y = list(FUN = "sum")) %>% climatology()
#' })
#' mg <- makeMultiGrid(season.list, skip.temporal.check = TRUE)
#' plotClimatology(mg, names.attr = names(seasons), at = seq(0, 525, 10),
#'                 backdrop.theme = "coastline",
#'                 main = "Mean Potential Evapotranspiration Thornthwaite (1981-2010, mm/season)")
#' # A typical operation is computing trends
#' # Here we are interested in the JJA PET trend in Iberia:
#' # We consider the annually aggregated PET series (n = 30 years)
#' jja.pet <- subsetGrid(thpet, season = 6:8) %>% aggregateGrid(aggr.y = list(FUN = "sum"))
#' trend.thpet <- climatology(jja.pet, clim.fun = list(FUN = "trend.1D",
#'                                                     dates = getRefDates(jja.pet),
#'                                                     method = "kendall"))
#' plotClimatology(trend.thpet, backdrop.theme = "countries",
#'                 main = "Mann-Kendall Tau PET trend (JJA, 1981-2010)")
#' pval.estimate <- climatology(jja.pet, clim.fun = list(FUN = "trend.1D",
#'                                                       dates = getRefDates(jja.pet),
#'                                                       method = "kendall",
#'                                                       return.pvalue = TRUE))
#' sig.points <- map.stippling(clim = pval.estimate, threshold = 0.05, condition = "LT", 
#'                             pch = 19, cex = .5, col = "purple")
#' plotClimatology(trend.thpet, backdrop.theme = "countries",
#'                 main = "Mann-Kendall Tau PET trend (JJA, 1981-2010)",
#'                 sp.layout = list(sig.points))
#' # See further examples in 'help(trend.1D)'


petGrid <- function(tasmin = NULL,
                    tasmax = NULL,
                    tas = NULL,
                    pr = NULL,
                    method = c("thornthwaite", "hargreaves"),
                    ...) {
    method <- match.arg(method, choices = c("thornthwaite", "hargreaves"))
    out <- switch(method,
           "thornthwaite" = petGrid.th(tas),
           "hargreaves" = petGrid.har(tasmin, tasmax, pr))
    ## Recover the grid structure -----------------------
    coords <- getCoordinates(out$ref.grid)
    pet.grid <- out$ref.grid
    pet.grid$Data <- lapply(out$et0.list, "mat2Dto3Darray", x = coords$x, y = coords$y) %>% abind(along = -1) %>% unname()
    attr(pet.grid$Data, "dimensions") <- getDim(out$ref.grid)
    pet.grid$Variable$varName <- paste("PET", method, sep = "_")
    attr(pet.grid$Variable, "longname") <- paste("potential_evapotranspiration", method, sep = "_")
    attr(pet.grid$Variable, "units") <- "mm.month-1"
    attr(pet.grid$Variable, "daily_agg_cellfun") <- "sum"
    attr(pet.grid$Variable, "monthly_agg_cellfun") <- "sum"
    attr(pet.grid$Variable, "time_resolution") <- "MM"
    attr(pet.grid, "origin") <- paste0("Calculated with R package 'SPEI' v",
                                       packageVersion("SPEI"), "using R package 'drought4R' v",
                                       packageVersion("drought4R"))
    attr(pet.grid, "URL") <- "https://github.com/SantanderMetGroup/drought4R"
    pet.grid <- redim(pet.grid, drop = TRUE)
    invisible(pet.grid)
}


#' @importFrom SPEI thornthwaite
#' @importFrom transformeR getCoordinates getSeason array3Dto2Dmat getTimeResolution redim getShape subsetGrid
#' @keywords internal    
#' @author J Bedia

petGrid.th <- function(tas, ...) {
    if (is.null(tas)) {
        stop("Mean temperature grid is required by Thornthwaite method", call. = FALSE)
    }
    if (!identical(1:12, getSeason(tas))) stop("The input grid must encompass a whole-year season")
    if (getTimeResolution(tas) != "MM") stop("A monthly input grid is required by the Thornthwaite method")
    tas <- redim(tas, member = TRUE)
    ref.grid <- tas
    coords <- getCoordinates(tas)
    lat <- expand.grid(coords$y, coords$x)[2:1][ ,2]
    n.mem <- getShape(tas, "member")
    et0.list <- lapply(1:n.mem, function(x) {
        aux <- subsetGrid(tas, members = x, drop = TRUE)
        Tave <- array3Dto2Dmat(aux$Data)
        et0 <- matrix(nrow = nrow(Tave), ncol = ncol(Tave))
        for (i in 1:ncol(et0)) {
            et0[,i] <- tryCatch(expr = thornthwaite(Tave = Tave[,i], lat = lat[i], ...),
                                error = function(er) return(rep(NA, nrow(et0)))
            )
        }
        return(et0)
    })
    return(list("et0.list" = et0.list, "ref.grid" = ref.grid))
}

#' @importFrom SPEI hargreaves
#' @importFrom transformeR getCoordinates getSeason array3Dto2Dmat getTimeResolution
#' @keywords internal    
#' @author J Bedia

petGrid.har <- function(tasmin, tasmax, pr, ...) {
    if (is.null(tasmin) || is.null(tasmax) || is.null(pr)) {
        stop("tasmin, tasmax and pr grids are required by Hargreaves method", call. = FALSE)
    }
    if (getTimeResolution(tasmin) != "MM") stop("A monthly input grid is required by the Hargreaves method")
    checkDim(tasmin, tasmax, pr)
    ref.grid <- tasmin
    coords <- getCoordinates(tasmin)
    lat <- expand.grid(coords$y, coords$x)[2:1][ ,2]
    n.mem <- getShape(tasmin, "member")
    et0.list <- lapply(1:n.mem, function(x) {
        aux.tasmin <- subsetGrid(tasmin, members = x, drop = TRUE)
        aux.tasmax <- subsetGrid(tasmax, members = x, drop = TRUE)
        aux.pr <- subsetGrid(pr, members = x, drop = TRUE)
        Tmin <- array3Dto2Dmat(aux.tasmin$Data)
        Tmax <- array3Dto2Dmat(aux.tasmax$Data)
        Pre <- array3Dto2Dmat(aux.pr$Data)
        et0 <- matrix(nrow = nrow(Tmin), ncol = ncol(Tmin))
        for (i in 1:ncol(et0)) {
            et0[,i] <- tryCatch(expr = hargreaves(Tmin = Tmin[,i],
                                                  Tmax = Tmax[,i],
                                                  lat = lat[i],
                                                  Pre = Pre[,i], ...),
                                error = function(er) return(rep(NA, nrow(et0)))
            )
        }
        return(et0)
    })
    return(list("et0.list" = et0.list, "ref.grid" = ref.grid))
}







 
