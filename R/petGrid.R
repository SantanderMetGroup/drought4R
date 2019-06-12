##     petGrid.R Compute Potential Evapotranspiration Grids
##
##     Copyright (C) 2018 Santander Meteorology Group (http://www.meteo.unican.es)
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
#' \code{"hargreaves"} methods are available (monthly), using the implementation of package \pkg{SPEI}.
#' In addition, \code{"hargreaves-samani"} is available for daily data. See details.
#' @param k.th Optional calibration coefficient for the Thornthwaite method. Unused by default. See Details.
#' @param what Optional character string, only applied for the Hargreaves-Samani method.
#' If set to \code{what = "rad"}, it returns the estimated radiation (it is a function of latitude and date).
#' Otherwise, by default, returns the estimated daily potential evapotranspiration.
#' @param ... Further arguments passed to the PET internals
#' @details
#' This function is a wrapper of the functions with the same name as given in \code{method}
#' from the \pkg{SPEI} package (Begueria and Vicente-Serrano). Monthly input data are thus required. 
#' The latitude of the sites/grid points is internally used.
#' In case of multimember grids (e.g. seasonal forecast data), the PET is calculated for each member sepparately. 
#' 
#' \strong{Calibration coefficient for Thornthwaite}
#' 
#' The use of a calibration coefficient (\code{k.th} argument) can provide better PET estimates under certain conditions. For instance,
#' Camargo \emph{et al.} (1999) found that a value of k=0.72 is the best for estimating monthly ET0, while Pereira and Pruitt (2004)
#' recommended k=0.69 for daily ET0. Trajkovic et al. (2019) propose an optimal calibration factor of 0.65 after an intercomparison of
#' several PET estimation methods and calibration coefficients in northern Serbia. 
#' 
#' \strong{Note:} the calibration factor for the Thorthwaite method requires minimum and maximum temperatures instead of mean temperature, 
#' as it also includes thermal amplitude in its formulation.
#' 
#' @importFrom transformeR array3Dto2Dmat mat2Dto3Darray checkDim getCoordinates getDim
#' @importFrom utils packageVersion
#' @importFrom abind abind
#' @importFrom magrittr %>% 
#' @author J Bedia
#' @export
#' @references 
#' \itemize{
#' \item Camargo AP, Marin FR, Sentelhas PC, Picini AG (1999) Adjust of the Thornthwaite’s method to estimate the potential evapotranspiration for 
#' arid and superhumid climates, based on daily temperature amplitude. Rev Bras Agrometeorol 7(2):251–257
#' \item Pereira, A.R., Pruitt, W.O., 2004. Adaptation of the Thornthwaite scheme for estimating daily reference evapotranspiration. Agricultural Water Management 66, 251–257. https://doi.org/10.1016/j.agwat.2003.11.003
#' \item Trajkovic, S., Gocic, M., Pongracz, R., Bartholy, J., 2019. Adjustment of Thornthwaite equation for estimating evapotranspiration in Vojvodina. Theor Appl Climatol. https://doi.org/10.1007/s00704-019-02873-1
#' }
#' @examples \donttest{
#' # Thorthwaite requires monthly mean temperature data as input:
#' data("tas.cru.iberia")
#' thpet <- petGrid(tas = tas.cru.iberia, method = "thornthwaite")
#' require(transformeR)
#' require(magrittr)
#' require(visualizeR)
#' # This is the climatology (by seasons)
#' seasons <- list("DJF" = c(12,1,2), "MAM" = 3:5, "JJA" = 6:8, "SON" = 9:11)
#' season.list <- lapply(1:length(seasons), function(x) {
#'     subsetGrid(thpet, season = seasons[[x]]) %>% 
#'                aggregateGrid(aggr.y = list(FUN = "sum")) %>% climatology()
#' })
#' mg <- makeMultiGrid(season.list, skip.temporal.check = TRUE)
#' spatialPlot(mg, names.attr = names(seasons), at = seq(0, 525, 10),
#'             backdrop.theme = "coastline",
#'             main = "Mean Potential Evapotranspiration Thornthwaite (1981-2010, mm/season)",
#'             rev.colors = TRUE)
#' # A typical operation is computing trends
#' # Here we are interested in the JJA PET trend in Iberia:
#' # We consider the annually aggregated PET series (n = 30 years)
#' jja.pet <- subsetGrid(thpet, season = 6:8) %>% aggregateGrid(aggr.y = list(FUN = "sum"))
#' trend.thpet <- climatology(jja.pet, clim.fun = list(FUN = "trend.1D",
#'                                                     dates = getRefDates(jja.pet),
#'                                                     method = "kendall"))
#' spatialPlot(trend.thpet, backdrop.theme = "countries",
#'             main = "Mann-Kendall Tau PET trend (JJA, 1981-2010)",
#'             rev.colors = TRUE)
#' pval.estimate <- climatology(jja.pet, clim.fun = list(FUN = "trend.1D",
#'                                                       dates = getRefDates(jja.pet),
#'                                                       method = "kendall",
#'                                                       return.pvalue = TRUE))
#' sig.points <- map.stippling(clim = pval.estimate, threshold = 0.05, condition = "LT", 
#'                             pch = 19, cex = .5, col = "purple")
#' spatialPlot(trend.thpet, backdrop.theme = "countries",
#'             rev.colors = TRUE,
#'             main = "Mann-Kendall Tau PET trend (JJA, 1981-2010)",
#'             sp.layout = list(sig.points))
#' # See further examples in 'help(trend.1D)'
#' }


petGrid <- function(tasmin = NULL,
                    tasmax = NULL,
                    tas = NULL,
                    pr = NULL,
                    method = c("thornthwaite", "hargreaves", "hargreaves-samani"),
                    what = c("PET", "rad"),
                    k.th = NULL,
                    ...) {
    method <- match.arg(method, choices = c("thornthwaite", "hargreaves", "hargreaves-samani"))
    message("[", Sys.time(), "] Computing PET-", method, " ...")
    out <- switch(method,
           "thornthwaite" = petGrid.th(tas, tasmin, tasmax, k.th, ...),
           "hargreaves" = petGrid.har(tasmin, tasmax, pr, ...),
           "hargreaves-samani" = petGrid.hs(tasmin, tasmax, what))
    message("[", Sys.time(), "] Done")
    ## Recover the grid structure -----------------------
    coords <- getCoordinates(out$ref.grid)
    pet.grid <- out$ref.grid
    pet.grid$Data <- lapply(out$et0.list, "mat2Dto3Darray", x = coords$x, y = coords$y) %>% abind(along = -1) %>% unname()
    attr(pet.grid$Data, "dimensions") <- getDim(out$ref.grid)
    pet.grid$Variable$varName <- paste("PET", method, sep = "_")
    attr(pet.grid$Variable, "longname") <- paste("potential_evapotranspiration", method, sep = "_")
    attr(pet.grid$Variable, "description") <- attr(pet.grid$Variable, "longname")
    attr(pet.grid$Variable, "daily_agg_cellfun") <- "sum"
    tres <- getTimeResolution(out$ref.grid)
    attr(pet.grid$Variable, "verification_timestep") <- tres
    if (tres == "MM") {
        attr(pet.grid$Variable, "units") <- "mm.month-1"    
        attr(pet.grid$Variable, "monthly_agg_cellfun") <- "sum"
    } else if (tres == "DD") {
        attr(pet.grid$Variable, "units") <- "mm.day-1"    
    }
    if (method == "thornthwaite" | method == "hargreaves") {
        attr(pet.grid, "origin") <- paste0("Calculated with R package 'SPEI' v",
                                           packageVersion("SPEI"), " using 'drought4R' v",
                                           packageVersion("drought4R"))
    } else {
        attr(pet.grid, "origin") <- paste0("Calculated with 'drought4R' v",
                                           packageVersion("drought4R"))
    }
    attr(pet.grid, "URL") <- "https://github.com/SantanderMetGroup/drought4R"
    pet.grid <- redim(pet.grid, drop = TRUE)
    invisible(pet.grid)
}


#' @importFrom SPEI thornthwaite
#' @importFrom transformeR getCoordinates getSeason array3Dto2Dmat getTimeResolution redim getShape subsetGrid gridArithmetics
#' @keywords internal    
#' @author J Bedia

petGrid.th <- function(tas, tasmin, tasmax, k.th, ...) {
    if (is.null(tas)) {
        if (is.null(tasmin) | is.null(tasmax)) {
            stop("\'tas\' has been set to NULL. Therefore, both \'tasmin\' and \'tasmax\' arguments are required by Thornthwaite method", call. = FALSE)
        }
        tas <- if (!is.null(k.th)) {
            ## Effective temp. See equation 5 Trajkovic et al. 2019
            gridArithmetics(tasmax, 3, tasmin, operator = c("*", "-")) %>% gridArithmetics(., 0.5 * k.th, operator = "*") 
        } else {
            gridArithmetics(tasmax, tasmin, 2, operator = c("+", "/"))
        }
    } else {
        if (!is.null(tasmin) | !is.null(tasmax)) {
            message("\'tas\' argument has been provided. Therefore, \'tasmin\' and \'tasmax\' arguments will be ignored.")
        }
    }
    if (!identical(1:12, getSeason(tas))) stop("The input grid must encompass a whole-year season")
    if (getTimeResolution(tas) != "MM") stop("A monthly input grid is required by the Thornthwaite method", call. = FALSE)
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
#' @importFrom transformeR getCoordinates getSeason array3Dto2Dmat getTimeResolution checkTemporalConsistency
#' @importFrom magrittr %<>% 
#' @keywords internal    
#' @author J Bedia

petGrid.har <- function(tasmin, tasmax, pr, ...) {
    if (is.null(tasmin) || is.null(tasmax) || is.null(pr)) {
        stop("tasmin, tasmax and pr grids are required by Hargreaves method", call. = FALSE)
    }
    if (getTimeResolution(tasmin) != "MM") stop("A monthly input grid is required by the Hargreaves method")
    tasmin %<>% redim(member = TRUE)
    tasmax %<>% redim(member = TRUE)
    pr %<>% redim(member = TRUE)
    suppressMessages(checkDim(tasmin, tasmax, pr))
    checkTemporalConsistency(tasmin, tasmax, pr)
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



#' @import transformeR
#' @author J. Bedia, M. Iturbide
#' @import transformeR
#' @importFrom magrittr extract extract2 %>% %<>% 
#' @keywords internal    
#' @note This function is intended for daily data only. 
#' For hourly or shorter periods, see FAO, eq. 28 p47 (not implemented yet).
#' Radiation will be returned in MJ.m-2.day-1. 
#' @references Allen, R.G., Pereira, L.S., Raes, D., Smith, M., 2006. FAO Irrigation and Drainage Paper. Crop Evapotranspiration (guidelines for computing crop water requirements) (No. 56). FAO.


petGrid.hs <- function(tasmin, tasmax, what) {
    tasmin %<>% redim(member = TRUE)
    tasmax %<>% redim(member = TRUE)
    suppressMessages(checkDim(tasmin, tasmax))
    checkTemporalConsistency(tasmin, tasmax)
    if (isFALSE(getTimeResolution(tasmin) == "DD")) {
        stop("Hargreaves-Samani method is meant for daily data", call. = FALSE)
    }
    if (typeofGrid(tasmin) == "rotated_grid") {
        stop("Rotated grids are unsupported. Please consider regridding", call. = FALSE)
    }
    what = match.arg(what, choices = c("PET", "rad"))
    lats <- getCoordinates(tasmin) 
    if (is.data.frame(lats)) {
        lats %<>% extract2("y")
    } else {
        lats %<>% expand.grid() %>% extract(,2)
    }
    lats <- lats * (pi / 180) ## Conversion decimal degrees --> radians
    refdates <- getRefDates(tasmin)
    J <- refdates %>% as.Date() %>% format("%j") %>% as.integer()
    ind <- getYearsAsINDEX(tasmin) %>% which.leap()
    year.len <- rep(365, length(refdates))
    year.len[ind] <- 366 ## Number of days per year (controlling for leap years)
    ds <- 0.409 * sin(2 * pi * J / year.len - 1.39) ## Solar declination, FAO, eq. 24 (p. 46)
    n.mem <- getShape(tasmin, "member")
    aux.lat <- matrix(lats, nrow = length(J), ncol = length(lats), byrow = TRUE) 
    ws <- acos(-tan(aux.lat) * tan(ds))   ## Sunset hour angle, FAO, eq. 25 (p. 46)
    dr <- 1 + 0.033 * cos(2 * pi * J / year.len) ## inverse relative distance Earth-Sun , FAO, eq. 23 (p. 46)
    Gsc <- 0.082 ## Solar constant (MJ.m-2.min-1)
    Ra <- (24 * 60 * Gsc * dr * (ws * sin(aux.lat) * sin(ds) + cos(aux.lat) * cos(ds) * sin(ws))) / pi ## FAO, eq. 21 (p. 46)
    et0.list <- lapply(1:n.mem, function(x) {
        if (what == "rad") {
            return(Ra)
        } else {
            tn <- subsetGrid(tasmin, members = x, drop = TRUE) %>% extract2("Data") %>% array3Dto2Dmat()
            tx <- subsetGrid(tasmax, members = x, drop = TRUE) %>% extract2("Data") %>% array3Dto2Dmat()
            trng <- tx - tn
            # Ensure that no negative temperature daily ranges exist
            if (any(trng < 0)) {
                trng[which(trng < 0)] <- 0
            }
            .0023 * (Ra / 2.45) * ((tx + tn) / 2 + 17.8) * sqrt(trng) %>% return() ## FAO, eq. 52 (p. 64), applying conversion of units to Ra (MJ.m-2.day-1 --> mm.day-1)
        }
    })
    return(list("et0.list" = et0.list, "ref.grid" = tasmin))
}




# #' @import transformeR
# #' @author J. Bedia, M. Iturbide
# #' @import transformeR
# #' @importFrom magrittr extract extract2 %>% %<>% 
# #' @keywords internal    
# #' Allen, R.G., Pereira, L.S., Raes, D., Smith, M., 2006. FAO Irrigation and Drainage Paper. Crop Evapotranspiration (guidelines for computing crop water requirements) (No. 56). FAO.
# 
# # load("~/workspace/gitRepo/fireDanger/ignore/fffdi_vignette/fffi_data.Rdata", verbose = TRUE)
# # tasmin <- tasmin.fin 
# # tasmax <- tasmax.fin 
# 
# petGrid.pm <- function(tasmin, tasmax, elev, u2, HRm, what) {
#     tasmin %<>% redim(member = TRUE)
#     tasmax %<>% redim(member = TRUE)
#     suppressMessages(checkDim(tasmin, tasmax))
#     checkTemporalConsistency(tasmin, tasmax)
#     if (isFALSE(getTimeResolution(tasmin) == "DD")) {
#         stop("Penman-Monteith method is meant for daily data", call. = FALSE)
#     }
#     if (typeofGrid(tasmin) == "rotated_grid") {
#         stop("Rotated grids are unsupported. Please consider regridding", call. = FALSE)
#     }
#     what = match.arg(what, choices = c("PET", "rad"))
#     lats <- getCoordinates(tasmin) 
#     if (is.data.frame(lats)) {
#         lats %<>% extract2("y")
#     } else {
#         lats %<>% expand.grid() %>% extract(,2)
#     }
#     lats <- lats * (pi / 180) # Conversion to radians
#     refdates <- getRefDates(tasmin)
#     J <- refdates %>% as.Date() %>% format("%j") %>% as.integer()
#     nmonthdays <- sapply(refdates, FUN = "ndays")
#     ind <- getYearsAsINDEX(tasmin) %>% which.leap()
#     year.len <- rep(365, length(refdates))
#     year.len[ind] <- 366 ## Number of days per year (controlling for leap years)
#     Gsc <- 0.082
#     ds <- 0.409 * sin(2 * pi * J / year.len - 1.39) ## Solar declination FAO, eq. 24 (p. 46)
#     n.mem <- getShape(tasmin, "member")
#     aux.lat <- matrix(lats, nrow = length(J), ncol = length(lats), byrow = TRUE) 
#     ws <- acos(-tan(aux.lat) * tan(ds)) ## Sunset hour angle, FAO eq. 25 (p. 46)  
#     N <- (24 / pi) * ws ## Daylight hours (maximum possible duration of sunshine hours), FAO, eq. 34 (p. 48)
#     dr <- 1 + 0.033 * cos(2 * pi * J / year.len)
#     Ra <- (24 * 60 * Gsc * dr * (ws * sin(aux.lat) * sin(ds) + cos(aux.lat) * cos(ds) * sin(ws))) / pi ## FAO, eq. 21 (p. 46)
#     et0.list <- lapply(1:n.mem, function(x) {
#         if (what == "rad") {
#            return(Ra)
#         } else {
#             P <- 101.3 * ((293 - 0.0065 * elev) / 293) ^ 5.26
#             Rs <- (0.25 + 0.5 * (n/N)) * Ra ## Solar radiation using Angstrom formula, FAO, eq. 35 (p. 50)
#             Rso <- (0.75 + 2/100000) * Ra
#             Rns <- 0.77 * Rs  
#             o <- 4.903/1000000000
#             tn <- subsetGrid(tasmin, members = x, drop = TRUE) %>% extract2("Data") %>% array3Dto2Dmat()
#             tx <- subsetGrid(tasmax, members = x, drop = TRUE) %>% extract2("Data") %>% array3Dto2Dmat()
#             tm <- (tx + tn) / 2
#             lam <- 2.501 - (2.361 / 1000) * tm
#             y <- 0.00163 * (P / lam)
#             es <- (0.6108 * exp(17.27 * tn / (tn + 237.3)) + 0.6108 * exp(17.27 * tx / (tx + 237.3))) / 2
#             ea <- es * HRm / 100
#             Rnl <- o * (tm + 273.16) * (0.34 - 0.14 * ea^0.5) * (1.35 * Rs / Rso - 0.35)
#             Rn <- Rns - Rnl
#             A <- (4098 * es) / (tm + 237.3)^2
#             ((0.408 * A * Rn) + (y * 900 * u2 * (es - ea)) / (tm + 273)) / (A + y * (1 + 0.34 * u2)) %>% return()
#         }
#     })
#     return(list("et0.list" = et0.list, "ref.grid" = tasmin))
# }

#' Calculate the number of days of the current month
#' @param d A date (character) in format YYYY-MM-DD...
#' @return The number of days of the current month
#' @references 
#' \url{http://stackoverflow.com/questions/6243088/find-out-the-number-of-days-of-a-month-in-r}
#' @keywords internal
#' @note This function is an internal helper of loadeR, from wehere it has been duplicated (bad practice...)
#' @importFrom utils tail

ndays <- function(d) {
    as.difftime(tail((28:31)[which(!is.na(as.Date(paste0(substr(d, 1, 8), 28:31), '%Y-%m-%d')))], 1), units = "days")
}


