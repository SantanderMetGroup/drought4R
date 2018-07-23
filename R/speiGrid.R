##     speiGrid.R Computation of Standardized Precipitation(-Evapotranspiration) Index Grids
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
##     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
##     GNU General Public License for more details.
## 
##     You should have received a copy of the GNU General Public License
##     along with this program. If not, see <http://www.gnu.org/licenses/>.

#' @title Computation of Standardized Precipitation(-Evapotranspiration) Index Grids
#' @description Returns SPI/SPEI grid using precipitation (and PET for SPEI) input grids
#' @param pr.grid Precipitation grid (monthly accumulated values, in mm)
#' @param et0.grid Potential evapotranspiration grid (same units as \code{pr.grid})
#' @param scale Integer. Time scale at which the SPEI/SPI are computed. Default to 3 (months)
#' @param ... Further arguments passed to \code{\link[SPEI]{spei}}
#' @details The function is a wrapper of function \code{\link[SPEI]{spei}} from package \pkg{SPEI} adapted
#' to \strong{climate4R} input grids
#' @return A \strong{climate4R} grid with SPEI/SPI data
#' @importFrom magrittr %>% extract2 %<>% 
#' @importFrom transformeR redim checkDim getSeason getTimeResolution getShape subsetGrid array3Dto2Dmat getRefDates checkTemporalConsistency
#' @importFrom SPEI spei
#' @importFrom abind abind
#' @importFrom stats ts
#' @export
#' @seealso \code{\link{petGrid}}, for PET calculation
#' @examples 
#' # By default, et0.grid is null, and SPI is computed from precipitation: 
#' data("pr.cru.iberia")
#' spi3 <- speiGrid(pr.grid = pr.cru.iberia, scale = 3, na.rm = TRUE)
#' ## If PET is used, then SPEI is calculated
#' data("tas.cru.iberia")
#' et0.grid <- petGrid(tas = tas.cru.iberia, method = "thornthwaite")
#' spei3 <- speiGrid(pr.cru.iberia, et0.grid = et0.grid, scale = 3, na.rm = TRUE)

speiGrid <- function(pr.grid, et0.grid = NULL, scale = 3, ...) {
    pr.grid <- redim(pr.grid, member = TRUE)
    if (!identical(1:12, getSeason(pr.grid))) stop("The input grid must encompass a whole-year season", call. = FALSE)
    if (getTimeResolution(pr.grid) != "MM") stop("A monthly input grid is required for SPI/SPEI computation", call. = FALSE)
    if (!is.null(et0.grid)) {
        et0.grid <- redim(et0.grid, member = TRUE)
        suppressMessages(checkDim(pr.grid, et0.grid))
        checkTemporalConsistency(pr.grid, et0.grid)
    }
    coords <- getCoordinates(pr.grid)
    dimNames <- getDim(pr.grid)
    n.mem <- getShape(pr.grid, "member")
    method <- paste(ifelse(is.null(et0.grid), "SPI", "SPEI"), scale, sep = "-")
    datelim <- getRefDates(pr.grid) %>% range() 
    yr.start <- substr(datelim[1], start = 1, stop = 4) %>% as.integer()  
    mon.start <- substr(datelim[1], start = 6, stop = 7) %>% as.integer()  
    yr.end <- substr(datelim[2], start = 1, stop = 4) %>% as.integer()  
    mon.end <- substr(datelim[2], start = 6, stop = 7) %>% as.integer()  
    grid.dates <- getRefDates(pr.grid) %>% substr(start = 1, stop = 7)
    all.dates <- seq(as.Date(datelim[1]), as.Date(datelim[2]), by = "month") %>% as.character() %>% substr(start = 1, stop = 7)
    naind <- which(!all.dates %in% grid.dates)
    message("[", Sys.time(), "] Computing ", method, " ...")
    spei.list <- lapply(1:n.mem, function(x) {
        pr <- subsetGrid(pr.grid, members = x, drop = TRUE) %>% extract2("Data") %>% array3Dto2Dmat()
        if (length(naind) > 0) {
            aux <- matrix(NA, nrow = length(all.dates), ncol = ncol(pr))
            aux[-naind, ] <- pr
            pr <- aux
            aux <- NULL
        }
        if (is.null(et0.grid)) {
            pet <- matrix(0, nrow = nrow(pr), ncol = ncol(pr))
        } else {
            pet <- subsetGrid(et0.grid, members = x, drop = TRUE) %>% extract2("Data") %>% array3Dto2Dmat() 
            if (length(naind) > 0) {
                aux <- matrix(NA, nrow = length(all.dates), ncol = ncol(pet))
                aux[-naind, ] <- pet
                pet <- aux
                aux <- NULL
            }
        }
        wbalance <- pr - pet
        pt <- pet <- NULL
        #
        index <- apply(wbalance, MARGIN = 2, FUN = function(i) {
            ts(data = i,
               start = c(yr.start, mon.start),
               end = c(yr.end, mon.end),
               frequency = 12) %>% spei(scale = scale, ...) %>% extract2("fitted")
        })
        if (length(naind) > 0) index <- index[-naind, ]
        mat2Dto3Darray(index, x = coords$x, y = coords$y)
    })
    message("[", Sys.time(), "] Done")
    ## Recover the grid structure -----------------------
    pr.grid$Data <- do.call("abind", c(spei.list, along = -1L)) %>% unname()
    attr(pr.grid$Data, "dimensions") <- dimNames
    pr.grid <- redim(pr.grid, drop = TRUE)
    pr.grid$Variable$varName <- method
    longname <- ifelse(grepl("SPI", method), "Standardized Precipitation Index", "Standardized Precipitation-Evapotranspiration Index")
    attr(pr.grid$Variable, "longname") <- paste(longname, scale)
    attr(pr.grid$Variable, "units") <- "dimensionless"
    attr(pr.grid$Variable, "daily_agg_cellfun") <- "sum"
    attr(pr.grid$Variable, "monthly_agg_cellfun") <- "sum"
    attr(pr.grid$Variable, "verification_timestep") <- "MM"
    attr(pr.grid, "origin") <- paste0("Calculated with R package 'SPEI' v",
                                       packageVersion("SPEI"), " using R package 'drought4R' v",
                                       packageVersion("drought4R"))
    attr(pr.grid, "URL") <- "https://github.com/SantanderMetGroup/drought4R"
    invisible(pr.grid)
}
