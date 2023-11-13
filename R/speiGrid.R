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
#' @param params A multi-member grid with the distribution parameter values for computing the spei. 
#' Each member corresponds to a parameter. The time dimension length must be 12 (months of the year). 
#' This grid is generated when the parameter \code{return.coefficients} is set as \code{TRUE}. 
#' @param return.coefficients Logical (Default to FALSE). If TRUE, the function returns the parameter values
#' of the distribution that can be further used for computing the index, thus avoiding parameter fitting
#' in subsequent applications of the function. 
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

speiGrid <- function(pr.grid, et0.grid = NULL, scale = 3, params = NULL, return.coefficients = FALSE, ...) {
  arg.list <- list(...)
  arg.list[["scale"]] <- scale
  if(!is.null(params) & isTRUE(return.coefficients)) {
    message("params provided, return.coefficients igonred.")
    return.coefficients <- FALSE
  }
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
  if(n.mem > 1 & (return.coefficients == TRUE | !is.null(params))) stop("When return.coefficients == TRUE or params is different from NULL,
                                                                        Multi-member grids are not accepted.
                                                                        Consider applying this function for each member separately.")
  
  if (!is.null(params)) {
    params.mat <- lapply(1:getShape(params, "member"), function(p) 
      subsetGrid(params, members = p, drop = TRUE) %>% 
        extract2("Data") %>% array3Dto2Dmat()) %>% c(along = 0L) %>% do.call(what = "abind")
  }
  
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
    if(isFALSE(return.coefficients) & is.null(params)) {
      index <- apply(wbalance, MARGIN = 2, FUN = function(i) {
        arg.list[["data"]] <- ts(data = i, 
                                 start = c(yr.start, mon.start),
                                 end = c(yr.end, mon.end),
                                 frequency = 12) 
        do.call("spei", arg.list) %>% extract2("fitted")
      })
      if (length(naind) > 0) index <- index[-naind, ]
      x <- mat2Dto3Darray(index, x = coords$x, y = coords$y)   
      
    } else if (isTRUE(return.coefficients)) {
      index <- lapply(1:ncol(wbalance), FUN = function(i) {
        arg.list[["data"]] <- ts(data = wbalance[,i], 
                                 start = c(yr.start, mon.start),
                                 end = c(yr.end, mon.end),
                                 frequency = 12) 
        do.call("spei", arg.list) %>% extract2("coefficients")
      }) %>% c(along = 2L) %>% do.call(what = "abind")
      x <- lapply(1:dim(index)[1], function(p) {
        mat.index <- t(index[p,,])
        if (length(naind) > 0) mat.index <- mat.index[-naind, ]
        mat2Dto3Darray(mat.index, x = coords$x, y = coords$y)
      }) %>% c(along = 0) %>% do.call(what = "abind")
    
    } else if (!is.null(params)) {
      index <- lapply(1:ncol(wbalance), FUN = function(i) {
        arg.list[["data"]] <- ts(data = wbalance[,i], 
                                 frequency = 12) 
        arg.list[["params"]] <- params.mat[,,i] %>% abind(along = 1.5)
        do.call("spei", arg.list) %>% extract2("fitted")
      }) 
      if (length(naind) > 0) index <- index[-naind, ]
      x <- mat2Dto3Darray(index, x = coords$x, y = coords$y) 
    } 
  })
  message("[", Sys.time(), "] Done")
  ## Recover the grid structure -----------------------
  pr.grid$Data <- do.call("abind", c(spei.list, along = -1L)) %>% unname()
  if (isTRUE(return.coefficients))   pr.grid <- redim(pr.grid, drop = TRUE)
  attr(pr.grid$Data, "dimensions") <- dimNames
  pr.grid <- redim(pr.grid, drop = TRUE)
  pr.grid$Variable$varName <- if (return.coefficients) paste(method, "params") else method
  longname <- ifelse(grepl("SPI", method), "Standardized Precipitation Index", "Standardized Precipitation-Evapotranspiration Index")
  description <- ifelse(return.coefficients, paste("Distribution parametters for period ", c(ref.start, "-", ref.end) %>% paste(collapse = "")), paste("Index using ref period", c(ref.start, "-", ref.end) %>% paste(collapse = "")))
  attr(pr.grid$Variable, "longname") <- if (return.coefficients) paste("params for", longname, scale) else paste(longname, scale)
  attr(pr.grid$Variable, "description") <- description
  attr(pr.grid$Variable, "units") <- "dimensionless"
  attr(pr.grid$Variable, "daily_agg_cellfun") <- "sum"
  attr(pr.grid$Variable, "monthly_agg_cellfun") <- "sum"
  attr(pr.grid$Variable, "verification_timestep") <- "MM"
  if (isTRUE(return.coefficients)) pr.grid[["Members"]] <- if(grepl("SPI", method)) c("alpha", "beta") else c("xi", "alpha", "kappa")
  attr(pr.grid, "origin") <- paste0("Calculated with R package 'SPEI' v",
                                    packageVersion("SPEI"), " using R package 'drought4R' v",
                                    packageVersion("drought4R"))
  attr(pr.grid, "URL") <- "https://github.com/SantanderMetGroup/drought4R"
  invisible(pr.grid)
}
