#' @title Calculation of Potential Evapotranspiration
#' @description Implementation of several PET methods 
#' @param tas Grid of mean monthly temperature (degC)
#' @param tasmax Grid of maximum monthly temperature (degC)
#' @param tasmin Grid of minimum monthly temperature (degC)
#' @param pr Grid of total monthly precipitation amount (mm/month)
#' @param method Potential evapotranspiration method. Currently \code{"thornthwaite"} and 
#' \code{"hargreaves"} methods are available. See details.
#' @details
#' This function is a wrapper of the functions with the same name as given in \code{method}
#' from the \pkg{SPEI} package (Begueria and Vicente-Serrano). Monthly input data are thus required.  
#' @importFrom transformeR array3Dto2Dmat mat2Dto3Darray checkDim getCoordinates getDim
#' @importFrom utils packageVersion
#' @importFrom abind abind
#' @importFrom magrittr %>% 
#' @author J Bedia
#' @export



petGrid <- function(tasmin = NULL,
                    tasmax = NULL,
                    tas = NULL,
                    pr = NULL,
                    method = c("thornthwaite", "hargreaves")) {
    method <- match.arg(method, choices = c("thornthwaite", "hargreaves"))
    out <- switch(method,
           "thornthwaite" = petGrid.th(tas),
           "hargreaves" = petGrid.har(tasmin, tasmax, pr))
    ## Recover the grid structure -----------------------
    coords <- getCoordinates(out$ref.grid)
    pet.grid <- out$ref.grid
    pet.grid$Data <- lapply(out$et0.list, "mat2Dto3Darray", x = coords$x, y = coords$y) %>% abind::abind(along = 1)
    attr(pet.grid$Data, "dimensions") <- getDim(out$ref.grid)
    pet.grid$Variable$varName <- paste("PET", substr(method, 1, 3), sep = "_")
    attr(pet.grid$Variable, "longname") <- paste("potential_evapotranspiration", method, sep = "_")
    attr(pet.grid$Variable, "units") <- "mm.month-1"
    attr(pet.grid$Variable, "daily_agg_cellfun") <- "sum"
    attr(pet.grid$Variable, "monthly_agg_cellfun") <- "sum"
    attr(pet.grid$Variable, "verification_time") <- "MM"
    attr(pet.grid, "origin") <- paste0("Calculated with R package 'SPEI' v",
                                       packageVersion("SPEI"), "using R package 'dRought' v")
    attr(pet.grid, "URL") <- "https://github.com/SantanderMetGroup/dRought"
    pet.grid <- redim(pet.grid, drop = TRUE)
    invisible(pet.grid)
}


#' @importFrom SPEI thornthwaite
#' @importFrom transformeR getCoordinates getSeason array3Dto2Dmat getTimeResolution redim getShape subsetGrid
#' @keywords internal    
#' @author J Bedia

petGrid.th <- function(tas) {
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
            et0[,i] <- tryCatch(expr = thornthwaite(Tave = Tave[,i], lat = lat[i]),
                                error = function(er) return(rep(NA, nrow(et0)))
            )
        }
        return(et0)
    })
    return(list("pet" = et0.list, "ref.grid" = ref.grid))
}

#' @importFrom SPEI hargreaves
#' @importFrom transformeR getCoordinates getSeason array3Dto2Dmat getTimeResolution
#' @keywords internal    
#' @author J Bedia

petGrid.har <- function(tasmin, tasmax, pr) {
    if (is.null(tasmin) || is.null(tasmax) || is.null(pr)) {
        stop("tasmin, tasmax and pr grids are required by Hargreaves method", call. = FALSE)
    }
    if (getTimeResolution(tasmin) != "MM") stop("A monthly input grid is required by the Hargreaves method")
    checkDim(tasmin, tasmax, pr)
    ref.grid <- tasmin
    coords <- getCoordinates(tasmin)
    lat <- expand.grid(coords$y, coords$x)[2:1][ ,2]
    Tmin <- array3Dto2Dmat(tasmin$Data)
    Tmax <- array3Dto2Dmat(tasmax$Data)
    Pre <- array3Dto2Dmat(pr$Data)
    et0 <- matrix(nrow = nrow(Tmin), ncol = ncol(Tmin))
    for (i in 1:ncol(et0)) {
        et0[,i] <- tryCatch(expr = hargreaves(Tmin = Tmin[,i],
                                              Tmax = Tmax[,i],
                                              lat = lat[i],
                                              Pre = Pre[,i]),
                            error = function(er) return(rep(NA, nrow(et0)))
        )
    }
    return(list("pet" = et0, "ref.grid" = ref.grid))
}

 
