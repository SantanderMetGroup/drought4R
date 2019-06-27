#' @title Mean monthly temperature data in the Iberian Peninsula 1981-2010
#' @description Gridded observations of monthly mean surface temperature for the Iberian Peninsula (1981-2010, in degC)
#' @name tas.cru.iberia
#' @docType data
#' @format A \pkg{climate4R} grid
#' @source Climate Research Unit, University of East Anglia (CRU TS v4.00), \url{https://crudata.uea.ac.uk/cru/data/hrg/}.
#' The R dataset has been created with package \pkg{loadeR}, of the \strong{climate4R bundle} (Cofino et al 2017)
#' @references \itemize{
#' \item Harris, I., Jones, P.D., Osborn, T.J., Lister, D.H., 2014. Updated high-resolution grids of monthly climatic observations - the CRU TS3.10 Dataset: UPDATED HIGH-RESOLUTION GRIDS OF MONTHLY CLIMATIC OBSERVATIONS. International Journal of Climatology 34, 623-642. doi:10.1002/joc.3711
#' \item Cofino, A.S., Bedia, J., Iturbide, M., Vega, M., Herrera, S., Fernandez, J., Frias, M.D., Manzanas, R., Gutierrez, J.M., 2017. The ECOMS User Data Gateway: Towards seasonal forecast data provision and reseacrh reproducibility in the era of climate services. Climate Services. doi:10.1002/joc.3711
#' }
#' @examples 
#' data("tas.cru.iberia")
#' require(visualizeR)
#' spatialPlot(climatology(tas.cru.iberia), backdrop.theme = "countries", rev.colors = TRUE)
NULL

#' @title Mean maximum monthly temperature data in the Iberian Peninsula 1981-2010
#' @description Gridded observations of monthly mean maximum surface temperature for the Iberian Peninsula (1981-2010, in degC)
#' @name tasmax.cru.iberia
#' @docType data
#' @format A \pkg{climate4R} grid
#' @source Climate Research Unit, University of East Anglia (CRU TS v4.00), \url{https://crudata.uea.ac.uk/cru/data/hrg/}.
#' The R dataset has been created with package \pkg{loadeR}, of the \strong{climate4R bundle} (Cofino et al 2017)
#' @references \itemize{
#' \item Harris, I., Jones, P.D., Osborn, T.J., Lister, D.H., 2014. Updated high-resolution grids of monthly climatic observations - the CRU TS3.10 Dataset: UPDATED HIGH-RESOLUTION GRIDS OF MONTHLY CLIMATIC OBSERVATIONS. International Journal of Climatology 34, 623-642. doi:10.1002/joc.3711
#' \item Cofino, A.S., Bedia, J., Iturbide, M., Vega, M., Herrera, S., Fernandez, J., Frias, M.D., Manzanas, R., Gutierrez, J.M., 2017. The ECOMS User Data Gateway: Towards seasonal forecast data provision and reseacrh reproducibility in the era of climate services. Climate Services. doi:10.1002/joc.3711
#' }
#' @examples 
#' data("tasmax.cru.iberia")
#' require(visualizeR)
#' spatialPlot(climatology(tasmax.cru.iberia), backdrop.theme = "countries", rev.colors = TRUE)
NULL

#' @title Mean minimum monthly temperature data in the Iberian Peninsula 1981-2010
#' @description Gridded observations of monthly mean minimum surface temperature for the Iberian Peninsula (1981-2010, in degC)
#' @name tasmin.cru.iberia
#' @docType data
#' @format A \pkg{climate4R} grid
#' @source Climate Research Unit, University of East Anglia (CRU TS v4.00), \url{https://crudata.uea.ac.uk/cru/data/hrg/}.
#' The R dataset has been created with package \pkg{loadeR}, of the \strong{climate4R bundle} (Cofino et al 2017)
#' @references \itemize{
#' \item Harris, I., Jones, P.D., Osborn, T.J., Lister, D.H., 2014. Updated high-resolution grids of monthly climatic observations - the CRU TS3.10 Dataset: UPDATED HIGH-RESOLUTION GRIDS OF MONTHLY CLIMATIC OBSERVATIONS. International Journal of Climatology 34, 623-642. doi:10.1002/joc.3711
#' \item Cofino, A.S., Bedia, J., Iturbide, M., Vega, M., Herrera, S., Fernandez, J., Frias, M.D., Manzanas, R., Gutierrez, J.M., 2017. The ECOMS User Data Gateway: Towards seasonal forecast data provision and reseacrh reproducibility in the era of climate services. Climate Services. doi:10.1002/joc.3711
#' }
#' @examples 
#' data("tasmin.cru.iberia")
#' require(visualizeR)
#' spatialPlot(climatology(tasmin.cru.iberia), backdrop.theme = "countries", rev.colors = TRUE)
NULL

#' @title Monthly accumulated precipitation data in the Iberian Peninsula (1981-2010)
#' @description Gridded observations of monthly accumulated precipitation for the Iberian Peninsula (1981-2010, in mm/month)
#' @name pr.cru.iberia
#' @docType data
#' @format A \pkg{climate4R} grid
#' @source Climate Research Unit, University of East Anglia (CRU TS v4.00), \url{https://crudata.uea.ac.uk/cru/data/hrg/}.
#' The R dataset has been created with package \pkg{loadeR}, of the \strong{climate4R bundle} (Cofino et al 2017)
#' @references \itemize{
#' \item Harris, I., Jones, P.D., Osborn, T.J., Lister, D.H., 2014. Updated high-resolution grids of monthly climatic observations - the CRU TS3.10 Dataset: UPDATED HIGH-RESOLUTION GRIDS OF MONTHLY CLIMATIC OBSERVATIONS. International Journal of Climatology 34, 623-642. doi:10.1002/joc.3711
#' \item Cofino, A.S., Bedia, J., Iturbide, M., Vega, M., Herrera, S., Fernandez, J., Frias, M.D., Manzanas, R., Gutierrez, J.M., 2017. The ECOMS User Data Gateway: Towards seasonal forecast data provision and reseacrh reproducibility in the era of climate services. Climate Services. doi:10.1002/joc.3711
#' }
#' @examples 
#' data("pr.cru.iberia")
#' require(visualizeR)
#' spatialPlot(climatology(pr.cru.iberia), backdrop.theme = "countries")
NULL




























#' @title Daily minimum temperature data in the Iberian Peninsula for the year 2000
#' @description Gridded observations of daily minimum surface temperature for the Iberian Peninsula (2000, in degC)
#' @name tasmin.eobs.iberia.daily
#' @docType data
#' @format A \pkg{climate4R} grid
#' @source E-OBS dataset v17 (Haylock \emph{et al.} 2008)
#' The R dataset has been created with the climate4R package \pkg{loadeR} (Iturbide \emph{et al.} 2019), reading directly from
#'  the OPeNDAP service at <http://opendap.knmi.nl/knmi/thredds/dodsC/e-obs_0.25regular/tn_0.25deg_reg_v17.0.nc>.
#' @references \itemize{
#' \item Haylock, M.R., Hofstra, N., Klein Tank, A.M.G., Klok, E.J., Jones, P.D., New, M., 2008. A European daily high-resolution gridded data set of surface temperature and precipitation for 1950–2006. Journal of Geophysical Research 113. https://doi.org/10.1029/2008JD010201
#' \item Iturbide, M., Bedia, J., Herrera, S., Baño-Medina, J., Fernández, J., Frías, M.D., Manzanas, R., San-Martín, D., Cimadevilla, E., Cofiño, A.S., Gutiérrez, J.M., 2019. The R-based climate4R open framework for reproducible climate data access and post-processing. Environmental Modelling & Software 111, 42–54. https://doi.org/10.1016/j.envsoft.2018.09.009
#' }
#' @examples 
#' data("tasmin.eobs.iberia.daily.rda")
#' require(visualizeR)
#' spatialPlot(climatology(tasmin.cru.iberia), backdrop.theme = "countries", rev.colors = TRUE)
NULL

#' @title Daily maximum temperature data in the Iberian Peninsula for the year 2000
#' @description Gridded observations of daily maximum surface temperature for the Iberian Peninsula (2000, in degC)
#' @name tasmax.eobs.iberia.daily
#' @docType data
#' @format A \pkg{climate4R} grid
#' @source E-OBS dataset v17 (Haylock \emph{et al.} 2008)
#' The R dataset has been created with the climate4R package \pkg{loadeR} (Iturbide \emph{et al.} 2019), reading directly from
#'  the OPeNDAP service at <http://opendap.knmi.nl/knmi/thredds/dodsC/e-obs_0.25regular/tx_0.25deg_reg_v17.0.nc>.
#' @references \itemize{
#' \item Haylock, M.R., Hofstra, N., Klein Tank, A.M.G., Klok, E.J., Jones, P.D., New, M., 2008. A European daily high-resolution gridded data set of surface temperature and precipitation for 1950–2006. Journal of Geophysical Research 113. https://doi.org/10.1029/2008JD010201
#' \item Iturbide, M., Bedia, J., Herrera, S., Baño-Medina, J., Fernández, J., Frías, M.D., Manzanas, R., San-Martín, D., Cimadevilla, E., Cofiño, A.S., Gutiérrez, J.M., 2019. The R-based climate4R open framework for reproducible climate data access and post-processing. Environmental Modelling & Software 111, 42–54. https://doi.org/10.1016/j.envsoft.2018.09.009
#' }
#' @examples 
#' data("tasmax.eobs.iberia.daily.rda")
#' require(visualizeR)
#' spatialPlot(climatology(tasmin.cru.iberia), backdrop.theme = "countries", rev.colors = TRUE)
NULL

#' @title Daily accumulated precipitation data in the Iberian Peninsula for the year 2000
#' @description Gridded observations of daily accumulated precipitation for the Iberian Peninsula (2000, in mm)
#' @name pr.eobs.iberia.daily
#' @docType data
#' @format A \pkg{climate4R} grid
#' @source E-OBS dataset v17 (Haylock \emph{et al.} 2008)
#' The R dataset has been created with the climate4R package \pkg{loadeR} (Iturbide \emph{et al.} 2019), reading directly from
#'  the OPeNDAP service at <http://opendap.knmi.nl/knmi/thredds/dodsC/e-obs_0.25regular/rr_0.25deg_reg_v17.0.nc>.
#' @references \itemize{
#' \item Haylock, M.R., Hofstra, N., Klein Tank, A.M.G., Klok, E.J., Jones, P.D., New, M., 2008. A European daily high-resolution gridded data set of surface temperature and precipitation for 1950–2006. Journal of Geophysical Research 113. https://doi.org/10.1029/2008JD010201
#' \item Iturbide, M., Bedia, J., Herrera, S., Baño-Medina, J., Fernández, J., Frías, M.D., Manzanas, R., San-Martín, D., Cimadevilla, E., Cofiño, A.S., Gutiérrez, J.M., 2019. The R-based climate4R open framework for reproducible climate data access and post-processing. Environmental Modelling & Software 111, 42–54. https://doi.org/10.1016/j.envsoft.2018.09.009
#' }
#' @examples 
#' data("pr.eobs.iberia.daily.rda")
#' require(visualizeR)
#' spatialPlot(climatology(tasmin.cru.iberia), backdrop.theme = "countries")
NULL

