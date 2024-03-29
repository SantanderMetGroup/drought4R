% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rdas.R
\docType{data}
\name{tasmax.eobs.iberia.daily}
\alias{tasmax.eobs.iberia.daily}
\title{Daily maximum temperature data in the Iberian Peninsula for the year 2000}
\format{
A \pkg{climate4R} grid
}
\source{
E-OBS dataset v17 (Haylock \emph{et al.} 2008)
The R dataset has been created with the climate4R package \pkg{loadeR} (Iturbide \emph{et al.} 2019), reading directly from
 the OPeNDAP service at <http://opendap.knmi.nl/knmi/thredds/dodsC/e-obs_0.25regular/tx_0.25deg_reg_v17.0.nc>.
}
\description{
Gridded observations of daily maximum surface temperature for the Iberian Peninsula (2000, in degC)
}
\examples{
data("tasmax.eobs.iberia.daily.rda")
require(visualizeR)
spatialPlot(climatology(tasmin.cru.iberia), backdrop.theme = "countries", rev.colors = TRUE)
}
\references{
\itemize{
\item Haylock, M.R., Hofstra, N., Klein Tank, A.M.G., Klok, E.J., Jones, P.D., New, M., 2008. A European daily high-resolution gridded data set of surface temperature and precipitation for 1950–2006. Journal of Geophysical Research 113. https://doi.org/10.1029/2008JD010201
\item Iturbide, M., Bedia, J., Herrera, S., Baño-Medina, J., Fernández, J., Frías, M.D., Manzanas, R., San-Martín, D., Cimadevilla, E., Cofiño, A.S., Gutiérrez, J.M., 2019. The R-based climate4R open framework for reproducible climate data access and post-processing. Environmental Modelling & Software 111, 42–54. https://doi.org/10.1016/j.envsoft.2018.09.009
}
}
