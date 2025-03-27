\name{M_stomatocyte_L12}
\alias{M_stomatocyte_L12}
\docType{data}
\title{
Stomatocyte shape data for spectral order l<=12.
}
\description{
The MemRBC object keeps the result of CNM for M.C0=-6.5.
}
\usage{data("M_stomatocyte_L12")}
\format{
The format is that of a list, class is MemRBC.
MemRBC-objects have own print, plot and update functions.

}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
$A contains current coefficients of shape

$grd : spatial grid data for plots and integration

$bas : basis function values precomputed for all spectral orders and spatial points

$LA : list of coefficients A preiviously recorded

$history : commands that created M4

and many other objects from previous appplications
}
\source{
 Frickenhaus, S. (2025). MembraneRBC - a numerical modeling laboratory for the stomatocyte-discocyte-echinocyte-transformation of Red Blood Cell shape, ZENODO, DOI: https://doi.org/10.5281/zenodo.13908340
}
\references{
 Frickenhaus, S. (2025). MembraneRBC - a numerical modeling laboratory for the stomatocyte-discocyte-echinocyte-transformation of Red Blood Cell shape, ZENODO, DOI: https://doi.org/10.5281/zenodo.13908340
}
\examples{
data(M_stomatocyte_L12)
plot(M_stomatocyte_L12)
M_stomatocyte_L12$history
}
\keyword{datasets}
