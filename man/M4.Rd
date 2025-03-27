\name{M4}
\alias{M4}
\docType{data}
\title{
Start data set for a discocyte
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("M4")}
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
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
library(MemRBC)
data("M4")
SetParams(M4) # set global model parameters for re-use
Energy(M4)
plot(M4)
M4$history
attributes(M4$A) # more details, like C0 value,
                 # are stored as attributes to A

sapply(M4$LA,attr,"E") # list energy from recorded list of A

names(M4) # list all objects in M4

C<-updateX_only(M$A,M$grd,M$bas) # make 3d coordinates
plot3d(C$X,aspect=FALSE) # plot spatial points

PlotA(M4)
}
\keyword{datasets}
