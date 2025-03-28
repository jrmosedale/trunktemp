% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/treetrunkcode.R
\name{runmodeltime}
\alias{runmodeltime}
\title{Runs tree trunk model to return a timeseries of trunk temperatures for a given layer and segment}
\usage{
runmodeltime(
  microclim,
  seg,
  lyr,
  treeradius,
  refl = 0.23,
  em = 0.97,
  surfwet = 1,
  kwood,
  cwood,
  rhowood,
  nseg = 16,
  nlyr = 8
)
}
\arguments{
\item{microclim}{- micropoint model outputs of beneath canopy conditions (may be spline interpolated using `splinmicroclimatevars`)}

\item{seg}{- integer of trunk segment for which outputs required}

\item{lyr}{- integer of trunk segment for which outputs required}

\item{treeradius}{- radius of trunk at height to be modelled in metres}

\item{refl}{- trunk reflectance in ??}

\item{em}{- trunk surface emmisivity in ??}

\item{surfwet}{- surface wetness of trunk (0:1)}

\item{kwood}{- array of wood thermal conductivities for each layer}

\item{cwood}{- array of wood specific heat capacities for each layer in J/Kg/K}

\item{rhowood}{- array of wood densities for each layer in kg/m^3}

\item{nseg}{- integer number of model trunk segments}

\item{nlyr}{- integer number of model trunk layers}
}
\value{
array of trunk temperatures (degree C) corresponding to location defined by `seg` and `lyr` for entire time period defined by microclim inputs
}
\description{
Runs tree trunk model to return a timeseries of trunk temperatures for a given layer and segment
}
\examples{
Outputs south-facing, outer layer trunk temperatures
tts<-runmodeltime(microclim, 9, 1, treeradius = 0.75, refl = 0.23, em = 0.97, surfwet = 1, kwood = rep(0.42,8), cwood = rep(3000,8), rhowood = rep(900,8), nseg = 16, nlyr = 8)
plot(tts,type="l",ylim=c(5,20),col="red")
par(new=T)
plot(microclim$tair,type="l",ylim=c(5,20),col="gray")
}
