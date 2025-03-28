% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/treetrunkcode.R
\name{runmodel}
\alias{runmodel}
\title{Runs tree trunk model and returns temperatures at time interval n for all trunk segments and layers}
\usage{
runmodel(
  microclim,
  n,
  treeradius,
  refl = 0.23,
  em = 0.97,
  surfwet = 1,
  kwood,
  cwood,
  rhowood,
  nseg = 16,
  nlyr = 8,
  timestep = 60
)
}
\arguments{
\item{microclim}{- micropoint model outputs of beneath canopy conditions (may be spline interpolated using `splinmicroclimatevars`)}

\item{n}{- time interval (seconds from start) for which outputs returned}

\item{treeradius}{- radius of trunk at height to be modelled in metres}

\item{refl}{- trunk reflectance in ??}

\item{em}{- trunk surface emmisivity in ??}

\item{surfwet}{- surface wetness of trunk (0:1)}

\item{kwood}{- array of wood thermal conductivities for each layer}

\item{cwood}{- array of wood specific heat capacities for each layer in J/Kg/K}

\item{rhowood}{- array of wood densities for each layer in kg/m^3}

\item{nseg}{- integer number of model trunk segments}

\item{nlyr}{- integer number of model trunk layers}

\item{timestep}{- input and model time step in seconds}
}
\value{
2D matrix (nseg x nlyr) of trunk temperatures at time interval ii
}
\description{
Runs tree trunk model and returns temperatures at time interval n for all trunk segments and layers
}
\examples{
ii <- which.max(microclim$Rdirdown)
treetemps<-runmodel(microclim, ii, treeradius = 0.75, refl = 0.23, em = 0.97, surfwet = 1, kwood = rep(0.42,8), cwood = rep(3000,8), rhowood = rep(900,8), nseg = 16, nlyr = 8, timestep = 60)
}
