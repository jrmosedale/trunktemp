#' Spline interpolates microclimate input and filters by month  - REMOVE - USE preparevars
#'
#' @param mout - micropoint model outputs of beneath canopy conditions
#' @param climdata - climate data used as input to micropoint modelling
#' @param mon - integer array of month number(s) (1-12) to include
#' @param lat - latitude
#' @param lon - longitude
#' @param timestep - in seconds for interpolation
#'
#' @return dataframe of splided variables for input into tree trunk model
#' @export
#'
#' @examples
#' mon=c(5,6)
#' microclim<-splinmicroclimatevars(microptout, climdata, mon = mon, lat = 49.96807, lon = -5.215668, timestep = 60)
splinmicroclimatevars<-function(mout,climdata,mon,lat,lon,timestep) {
  tme <- as.POSIXlt(climdata$obs_time,tz="UTC")
  s<-which((tme$mon+1) %in% mon)
  s<-c(s[1]-1,s,s[length(s)]+1)
  mout2 <- mout[s,]
  climdata2 <- climdata[s,]

  # spline microclimate variables
  n<-(length(s)-1)*3600/timestep+1
  tt <- as.numeric(as.POSIXlt(mout2$obs_time,tz="UTC"))
  obs_time<-as.POSIXlt(spline(tt,n=n)$y,origin="1970-01-01 00:00", tz="UTC")
  tair<-spline(mout2$tair,n=n)$y
  relhum<-spline(mout2$relhum,n=n)$y
  windspeed<-spline(mout2$windspeed,n=n)$y
  Rdirdown<-spline(mout2$Rdirdown,n=n)$y
  Rdifdown<-spline(mout2$Rdifdown,n=n)$y
  Rlwdown<-spline(mout2$Rlwdown,n=n)$y
  Rswup<-spline(mout2$Rswup,n=n)$y
  Rlwup<-spline(mout2$Rlwup,n=n)$y
  # cap at zero
  Rdirdown[Rdirdown<0]<-0
  Rdifdown[Rdifdown<0]<-0
  Rswup[Rswup<0]<-0
  windspeed[windspeed<0]<-0 # NEW
  # spline pressure
  pk<-spline(climdata2$pres,n=n)$y
  # Convert dates
  year<-obs_time$year+1900
  month<-obs_time$mon+1
  day<-obs_time$mday
  hour<-obs_time$hour+obs_time$min/60+obs_time$sec/3600
  # Calculate solar zenith
  sp<-solpositionCppv(lat,lon,year,month,day,hour)
  s2<-which((obs_time$mon+1) %in% mon)
  dfout<-data.frame(year=year[s2],
                    month=month[s2],
                    day=day[s2],
                    hour=hour[s2],
                    tair=tair[s2],
                    relhum=relhum[s2],
                    windspeed=windspeed[s2],
                    Rdirdown=Rdirdown[s2],
                    Rdifdown=Rdifdown[s2],
                    Rlwdown=Rlwdown[s2],
                    Rswup=Rswup[s2],
                    Rlwup=Rlwup[s2],
                    pk=pk[s2],
                    zen=sp$zen[s2],
                    azi=sp$azi[s2])
  return(dfout)
}

#' Calculates physical and thermal parameters of tree trunk
#'
#' @param microclim - micropoint model outputs of beneath canopy conditions
#' @param ii
#' @param treeradius - radius of trunk at height to be modelled in metres
#' @param cwood - array of specific heat capacities for each layer of wood in J/Kg/K
#' @param rhowood - array of wood density for each layer of wood in kg/metre^3
#'
#' @param nseg - integer number of model trunk segments
#' @param nlyr - integer number of model trunk layers
#'
#' @return
#' @export
#'
#' @examples
initializevars <- function(microclim, ii, treeradius, cwood, rhowood, nseg, nlyr) {
  # Calculate node distances between layers
  # ** Calculate node thickness
  thick <- sqrt(c(1:nlyr))
  thick <- (thick/ sum(thick)) * treeradius
  # ** Calcuate node position from outer to centre
  p <- 0.5*thick[1]
  for (i in 2:(nlyr-1)) p[i]<-p[i-1]+0.5*thick[i]+0.5*thick[i+1]
  p[nlyr]<-treeradius
  # ** Calcuate distance between layer nodes
  ldist <- p[2:nlyr]-p[1:(nlyr-1)]
  # Calculate cicumference of each layer
  radinner<-0.75-cumsum(thick)
  radouter<-radinner+thick
  circumference<-2*pi*(radinner+radouter)/2
  # Calculate node distances between segments
  sdist<-circumference/nseg
  # Calculate volume of a seqment within each layer
  volume <- sdist * thick
  # Calculate heat capacity of a segment within each layer
  heatcap <- volume*cwood*rhowood # (J/K)
  ptemps<-matrix(microclim$tair[ii],nrow=nlyr,ncol=nseg)
  return(list(ldist=ldist,sdist=sdist,heatcap=heatcap,ptemps=ptemps))
}

#' Prepares trunk model
#'
#' @param microclim - micropoint model outputs of beneath canopy conditions (may be spline interpolated using `splinmicroclimatevars`)
#' @param reps
#' @param treeradius - radius of trunk at height to be modelled in metres
#' @param refl - trunk reflectance in ??
#' @param em - trunk surface emmisivity in ??
#' @param surfwet - surface wetness of trunk (0:1)
#' @param kwood - array of wood thermal conductivities for each layer
#' @param cwood - array of wood specific heat capacities for each layer in J/Kg/K
#' @param rhowood - array of wood densities for each layer in kg/m^3
#' @param nseg - integer number of model trunk segments
#' @param nlyr - integer number of model trunk layers
#' @param timestep - model timestep in seconds
#'
#' @return
#' @export
#'
#' @examples
burnin<-function(microclim, reps = 5, treeradius, refl = 0.23, em = 0.97, surfwet = 1,
                 kwood, cwood, rhowood, nseg = 16, nlyr = 8,
                 timestep = 60)  {
  # initialize model
  initvars<-initializevars(microclim, 1, treeradius, cwood, rhowood, nseg, nlyr)
  temps <- burnincpp(microclim,reps,refl,em,treeradius,surfwet,initvars$ptemps,kwood,
                     initvars$ldist,initvars$sdist,initvars$heatcap,timestep)
  return(temps)
}
#' Runs tree trunk model and returns temperatures at time interval n for all trunk segments and layers
#'
#' @param microclim - micropoint model outputs of beneath canopy conditions (may be spline interpolated using `splinmicroclimatevars`)
#' @param n - time interval (seconds from start) for which outputs returned
#' @param treeradius - radius of trunk at height to be modelled in metres
#' @param refl - trunk reflectance in ??
#' @param em - trunk surface emmisivity in ??
#' @param surfwet - surface wetness of trunk (0:1)
#' @param kwood - array of wood thermal conductivities for each layer
#' @param cwood - array of wood specific heat capacities for each layer in J/Kg/K
#' @param rhowood - array of wood densities for each layer in kg/m^3
#' @param nseg - integer number of model trunk segments
#' @param nlyr - integer number of model trunk layers
#' @param timestep - input and model time step in seconds
#'
#' @return 2D matrix (nseg x nlyr) of trunk temperatures at time interval ii
#' @export
#'
#' @examples
#' ii <- which.max(microclim$Rdirdown)
#' treetemps<-runmodel(microclim, ii, treeradius = 0.75, refl = 0.23, em = 0.97, surfwet = 1, kwood = rep(0.42,8), cwood = rep(3000,8), rhowood = rep(900,8), nseg = 16, nlyr = 8, timestep = 60)
runmodel<-function(microclim, n, treeradius, refl = 0.23, em = 0.97, surfwet = 1,
                 kwood, cwood, rhowood, nseg = 16, nlyr = 8,
                 timestep = 60)  {
   # initialize model
  initvars<-initializevars(microclim, 1, treeradius, cwood, rhowood, nseg, nlyr)
  # burn in model
  temps<-burnin(microclim,5,treeradius,refl,em,surfwet,kwood,cwood,rhowood,nseg,nlyr,timestep)
  # run model
  temps <- runmodelcpp(microclim,n,refl,em,treeradius,surfwet,temps,kwood,
                     initvars$ldist,initvars$sdist,initvars$heatcap,timestep)
  return(temps)
}

#' Runs tree trunk model to return a timeseries of trunk temperatures for a given layer and segment
#'
#' @param microclim - micropoint model outputs of beneath canopy conditions (may be spline interpolated using `splinmicroclimatevars`)
#' @param seg - integer of trunk segment for which outputs required
#' @param lyr - integer of trunk segment for which outputs required
#' @param treeradius - radius of trunk at height to be modelled in metres
#' @param refl - trunk reflectance in ??
#' @param em - trunk surface emmisivity in ??
#' @param surfwet - surface wetness of trunk (0:1)
#' @param kwood - array of wood thermal conductivities for each layer
#' @param cwood - array of wood specific heat capacities for each layer in J/Kg/K
#' @param rhowood - array of wood densities for each layer in kg/m^3
#' @param nseg - integer number of model trunk segments
#' @param nlyr - integer number of model trunk layers
#'
#' @return array of trunk temperatures (degree C) corresponding to location defined by `seg` and `lyr` for entire time period defined by microclim inputs
#' @export
#'
#' @examples
#' Outputs south-facing, outer layer trunk temperatures
#' tts<-runmodeltime(microclim, 9, 1, treeradius = 0.75, refl = 0.23, em = 0.97, surfwet = 1, kwood = rep(0.42,8), cwood = rep(3000,8), rhowood = rep(900,8), nseg = 16, nlyr = 8)
#' plot(tts,type="l",ylim=c(5,20),col="red")
#' par(new=T)
#' plot(microclim$tair,type="l",ylim=c(5,20),col="gray")
runmodeltime <- function(microclim, seg, lyr, treeradius, refl = 0.23, em = 0.97, surfwet = 1,
                 kwood, cwood, rhowood, nseg = 16, nlyr = 8) {
  timestep<-round((microclim$hour[2]-microclim$hour[1])*3600,0)
  # initialize model
  initvars<-initializevars(microclim, 1, treeradius, cwood, rhowood, nseg, nlyr)
  # burn in model
  temps<-burnin(microclim,5,treeradius,refl,em,surfwet,kwood,cwood,rhowood,nseg,nlyr,timestep)
  # run model
  tss <- runmodelthroughtime(microclim,seg,lyr,refl,em,treeradius,surfwet,temps,kwood,initvars$ldist,
    initvars$sdist, initvars$heatcap, timestep)
  return(tss)
}

# Get POSIXct timeseries from numeric year,month,day and decimal hour
get_datetime<-function(year,month,day,hour){
  as.POSIXct(paste(year,sprintf("%02d",month),sprintf("%02d",day),sep='-'))+(3600*hour)
}


