#' Prepares microclimate and climate data for treetrunk model with option to interpolate to smaller timestep and filter by month(s)
#' REPLACES splinemicroclimate vars AND .sortdata functions
#' @param mout - micropoint model outputs of beneath canopy conditions
#' @param climdata - climate data used as input to micropoint modelling thst produced mout
#' @param lat - latitude
#' @param lon - longitude
#' @param mon - integer array of month number(s) (1-12) to include. If NA all months included.
#' @param timestep - in seconds for interpolation - if NA then no interpolation takes place.
#'
#' @return dataframe of splided variables for input into tree trunk model
#' @export
#'
#' @examples
#' mon=c(5,6)
#' microclim<-splinmicroclimatevars(microptout, climdata, mon = mon, lat = 49.96807, lon = -5.215668, timestep = 60)
preparevars<-function(mout,climdata,lat,lon,mon=NA,timestep=NA) {
  tme <- as.POSIXlt(climdata$obs_time,tz="UTC")

  # Interpolate microclimate variables to smaller timestep if requested
  if(is.na(timestep)==FALSE){
    # Filter by month if requested
    if(anyNA(mon)==FALSE){
      s<-which((tme$mon+1) %in% mon)
      s<-c(s[1]-1,s,s[length(s)]+1)
      mout <- mout[s,]
      climdata <- climdata[s,]
    }
    # Interpolate
    n<-(length(s)-1)*3600/timestep+1
    tt <- as.numeric(as.POSIXlt(mout$obs_time,tz="UTC"))
    obs_time<-as.POSIXlt(spline(tt,n=n)$y,origin="1970-01-01 00:00", tz="UTC")
    tair<-spline(mout$tair,n=n)$y
    relhum<-spline(mout$relhum,n=n)$y
    windspeed<-spline(mout$windspeed,n=n)$y
    Rdirdown<-spline(mout$Rdirdown,n=n)$y
    Rdifdown<-spline(mout$Rdifdown,n=n)$y
    Rlwdown<-spline(mout$Rlwdown,n=n)$y
    Rswup<-spline(mout$Rswup,n=n)$y
    Rlwup<-spline(mout$Rlwup,n=n)$y
    # cap at zero
    Rdirdown[Rdirdown<0]<-0
    Rdifdown[Rdifdown<0]<-0
    Rswup[Rswup<0]<-0
    windspeed[windspeed<0]<-0 # NEW
    # spline pressure
    pk<-spline(climdata$pres,n=n)$y
  } else{
    obs_time<-tme
    tair<-mout$tair
    relhum<-mout$relhum
    windspeed<-mout$windspeed
    Rdirdown<-mout$Rdirdown
    Rdifdown<-mout$Rdifdown
    Rlwdown<-mout$Rlwdown
    Rswup<-mout$Rswup
    Rlwup<-mout$Rlwup
    pk<-climdata$pres
 }

  # Calculate time & solar variables
  year<-obs_time$year+1900
  month<-obs_time$mon+1
  day<-obs_time$mday
  hour<-obs_time$hour+obs_time$min/60+obs_time$sec/3600

  # Get solar angles
  sp<-solpositionCppv(lat,lon,year,month,day,hour)

  # Prepare output
  if(anyNA(mon)==FALSE) s2<-which((obs_time$mon+1) %in% mon) else s2<-c(1:length(obs_time))
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

# Prepares model inputs - REMOVE - USE preparevars
.sortdata<-function(microclim, climdata, lat, long) {
  tme<-as.POSIXlt(climdata$obs_time, tz = "UTC")
  hour <- tme$hour + tme$min/60 + tme$sec/3600
  sp<-solpositionCppv(lat, long, tme$year + 1900, tme$mon + 1, tme$mday, hour)
  dfout<-data.frame(year=tme$year+1900,month=tme$mon+1,day=tme$mday,hour=hour,
                    tair=microclim$tair,relhum=microclim$relhum,windspeed=microclim$windspeed,
                    Rdirdown=microclim$Rdirdown,Rdifdown=microclim$Rdifdown,
                    Rlwdown=microclim$Rlwdown,Rswup=microclim$Rswup,
                    Rlwup=microclim$Rlwup,pk=climdata$pres,
                    zen=sp$zen,azi=sp$azi)
  return(dfout)
}
#' Calculates physical and thermal parameters of tree trunk
#'
#' @param microclim - micropoint model outputs of beneath canopy conditions
#' @param tradius - radius of trunk at height to be modelled in metres
#' @param cs - array of specific heat capacities for each layer of wood in J/Kg/K
#' @param rho - array of wood density for each layer of wood in kg/metre^3
#' @param depth - depth of interest for modelling (metres)
#'
#' @return list of trunk parameters for modelling temperature across 8 layers and 16 segments,
#' ensuring depth of interest corresponds to a layer node
#' @export
#'
.initializevars <- function(microclim, tradius, cs, rho, depth) {
  # Calculate node distances between layers
  # ** Calculate node thickness
  thick <- sqrt(c(1:8))
  thick <- (thick/ sum(thick)) * tradius
  # ** Calcuate node position from outer to centre
  p <- 0.5*thick[1]
  for (i in 2:7) p[i]<-p[i-1]+0.5*thick[i]+0.5*thick[i+1]
  p[8]<-tradius # NO - tradius corresponds to inner edge not node???????
  #p[8]<-tradius-(0.5*thick[7]) # if node !=radius

  # ** Correct node based on depth
  d <- abs(p - depth)
  s <- which.min(d)
  p[s] <- depth

  # ** correct thickness based on depth
  oedge <- 0
  #for (i in 2:8) oedge[i] <- p[i-1] + (p[i-1] - oedge[i - 1]) # REPLACED by BELOW
  for (i in 2:8) oedge[i] <-p[i]-(p[i] - p[i-1])/2
  iedge <- c(oedge[2:8], tradius)
  thick <- iedge - oedge

  # ** Calcuate distance between layer nodes
  ldist <- p[2:8]-p[1:7]
  # Calculate cicumference of each layer
  radinner<-tradius-iedge
  radouter<-tradius-oedge
  circumference<-2*pi*(radinner+radouter)/2
  # Calculate node distances between segments
  sdist<-circumference/8
  # Calculate volume of a seqment within each layer
  volume <- sdist * thick
  # Calculate heat capacity of a segment within each layer
  heatcap <- volume*cs*rho # (J/K)
  ptemps<-matrix(microclim$tair[1],nrow=8,ncol=16) # correct to 16 cols as 16 segs
  return(list(ldist=ldist,sdist=sdist,heatcap=heatcap,nodes=p,ptemps=ptemps))
}

.initializevars_v2 <- function(microclim, tradius, thick, cs, rho, depth) {
  # Calculate node distances between layers
  # ** Calcuate node position from outer to centre
  oedge <- 0
  for (i in 2:8) oedge[i] <- sum(thick[1:(i-1)])# REPLACED by BELOW
  iedge <- c(oedge[2:8], tradius)

  # Calculate nodes
  p <- 0.5*thick[1]
  for (i in 2:8) p[i]<-sum(thick[1:(i-1)]) + (0.5*thick[i])
  p[8]<-tradius # ????tradius corresponds to inner edge not node???????

  # ** Calcuate distance between layer nodes
  ldist <- p[2:8]-p[1:7]
  # Calculate cicumference of each layer
  radinner<-tradius-iedge
  radouter<-tradius-oedge
  circumference<-2*pi*(radinner+radouter)/2
  # Calculate node distances between segments
  sdist<-circumference/8
  # Calculate volume of a seqment within each layer
  volume <- sdist * thick
  # Calculate heat capacity of a segment within each layer
  heatcap <- volume*cs*rho # (J/K)
  ptemps<-matrix(microclim$tair[1],nrow=8,ncol=16) # correct to 16 cols as 16 segs
  return(list(ldist=ldist,sdist=sdist,heatcap=heatcap,nodes=p,ptemps=ptemps))
}


# QUESTION - how do inputs related to splineing - CORRECTED seg number
#' Title
#'
#' @param microclim - prepared microclimate/weather data output by `preparevars` function
#' @param climdata - REMOVED!!! Use preparevars prior
#' @param depth - from trunk surface of location to model (metres)
#' @param aspect - of trunk location to model in degrees (north = 1, south = 180)
#' @param treeparams - named list of tradius=tree radius, cs=specific heat capacity of each layer, and rho = density of each layer
#' @param lat
#' @param long
#'
#' @return
#' @export
#'
#' @examples
runtreeemulator <- function(treeclimdata, depth, aspect, treeparams) {
  # Prepare input climate data
  #treeclimdata <- .sortdata(microclim, climdata, lat, long) # USE sep call to preparevars to allow for splining etc
  # Prepare input tree data - internal function assumes 16 segs and 8 layers
  tvars <- .initializevars(treeclimdata, treeparams$tradius, treeparams$cs, treeparams$rho, depth)
  # find which layer and segment
  dif <- abs(tvars$nodes-depth)
  lyr <- which.min(dif)
  aspts <- c(0:15)*360/16 # corrected 8 to 16 segments
  dif <- abs(aspect-aspts)
  seg <- which.min(dif)
  # burn in model
  temps <- burnincpp(treeclimdata,5,treeparams$refl,treeparams$em,treeparams$tradius,
                     treeparams$surfwet,tvars$ptemps,treeparams$k,
                     tvars$ldist,tvars$sdist,tvars$heatcap,3600)
  # run model
  tss <- runmodelthroughtime(treeclimdata,seg,lyr,treeparams$refl,treeparams$em,
                             treeparams$tradius,treeparams$surfwet,temps,treeparams$k,
                             tvars$ldist,tvars$sdist,tvars$heatcap,3600)
  return(tss)
}

#' Title Run treetrunk emulator with pre-defined layer widths
#'
#' @param microclim - prepared microclimate/weather data output by `preparevars` function
#' @param climdata - REMOVED!!! Use preparevars prior
#' @param depth - from trunk surface of location to model (metres)
#' @param aspect - of trunk location to model in degrees (north = 1, south = 180)
#' @param widths - vector of 8 layer widths
#' @param treeparams - named list of tradius=tree radius, cs=specific heat capacity of each layer, and rho = density of each layer
#' @param lat
#' @param long
#'
#' @return
#' @export
#'
#' @examples
runtreeemulator_v2 <- function(treeclimdata, depth, aspect, treeparams) {
  # Prepare input climate data
  #treeclimdata <- .sortdata(microclim, climdata, lat, long) # USE sep call to preparevars to allow for splining etc
  # Prepare input tree data - internal function assumes 16 segs and 8 layers
  tvars <- .initializevars_v2(treeclimdata, treeparams$tradius, treeparams$layer_widths,treeparams$cs, treeparams$rho, depth)
  # find which layer and segment correspond to aspect and depth
  dif <- abs(tvars$nodes-depth)
  lyr <- which.min(dif)
  aspts <- c(0:15)*360/16 # corrected 8 to 16 segments
  dif <- abs(aspect-aspts)
  seg <- which.min(dif)
  # burn in model
  temps <- burnincpp(treeclimdata,5,treeparams$refl,treeparams$em,treeparams$tradius,
                     treeparams$surfwet,tvars$ptemps,treeparams$k,
                     tvars$ldist,tvars$sdist,tvars$heatcap,3600)
  # run model
  tss <- runmodelthroughtime(treeclimdata,seg,lyr,treeparams$refl,treeparams$em,
                             treeparams$tradius,treeparams$surfwet,temps,treeparams$k,
                             tvars$ldist,tvars$sdist,tvars$heatcap,3600)
  return(tss)
}


.is <- function(r) {
  if (class(r)[1] == "PackedSpatRaster") r<-rast(r)
  if (class(r)[1] != "matrix") {
    if (dim(r)[3] > 1) {
      y<-as.array(r)
    } else y<-as.matrix(r,wide=TRUE)
  } else y<-r
  y
}
.rast <- function(m,tem) {
  r<-rast(m)
  ext(r)<-ext(tem)
  crs(r)<-crs(tem)
  r
}
#' Latitudes from SpatRaster object
.latsfromr <- function(r) {
  e <- ext(r)
  lts <- rep(seq(e$ymax - res(r)[2] / 2, e$ymin + res(r)[2] / 2, length.out = dim(r)[1]), dim(r)[2])
  lts <- array(lts, dim = dim(r)[1:2])
  lts
}
#' Longitudes from SpatRaster object
.lonsfromr <- function(r) {
  e <- ext(r)
  lns <- rep(seq(e$xmin + res(r)[1] / 2, e$xmax - res(r)[1] / 2, length.out = dim(r)[2]), dim(r)[1])
  lns <- lns[order(lns)]
  lns <- array(lns, dim = dim(r)[1:2])
  lns
}
#' lats and longs from SpatRaster object, including reprojection
.latslonsfromr <- function(r) {
  lats<-.latsfromr(r)
  lons<-.lonsfromr(r)
  xy<-data.frame(x=as.vector(lons),y=as.vector(lats))
  xy <- sf::st_as_sf(xy, coords = c('x', 'y'), crs = crs(r))
  ll <- sf::st_transform(xy, 4326)
  ll <- data.frame(lat = sf::st_coordinates(ll)[,2],
                   long = sf::st_coordinates(ll)[,1])
  lons<-array(ll$long,dim=dim(lons))
  lats<-array(ll$lat,dim=dim(lats))
  return(list(lats=lats,lons=lons))
}



era5_process <- function(filein, rte) {
   # get time - method varies as era5 file format changed
   nc_file <- ncdf4::nc_open(filein)
   era5vars <- names(nc_file$var)
   if ("expver" %in% era5vars)
     tme <- as.POSIXlt(nc_file$var$expver$dim[[1]]$vals, tz = "GMT",
                       origin = "1970-01-01")
   if (!"expver" %in% era5vars) tme <- as.POSIXlt(time(rast(filein, "t2m")), tz = "UTC")
   ncdf4::nc_close(nc_file)

   # extract and reproject variables
   varn <- c("t2m", "d2m", "sp", "u10" , "v10",  "tp", "msdwlwrf", "fdir", "ssrd")
   rlst <- list()
   for (i in 1:9) {
     r<-rast(filein, subds=varn[i])
     rlst[[i]]<-project(r,rte)
   }
   # convert variables
   # ~~ temperature
   rlsto<-list()
   rlsto[[1]]<-rlst[[1]]-273.15
   # ~~ relative humidity
   rh <- .rast((satvapCpp(.is(rlst[[2]])-273.15) /  satvapCpp(.is(rlsto[[1]]))) * 100, rte)
   rh[rh > 100]<-100
   rlsto[[2]]<-rh
   # ~~ pressure
   rlsto[[3]] <- rlst[[3]]/1000
   # ~~ swdown
   rlsto[[4]] <- rlst[[9]]/3600
   # ~~ difrad
   dni <- rlst[[8]]/3600
   ll<-.latslonsfromr(rte)
   si <- .rast(solarindexarray(tme$year+1900, tme$mon+1, tme$mday, tme$hour, ll$lats, ll$lons), rte)
   rlsto[[5]]  <- rlsto[[4]]  - (si * dni)
   # ~~ lwdown
   rlsto[[6]]  <- rlst[[7]]
   # wind speed
   rlsto[[7]]<-sqrt(rlst[[4]]^2+rlst[[5]]^2)*0.7477849 # Wind speed (m/s)
   # wind direction
   rlsto[[8]]<-(atan2(rlst[[4]],rlst[[5]])*180/pi+180)%%360
   # precipitation
   rlsto[[9]]<- rlst[[6]] * 1000
   names(rlsto)<-c("temp","relhum","pres","swdown","difrad","lwdown","windspeed","winddir","precip")
   for (i in 1:9)  terra::time(rlsto[[i]])<-as.POSIXct(tme)
   return(rlsto)
}
createtemplate<-function(metofficefilein, e) {
  rte<-rast(metofficefilein)[[1]]
  rte<-crop(rte,e)
  rte<-aggregate(rte,25,fun="mean",na.rm=TRUE)
  return(rte)
}

createclimdf <- function(era5data, i, j) {
  dfout <- data.frame(obs_time=as.POSIXlt(time(era5data[[1]]), tz="UTC"),
        temp = .is(era5data[[1]])[i,j,],
        relhum = .is(era5data[[2]])[i,j,],
        pres = .is(era5data[[3]])[i,j,],
        swdown = .is(era5data[[4]])[i,j,],
        difrad = .is(era5data[[5]])[i,j,],
        lwdown = .is(era5data[[6]])[i,j,],
        windspeed = .is(era5data[[7]])[i,j,],
        winddir = .is(era5data[[8]])[i,j,],
        precip = .is(era5data[[9]])[i,j,])
  return(dfout)
}
getmetofficedata <- function(path, year, month, rte, what = "tasmin") {
  mtxt<-ifelse(month<10,paste0("0",month),paste0("",month))
  dm <- c(31,28,31,30,31,30,31,31,30,31,30,31)
  if (year%%4==0) dm[2]<-29
  if(what == "tasmin") {
    fi<-file.path(path,paste0("tasmin_hadukgrid_uk_1km_day_",year,mtxt,"01-",
              year,mtxt,dm[month],".nc"))
  } else {
    fi<-file.path(path,paste0("tasmax_hadukgrid_uk_1km_day_",year,mtxt,"01-",
              year,mtxt,dm[month],".nc"))
  }
  tc<-rast(fi)
  tc<-crop(tc,ext(rte))
  return(tc)
}


#' Run treetrunk model at 25km resolution using era5 raster data previously reprojected to crs 27700
#'
#' @param era5data - era5data output by `era5_process`
#' @param reqhgt - height above ground to model (metres)
#' @param depth - depth into trunk to model (metres from surface)
#' @param aspect - orientation of trunk location to monitor (degrees where 180 = south-facing aspect)
#' @param forestparams - vegetation parameters used by micropoint
#' @param groundparams - soil and ground parameters used by micropoint
#' @param treeparams - list of trunk parameters tradius = radiuse(metres),cs = spec heat capacity for each layer (J/Kg/K),
#' rho = denisty for each layer (kf/m^3), refl = trunk reflectance (0:1), trunk emmisivity = 0.97,
#' k = heat conductivity for each layer, surfwetness (0:1)
#' @param rte
#'
#' @return
#' @export
#'
#' @examples
runmodelgrid <- function(era5data, reqhgt = 1, depth = 0.02, aspect = 180,
                         forestparams, groundparams, treeparams, rte, output_airtemp=TRUE) {
  # create template temperature dataset for storing data and times
  tempout <- era5data[[1]] * 0
  tme<-terra::time(tempout)
  tempout <- .is(mask(tempout, rte))
  airtout <- tempout
  dms <- dim(tempout)
  # get northing/southing of each 1km grid cell
  ll<-.latslonsfromr(rte)
  for (i in 1:dms[1]) {
    for (j in 1:dms[2]) {
      v<-.is(rte)[i,j]
      if (is.na(v) == FALSE) {
        # create climate data.frame for grid cell
        weather <- createclimdf(era5data, i, j)
        weather$precip<-ifelse(weather$precip<0,0,weather$precip)
        # run microclimate model of canopy conditions from micropoint package
        microclim<-suppressWarnings(runpointmodel(weather,reqhgt,forestparams,
                                       paii=NA,groundparams,ll$lats[i,j],ll$lons[i,j]))
        airtout[i,j,]<-microclim$tair
        # prepare inputs for tree trunk model - at original hourly timesteps
        modelin<-preparevars(microclim,weather,ll$lats[i,j],ll$lons[i,j],mon=NA,timestep=NA)
        # run tree trunk model
        tempout[i,j,] <- runtreeemulator_v2(modelin,depth,aspect,treeparams)
      } # end if
    } # end j
  } # end i
  # Get raster of tree temps
  tempout <- .rast(tempout, rte)
  terra::time(tempout)<-tme
  # Get raster of air temps
  airtout<-.rast(airtout, rte)
  terra::time(airtout)<-tme
  if(output_airtemp==TRUE) tempout<-list('treetemp'=tempout,'airtemp'=airtout) else(tempout<-list('treetemp'=tempout))
  return(tempout)
}



### Run micropoint model of beneath canopy conditions in timesteps to allow for LAI variation
# Calls get_lai to model LAI through the year
# Runs. model in timesteps of 5 days when LAI is varying
# Runs treetrunk model for two depths inner and outer

runmodelgrid_pheno <- function(era5data, reqhgt = 1, depth_outer = 0.01, depth_inner=NA, aspect = 180,
                         forestparams, groundparams, treeparams, rte, output_airtemp=TRUE, MinLAI=0.5){

  # Get max min LAI
  MaxLAI<-forestparams$pai

  # create template temperature dataset for storing data and times
  tempout <- era5data[[1]] * 0
  tme<-terra::time(tempout)
  tempout <- .is(mask(tempout, rte))
  tempinr<-tempout
  airtout <- tempout
  dms <- dim(tempout)
  # get northing/southing of each 1km grid cell
  ll<-.latslonsfromr(rte)
  for (i in 1:dms[1]) {
    for (j in 1:dms[2]) {
      v<-.is(rte)[i,j]
      if (is.na(v) == FALSE) {
        ### 1 Microclimate modelling
        weather <- createclimdf(era5data, i, j)
        weather$precip<-ifelse(weather$precip<0,0,weather$precip)
        lai_seq<-get_lai(T=weather$temp,tme=weather$obs_time,lat=ll$lats[i,j], MaxLAI, MinLAI)

        # Calculate series of time steps to model microclimate
        # Get start/end of spring and fall from lai sequence
        if(length(unique(lai_seq))>1){
          lai_seq<-round(lai_seq,3)
          rlvals<-rle(lai_seq)
          sprstart<-sum(rlvals$lengths[1:(which(rlvals$values>MinLAI)[1]-1)])+1
          sprend<-sum(rlvals$lengths[1:(which(rlvals$values==MaxLAI)[1]-1)])
          fallstart<-length(lai_seq)-sum(rev(rlvals$lengths)[1:(which(rev(rlvals$values==MaxLAI))-1)])+1
          fallend<-length(lai_seq)-rlvals$lengths[length(rlvals$lengths)]

          # Calculate timesteps
          timestep<-24*5
          sprsteps<-(sprend-sprstart)/timestep
          spring_sseq<-seq(sprstart,sprstart+ceiling(sprsteps)*timestep, timestep)
          fallsteps<-(fallend-fallstart)/timestep
          fall_sseq<-seq(fallstart,fallstart+ceiling(fallsteps)*timestep, timestep)
          startseq<-c(1,spring_sseq,fall_sseq)
          endseq<-c(startseq[2:length(startseq)]-1,length(lai_seq))
        } else{ # No variation in LAI/PAI
          # timestep<-24*5
          # startseq<-seq(1,length(lai_seq),timestep)
          # endseq<-c(startseq[2:length(startseq)]-1,length(lai_seq))
          startseq<-1
          endseq<-length(lai_seq)
        }
        # timestep loop for changing phenology
        microclim_list<-list()
        for (t in 1:length(startseq)){
          start<-startseq[t]
          end<-endseq[t]
          steppai<-mean(lai_seq[start:end])
          forestparams$pai<- steppai
          microclim_list[[t]]<-suppressWarnings(runpointmodel(weather[start:end,],reqhgt,forestparams,
                                                    paii=NA,groundparams,ll$lats[i,j],ll$lons[i,j]))
        } # end timestep loop
        microclim<-do.call(rbind, microclim_list)
        airtout[i,j,]<-microclim$tair

        ### 2 Tree trunk Modelling
        # prepare inputs for tree trunk model - at original hourly timesteps
        modelin<-preparevars(microclim,weather,ll$lats[i,j],ll$lons[i,j],mon=NA,timestep=NA)
        # Run outer depth
        tempout[i,j,] <- runtreeemulator_v2(modelin,depth_outer,aspect,treeparams)
        # Run inner depth
        tempinr[i,j,] <- runtreeemulator_v2(modelin,depth_inner,aspect,treeparams)
      } # end if
    } # end j
  } # end i
  # Get raster of tree temps
  tempout <- .rast(tempout, rte)
  terra::time(tempout)<-tme
  tempinr <- .rast(tempinr, rte)
  terra::time(tempinr)<-tme

  # Get raster of air temps
  if(output_airtemp==TRUE){
    airtout<-.rast(airtout, rte)
    terra::time(airtout)<-tme
  } else airtout<-NA
  results_out<-list('outtemp'=tempout,'inrtemp'=tempinr,'airtemp'=airtout)
  return(results_out)
}




# Might need to add air temperature as well
postprocess <- function(era5data, tasmin, tasmax, treetemps) {
  # resample era5 data
  tempc<-resample(era5data[[1]], tasmin[[1]])
  print('Resampling era5 to 1km complete...')
  # convert met office data to hourly
  temp<-.rast(blendtempCpp(.is(tasmin),.is(tasmax),.is(tempc)), tempc)
  print('HadUK data converted to hourly...')
  # calculate anomoly
  dT <- temp - tempc
  # resample tree trunk temperature
  treetempsf <- resample(treetemps, dT)
  print('Resampling of model temperature to 1km complete...')
  treetempsf <- treetempsf + dT
  return(treetempsf)
}

# Calculate hourly anomaly between era5 and haduk at 1km resolution
calc_temp_anomaly <- function(era5data, tasmin, tasmax) {
  # resample era5 data
  tempc<-resample(era5data[[1]], tasmin[[1]])
  print('Resampling era5 to 1km complete...')
  # convert met office data to hourly
  temp<-.rast(blendtempCpp(.is(tasmin),.is(tasmax),.is(tempc)), tempc)
  print('HadUK data converted to hourly...')
  # calculate anomoly
  dT <- temp - tempc
  out<-list("dT"=dT,"weather_temp"=temp)
  return(out)
}

correct_temp<-function(temps,dT){
  # resample tree trunk temperature
  temps1km <- resample(temps, dT)
  print('Resampling of model tree temperature to 1km complete...')
  temps1km <- temps1km + dT
  # resample air temperature
  return(temps1km)
}


# Function to filter era5 inputs by months
filter_era5<-function(era5data,smonth,emonth){
  tme<-time(era5data[[1]])
  sel<-which(month(tme)>=smonth & month(tme)<=emonth)
  era5out<-lapply(era5data, function(x) {x[[c(sel)]]})
  return(era5out)
}

