#' Simple budburst and senescance model using daily dd
#' Calculate change in PAI/LAI
#' Depend on temperature & day length (this version ignores water availability as parameter)
#' Canopy evolution in spring and fall described by either a sigmoid curve or linear line
#' @param tme - time series
#' @param lat - latitude of location
#' @param MaxLAI - maximum LAI/PAI
#' @param MinLAI - minimum LAI/PAI
#' @param sprg_sigmoid - TRUE if sigmoid curve to be fitted to spring change (otherwise linear)
#' @param fall_sigmoid - TRUE if sigmoid curve to be fitted to autumn change (otherwise linear)
#'
#' @return hourly timeseries of LAI/PAI values ranging between MinLAi and MAxLAI
#' Only daily change is modelled so values within the same day are constant
#' @import sigmoid
#' @export
#'
#' @examples
#' lai<-get_lai(T=weather$temp,tme=weather$obs_time,lat=51, MaxLAI=3,MinLAI=0.5, sprg_sigmoid=TRUE)
get_lai<-function(T,tme,lat, MaxLAI=3, MinLAI=0.5,sprg_sigmoid=FALSE,fall_sigmoid=FALSE){
  # 1 Calculate budburst (dd model) Fu model) using gdd from 1 Jan
  Tb<- -5
  ADDcrit<-591
  bb_to_fullleaf<-85 # days between budburst and full leaf

  dayTmax<-tapply(T,INDEX=yday(tme),FUN=max)
  dayTmin<-tapply(T,INDEX=yday(tme),FUN=min)
  dd<-cumsum(((dayTmax-dayTmin)/2)-Tb)
  bb<-which(dd>ADDcrit)[1]

  # 2 Calculate senecance using Delpierre model
  # Parameters - Quercus
  Pstart<-14 # 14.5 max day length at which senscence are effective DBF=1.5
  Tb<-26.5 # 26.5 maximum temperature at which senescence processes are effective DBF=28.5
  x<-2
  y<-0 # DBF=2
  Ycrit<-10178 #threshold for sum(Rsen) to reach DBF=8268

  # Calculate rates of senescence
  dayTmean<-tapply(T,INDEX=yday(tme),FUN=max)
  jdays<-mesoclim:::.jday(as.POSIXlt(tme[c(seq(1,length(tme),24))]))
  daylength<-mesoclim::daylength(jdays,lat)
  doy<-c(1:length(daylength))
  Rsen<-ifelse(doy>180 & daylength<Pstart & dayTmean<Tb , (Tb-dayTmean)^x * (1-(daylength/Pstart))^y, 0)
  Ssen<-cumsum(Rsen)
  # Get start and end days of senescence
  Dstart<-which(Rsen>0)[1]
  Y90<-which(Ssen>Ycrit)[1]

  # Calculate increase/decrease in LAI -  linear or sigmoid?
  lai<-rep(MinLAI,length(doy))
  if(sprg_sigmoid==FALSE){
    LAIinc<-as.numeric((MaxLAI-MinLAI)/(bb_to_fullleaf+1)) # assume increases over 85 days to max
    for(n in bb:(Y90-1)) lai[n]<-ifelse(lai[n-1]<MaxLAI,min(MaxLAI,lai[n-1]+LAIinc),lai[n-1])
  }
  if(sprg_sigmoid==TRUE){
    fullleaf<-bb+bb_to_fullleaf
    x<-seq((bb-(bb+8)),(fullleaf-(fullleaf-8)),length=fullleaf-bb)
    lai[bb:(fullleaf-1)]<-((MaxLAI-MinLAI)*sigmoid(x))+MinLAI
    lai[fullleaf:(Dstart-1)]<-MaxLAI
  }
  if(fall_sigmoid==FALSE){
    LAIdec<-as.numeric((MaxLAI-MinLAI)/(Y90-Dstart+1))
    for(n in Dstart:length(lai)) lai[n]<-lai[n]<-ifelse(lai[n-1]>MinLAI,max(MinLAI,lai[n-1]-LAIdec),MinLAI)
  }
  if(fall_sigmoid==TRUE){
    x<-seq((Dstart-(Dstart+8)),(Y90-(Y90-8)),length=Y90-Dstart+10)
    lai[Dstart:(Y90-1+10)]<-rev(((MaxLAI-MinLAI)*sigmoid(x))+MinLAI)
  }
  # plot(lai)

  # Convert back top hourly lai
  lai_hrly<-rep(lai,each=24)
return(lai_hrly)
}






