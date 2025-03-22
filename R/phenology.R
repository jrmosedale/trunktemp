# Phenology models for application
# Depend on temperature, day length (and water availabilityy - ignored)
# sigmoidal calculation
# where start pt x=0 b=(1/y)-1
#b<-(1/MinLAI)-1
# at end x<-100 and y>3
#x<-100
#y<-0.99
#k<- -log( (1/y-1)/b) / x
#x<-seq(0,100,1)
#y<-1 / (1+(b*exp(-k*x)))

# Thresholds / Variables
# Max LAI (PAI - something?)
# Day length 10.5hrs
# Temperature 8-10C
# Intitial LAI growth (linera or exponential?)
# Spatial variation?

# References
# Heiskanen et al 2012 Fig 2a ~ 0.5 LAI variation In coniferous stands, the standard deviation was 5–11% during the study period Finland

################### Debbie's model ###################
# Calculate budburst
# Calculate sensecence
# Calculate changing lai

#weather <- createclimdf(era5data, 10, 10)
#lai<-get_lai(T=weather$temp,tme=weather$obs_time,lat=51, MaxLAI=3,MinLAI=0.5, sprg_sigmoid=TRUE)
# plot(lai)
#tme[which(lai>MinLAI)[1]]
#pklai<-tme[which(lai==MaxLAI)[1]]
#tme[which(lai<MaxLAI & tme>pklai)[1]]
#tme[which(lai==MinLAI & tme>pklai)[1]]

# Simple budburst and senescance model using daily dd
get_lai<-function(T,tme,lat, MaxLAI=3, MinLAI=0.5,sprg_sigmoid=FALSE,fall_sigmoid=FALSE){
  # 1 Calculate budburst (dd model)Fu model) using gdd from 1 Jan
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

t0<-now()
budburst<-gb*0
LFstart<-gb*0
LFend<-gb*0
for (i in 1:dms[1]) {
  for (j in 1:dms[2]) {
    v<-.is(rte)[i,j]
    if (is.na(v) == FALSE) {
      # create climate data.frame for grid cell
      weather <- createclimdf(era5data, i, j)
      lai<-get_lai(T=weather$temp,tme=weather$obs_time,lat=ll$lats[i,j], MaxLAI=3,MinLAI=0.5)
      budburst[i,j]<-yday(tme[which(lai>MinLAI)[1]])
      pklai<-yday(tme[which(lai==MaxLAI)[1]])
      LFstart[i,j]<-yday(tme[which(lai<MaxLAI & yday(tme)>pklai)[1]])
      LFend[i,j]<-yday(tme[which(lai==MinLAI & yday(tme)>pklai)[1]])
    }
  }}
print(now()-t0)
plot(budburst,main='budburst')
plot(LFstart,main='LFstart')
plot(LFend,main='LFend')





