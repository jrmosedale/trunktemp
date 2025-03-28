# Thermal conductivity cannot be reliably calculated at MC>25% therefore use 25% +?
# Gb = basic spec gravity
# Average moisture contents
# Oak: heatwood 64, sapwood 78
# Sitka spruce heartwood 41 sapwood 142
# MCfs fibre sat point moisture content

# Get trunk parameters from parameters table by species
trunk_parameters<-function(spp_params,sp,tradius, outerlayers=c(0.015,0.02,0.03,0.035)){
  # Calculate thermal properties from moisture content etc
  sw_mc<-spp_params['sw_mc',sp]
  hw_mc<-spp_params['hw_mc',sp]
  bk_mc<-spp_params['bk_mc',sp]
  fsp<-spp_params['fsp',sp]
  bk_fsp<-spp_params['bk_fsp',sp]

  # Calculate thermal conductivity
  swK<-calc_thcond(Gb=spp_params['Gb',sp],x=sw_mc)
  hwK<-calc_thcond(Gb=spp_params['Gb',sp],x=hw_mc)
  bkK<-calc_thcond(Gb=spp_params['bk_Gb',sp],x=bk_mc)

  # Calculate spec heat capacity (sapwood outer and inner other layers) kJ/Kg/K
  swcph<-calc_sph(sw_mc,cp0=NA,T=283,fsat=fsp)
  hwcph<-calc_sph(hw_mc,cp0=NA,T=283,fsat=fsp)
  bkcph<-calc_sph(bk_mc,cp0=NA,T=283,fsat=bk_fsp)

  # Calculate wood density (green)
  swrho<-calc_density(Gb=spp_params['Gb',sp],x=sw_mc)
  hwrho<-calc_density(Gb=spp_params['Gb',sp],x=hw_mc)
  bkrho<-calc_density(Gb=spp_params['bk_Gb',sp],x=bk_mc)

  # Calculate layer parameters etc
  nlyr<-8
  nsegs<-16

  # Width of layers - probes shallow = 2 (2.5 cm) deep = layer 3 (5cm)
  # Fixed outer layers then calculate others
  #outerlayers<-c(0.015,0.02,0.03,0.035)
  #outerlayers<-c(0.005,0.025,0.03,0.035)
  thick<-sqrt(c(5:8))
  thick <- (thick/ sum(thick)) * (tradius-sum(outerlayers))
  layer_widths<-c(outerlayers,thick) # sum to tradius
  if(round(sum(layer_widths),2)!=round(tradius,2)) stop("Layer widths do NOT sum to tree radius!!!")

  if(sp=='nsp'){
    cs <- c(bkcph,swcph,swcph,rep(hwcph,(nlyr-3)))
    rho <- c(bkrho,swrho,swrho,rep(hwrho,(nlyr-3)))
    k <- c(bkK,swK,swK,rep(hwK,(nlyr-3)))
  }
  if(sp %in% c('oak','sycamore')){
    cs <- c(bkcph,swcph,rep(hwcph,(nlyr-2)))
    rho <- c(bkrho,swrho,rep(hwrho,(nlyr-2)))
    k <- c(bkK,swK,rep(hwK,(nlyr-2)))
  }

  treeparams <- list(
    tradius = tradius,
    layer_widths=layer_widths,
    cs = cs,
    rho = rho,
    refl = spp_params[["bk_swref",sp]],
    em = spp_params[["bk_tem",sp]],
    k = k,
    surfwet = 1 )

  return(treeparams)
}

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
#library(sigmoid)
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



# Spec gravity 5g/cm3 density = ~ 5 spec g
# p0=density of oven dry wood (320-720kg/m3)
# Oak Gb0.6 green, 0.49 12%
# Sitka spruce Gb 0.37 green, 0.4 12%
# calc_thcond(Gb=0.6,x=25)
# calc_thcond(Gb=0.37,x=25)
calc_thcond<-function(Gb,x=25,MCfs=30){
  if(x>25) warning("Calculation of thermal conductivity unreliable at moisture contents > 25%")
  A<-0.01864
  B<-0.1941
  C<-0.004064
  # Convert Spec Gravity to moisture content x
  # Assumes total vol shrinkage S0 estimated from Gb
  Gx<-Gb/(1-0.265*Gb*(1-x/MCfs)) #Eq 4.13
  # Calc conductivity from Spec gravity at MCx
  k<-Gx*(B+C*x)+A
  return(k)
}

# Calculate Specific heat capacity at different moisture contents using
# The heat capacity of wood depends on the temperature and moisture content
# of the wood but is practically independent of density or species. (Glass & Zelinka)
# Output units:
# x = Moisture content T = temp, cp0 if provide = spec heat cap of dry wood, fsat = moisture content at fibre saturation
# for (n in seq(10,90,10)) print(paste(n,calc_sph(n)))

calc_sph<-function(x,cp0=NA,T=283,fsat=30){
  if (is.na(cp0)) cp0<-0.1031+0.003867*T
  # Calculate correction factor Ac for <fibre sat moisture content
  cpw<-4.18 # Water Sp H Cap
  b1<- -0.06191
  b2<-0.000236
  b3<- -0.000133
  # Moisture content < fibre sat
  if(x<=fsat){
    Ac<-x*(b1+b2*T+b3*x) # Correction factor for <fibre sat moisture content
    cpx<-(cp0+cpw*x/100) / (1+x/100) + Ac
  }
  if(x>fsat){
    Ac<-fsat*(b1+b2*T+b3*fsat)
    cpfs<-(cp0+cpw*fsat/100) / (1+fsat/100) + Ac
    cpx<- cpfs*(fsat/x) + (cpw*(1-fsat/x))
  }
  return(cpx)
}

# Return density of wood at given moisture content in kg/m3 units
# rhow in cm/g3 =1
# Gb wood basic specific gravity (oven dry mass, green volume in g/cm3?)
# fs = moisture content at fibre sat (%)
calc_density<-function(Gb,x,MCfs=30, rhow=1){
  Gx<-Gb/(1-0.265*Gb*(1-x/MCfs))
  rho<-rhow*Gx*(1+x/100)
  return(rho*1000)
}
# calc_density(Gb,x=90)

