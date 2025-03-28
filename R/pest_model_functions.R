########################################################################
### General functions
########################################################################
# calculate day of year from hourly (daily min/max calculated) temperature data using a base temp and threshold
# temps - hourly vector of temperatures
# base_temp - the temperature from which accumulated degrees should be calculated
# acc_gdd - the threshold of growing degree days at which indicator is reached
# offset = number of hours from start to skip
gdd_micro <- function(temps, tme, Tbase, acc_gdd=NA, offset = 0) {
  out<-NA
  if(!all(is.na(temps))){
    if(offset>0) {
      temps[1:offset] <- 0 # turn offsetted values to 0
    }
    tstep<-as.numeric(tme[2])-as.numeric(tme[1]) # in secs
    tdiv<-60*60*24/tstep
    temprate<-(temps-Tbase)
    temprate<-ifelse(temprate<=0,0,temprate / tdiv)
    gdd<-cumsum(temprate)
    # If threshold dd provided returns Day of Year it is attained
    if(!is.na(acc_gdd)){
      sel<- which(gdd > acc_gdd)
      if(length(sel)>0) out <- yday(tme[min(sel)])
    }
    # If no threshold provided returns max gdd
    if(is.na(acc_gdd)) out<-max(gdd)
  }
  return(out)
}

########################################################################
### Asian longhorn GDD model - inner and outer barklocations
# No variation in DD requirements
# Sets a max length of generation = 4 yrs - does not continue beyond this point
# Just follows the first cohort of eggs laid - not G2, G3 etc
# Outputs proportion of beetles with different generation times in years: 1,2,3,4,>4
# Vary number of instars based on parameterization molt value data
# Min L instars = 6 (stages 1-7 +3), Max L instar =11 (stages 1-12 +3)
# therfore 6 different life cycle groups based on number of instars
# UCT and LCT = absolute limits to development
########################################################################
# Convert Degree hr data to Degree Day
dhr_to_dd<-function(dhr,tme){
  if(length(dhr)%%24!=0) stop("Length of data in dhr_to_dd NOT multiple of 24!!!")
  if(hour(tme[1])!=0 || hour(tme[length(tme)])!=23) warning("Check datetime data in dhr_to_dd - may not be whole timeseries of days!!!")
  #dd<-zoo::rollsum(dhr,24,align="right",fill=NA)
  dd<-as.numeric(by(dhr, yday(tme), FUN=sum))
  doy<-unique(yday(tme))
  return(list(dd=dd,doy=doy))
}

# Calculate pupahold time from 5 day running mean - use outer temp as pupates in nner cambium - return Day of Year
get_pupagate<-function(Tout,tme,pupahold=14){
  doys<-yday(tme)
  dayTmax<-tapply(Tout,INDEX=doys,FUN=max)
  roll5mean<-as.numeric(zoo::rollmean(dayTmax,5,align = 'right', na.pad = TRUE))
  pktempday<-which(roll5mean==max(roll5mean,na.rm=TRUE))
  pupagate<-max(pktempday,200) # sets to 19 July if peak temperature before this??
  pupagate<-pupagate + pupahold
  return(pupagate)
}
# Calculate dd for a specific stage - returns daily vector of values matching tme
calc_stage_dd<-function(s,parameters,Tout,Tinr,tme){
  tstep<-as.numeric(tme[2])-as.numeric(tme[1])
  daysecs<-24*60*60
  tmerate<-tstep/daysecs

  LCT<-parameters$LCT[which(parameters$stage_number==s)]
  UCT<-parameters$UCT[which(parameters$stage_number==s)]
  depth<-parameters$depth[which(parameters$stage_number==s)]
  if(depth=="outer") T<-Tout
  if(depth=="inner") T<-Tinr

  # Calc hourly rT
  rT<-calc_rT(T,LCT,UCT,tmerate)
  # Get day total (degreeday)
  dd<-dhr_to_dd(rT,tme)["dd"]
  return(dd)
}

# FUNCTION to calculate development rate
calc_rT<-function(T,LCT,UCT,tmerate){
  # Calculate dev rate
  rT<-ifelse(T<UCT & T>LCT, (T-LCT) * tmerate, 0)
  # Accumulative dev rate and completion day
  # gdd<-cumsum(rT)
  return(rT)
}

# Model group development over a year
year_development<-function(gstage,gdd,mstage,parameters,dd_matrix,pgate){
  # Continue until no more development in the year
  stgend_v<-c()
  complete_v<-c()
  while(cumsum(dd_matrix[,min(15,gstage)])[nrow(dd_matrix)]>0 & gstage<=15){
    dd<-gdd+cumsum(dd_matrix[,gstage])
    # Does it reach end of stage by year end?
    reqDD<-parameters$DD[which(parameters$stage_number==gstage)]

    if(dd[length(dd)] >=reqDD) completes_stage<-TRUE else completes_stage<-FALSE

    if(completes_stage){           # Completes stage
      # Does it pupate?
      if(gstage!=mstage){           # Does not pupate
        stgend_doy<-which(dd>=reqDD)[1]
        complete_v<-c(complete_v,gstage)
        stgend_v<-c(stgend_v,stgend_doy)
        #print(paste("Stage",gstage,"completed on",stgend_doy,"day of year"))
        gstage<-gstage+1
        if(gstage<=15) dd_matrix[1:stgend_doy,gstage]<-0
        if(gstage<=15) gdd<-dd[stgend_doy]-reqDD

      } else {                        # Pupates/ Diapause
        #print("Time to pupate!")
        stgend_doy<-which(dd>=reqDD)[1]
        complete_v<-c(complete_v,gstage)
        stgend_v<-c(stgend_v,stgend_doy)
        #print(paste("Stage",gstage,"completed on",stgend_doy,"day of year"))
        gstage<-13
        dd_matrix[1:stgend_doy,gstage]<-0
        gdd<-dd[stgend_doy]-reqDD
        # Season ends if later than pupagatepate
        if(!is.na(pgate) & stgend_doy>=pgate){  # Enters diapause rather than pupate
          #print("ready to pupate but wait till next year!")
          dd_matrix[,gstage]<-0 # Keep gdd as is so pupates at start of next year
        }
      }
    }

    if(!completes_stage){             # Does not complete stage
      #print(paste("Does not complete stage",gstage))
      gdd<-dd[length(dd)]
      dd_matrix[,gstage]<-0
    }
  } # end while

  return(list(gstage=gstage,gdd=gdd,complete_v=complete_v,stgend_v=stgend_v))
} # end function


run_alh<-function(Tout_list,Tinr_list,tme_list,parameter_file,numgrps=6,use_pgate=TRUE){
  # Checks first year of Tout_list - if all NA then skips calculations
  if(!all(is.na(Tout_list[[1]]))){

    # Set up parameters
    tstep<-as.numeric(tme[2])-as.numeric(tme[1])
    daysecs<-24*60*60
    tmerate<-tstep/daysecs
    startyr<-as.numeric(names(tme_list)[1])

    parameters<-read.csv(parameter_file)
    parameters[,3:7]<-round( parameters[,3:7],3)
    numstages<-nrow(parameters)
    startyr<-as.numeric(names(Tout_list)[1])
    e20<-16 # Days from emergence to oviposition
    pupahold<-14  # Days after max rolling 5 day summer temperature for initiation of pupal development (otherwise diapause)
    eggstart<-230 # 18 August day for starting eggs laid
    maxflightdoy<-345 # if. used...
    year_doys<-as.numeric(unlist(lapply(tme_list,function(x){length(x)/24})))

    # Vectors of keeping track of each group development	 - rtest with single group
    grp_stage<-rep(1,numgrps) # What is currenmt stage of each group
    grp_stagedd<-rep(0,numgrps) # store end of year dd of current stage for each group
    grp_moultstg<-which(parameters$moult<1 & parameters$moult>0)
    stageend_m<-matrix(data=rep(NA,numgrps*numstages),ncol=numgrps,nrow = numstages) # store time when stage completed for each group

    # For each year
    for(yr in names(tme_list)){
      # print(yr)
      # Get year's data
      tme<-tme_list[[yr]]
      Tout<-Tout_list[[yr]]
      Tinr<-Tinr_list[[yr]]
      days_in_yr<-length(tme)/24
      if(use_pgate) pgate<-get_pupagate(Tout,tme,pupahold) else pgate<-NA # day of year

      # Calculate DD development for every possible stage (15 stages)
      dd_matrix<-matrix(data = rep(0,numstages*days_in_yr), nrow = days_in_yr, ncol = numstages)
      for (s in 1:15){
        dd<-calc_stage_dd(s,parameters,Tout,Tinr,tme)
        dd_matrix[,s]<-dd$dd
      }

      # If first year - no development until eggstart doy
      if(as.numeric(yr)==startyr) dd_matrix[1:eggstart-1,]<-0

      # Run development of each group for the year
      for(g in 1:numgrps){
        #print(g)
        gstage<-grp_stage[g]
        gdd<-grp_stagedd[g]
        mstage<-grp_moultstg[g]
        if(gstage<=15){
          yrdev<-year_development(gstage,gdd,mstage,parameters,dd_matrix,pgate)
          #print(paste(yrdev))
          grp_stage[g]<-yrdev$gstage
          grp_stagedd[g]<-yrdev$gdd
          # Records total doys since oviposition
          if(!is.null(yrdev$stgend_v)){
            if (yr!=startyr) prevdoy<-sum(year_doys[1:(as.numeric(yr)-startyr)]) else prevdoy<-0
            stageend_m[yrdev$complete_v,g]<-yrdev$stgend_v - eggstart + prevdoy
          }
        }
      } # for group
    } # end year


    #print(stageend_m)

    # Vector of outputs for later conversion to raster
    if(all(is.na(stageend_m[15,]))) grps_completed<-0 else grps_completed<-length(which(!is.na(stageend_m[15,])))
    if(grps_completed==0) pop_completed<-0 else pop_completed<-max(parameters$moult[grp_moultstg][grps_completed])
    days_to_complete<-stageend_m[15,]
    days_to_pupa<-rep(NA,numgrps)
    for(g in 1:numgrps){days_to_pupa[g]<-stageend_m[grp_moultstg[g],g] }

  } else{ ################ if all NA then return NA values. { ################
    pop_completed<-NA
    days_to_complete<-rep(NA,numgrps)
    days_to_pupa<-rep(NA,numgrps)
  }

  output<-c(pop_completed,days_to_complete,days_to_pupa)
  names(output)<-c("pop_completed",
                   paste0("days_to_complete_",seq(1,6,1)),
                   paste0("days_to_pupa_",seq(1,6,1)))

  return(output)
} # end function


########################################################################
### Ips typographus non-linear egg to pupae model - outer bark
# see IpsTypographus_test.R
########################################################################
library(mesoclim)
run_ips<-function(T,Tair, lat, tme){
  # Set up
  emerge<-NA
  lay<-NA
  adult<-NA
  lay_g2<-NA
  lay_g3<-NA
  adult_g2<-NA
  adult_g3<-NA
  score<-NA
  emerge_doy<-NA
  lay_doy<-NA
  adult_doy<-NA
  lay_g2_doy<-NA
  adult_g2_doy<-NA
  lay_g3_doy<-NA
  adult_g3_doy<-NA

  if(!all(is.na(T))){

    # Calculate when latest diapause can occur - >1 October & rolling 7 day mean of daily Tmax  < 17.5C COPIED from Dmicans
    # Should be related to daylength?
    doys<-yday(tme)
    dayTmax<-tapply(T,INDEX=doys,FUN=max)
    sdoy<-274-6 # 1 Oct - 7 days
    roll7mean<-as.numeric(zoo::rollmean(dayTmax[sdoy:length((dayTmax))],7,align = 'right', na.pad = TRUE))
    diapause_doy<-which(roll7mean<17.5)[1]+sdoy-1

    # Impose diapause as end of developmental year
    tme<-tme[which((doys<diapause_doy))]
    T<-T[which((doys<diapause_doy))]
    Tair<-Tair[which((doys<diapause_doy))]

    # Hourly seq
    tstep<-as.numeric(tme[2])-as.numeric(tme[1])
    daysecs<-24*60*60
    tmerate<-tstep/daysecs
    tmedoys<-yday(tme) # day of year
    tmelength<- mesoclim::daylength(tmedoys,lat) # length of day
    # Daily seq
    Tdaymax<-tapply(Tair,INDEX=tmedoys,FUN=max)
    doys<-as.numeric(names(Tdaymax))
    # Hourly seq of day maxima
    tmeTairmx<-rep(as.numeric(Tdaymax),each=24)

    ## 1st Generation
    # Adult emergence - depends on max daily Tair calculated in daily steps
    Tair_emg<-14.5
    Tbase<-8.3
    ddend<-53
    start_doy<-66 # 7 March

    Tdaymax<-tapply(Tair,INDEX=tmedoys,FUN=max)
    doys<-as.numeric(names(Tdaymax))
    rT<-ifelse(Tdaymax>Tbase & doys>=start_doy, (Tdaymax-Tbase), 0)
    gdd<-cumsum(rT)
    emerge_doy<-doys[which(gdd>ddend & Tdaymax>Tair_emg)][1]
    #   emerge<-tme[which(Tair>Tair_emg)[1]]

    # Egg lay oviposition day- dependent on max daily Tair
    Tbase<-8.3
    ddend<-103.6
    rT<-ifelse(Tdaymax>Tbase & doys>=emerge_doy, (Tdaymax-Tbase), 0)
    gdd<-cumsum(rT)
    lay_doy<-doys[which(gdd>ddend)][1]

    # Egg to Pupa parameters
    a<- 0.10377
    b<- 9.60332
    g<- -0.02245
    Tmx<- 40
    #Topt<-30.4
    #Tlow<-5.8
    #Thigh<-38.9

    # Non-linear rate model
    #   ifelse(x<Tlow|x>Thigh,0,exp(a*T)-exp(a*Tmx-(Tmx-T)/b) + g )
    rT<-exp(a*T)-exp(a*Tmx-(Tmx-T)/b) + g
    rT<-ifelse(rT<0 ,0,rT)
    rT<-ifelse(tmedoys<=lay_doy,0,rT)
    rT<-rT*tmerate
    adult<-tme[which(cumsum(rT)>=1)[1]]
    adult_doy<-as.numeric(strftime(adult, format = "%j"))

    ## 2nd Generation
    # Emerge & lay G2
    lay_g2<-tme[which(cumsum(rT)>=1.6 & tmelength>=14.5 & tmeTairmx>=14.5 & tme>adult)][1]
    lay_g2_doy<-as.numeric(strftime(lay_g2, format = "%j"))
    # Egg to Pupa
    rT<-exp(a*T)-exp(a*Tmx-(Tmx-T)/b) + g
    rT<-ifelse(rT<0 ,0,rT)
    rT<-ifelse(tme<=lay_g2,0,rT)
    rT<-rT*tmerate
    adult_g2<-tme[which(cumsum(rT)>=1)[1]]
    adult_g2_doy<-as.numeric(strftime(adult_g2, format = "%j"))

    ## 3rd Generation
    lay_g3<-tme[which(cumsum(rT)>=1.6 & tmelength>=14.5 & tmeTairmx>=14.5 & tme>adult_g2)][1]
    lay_g3_doy<-as.numeric(strftime(lay_g3, format = "%j"))

    # Egg to Pupa
    rT<-exp(a*T)-exp(a*Tmx-(Tmx-T)/b) + g
    rT<-ifelse(rT<0 ,0,rT)
    rT<-ifelse(tme<=lay_g3,0,rT)
    rT<-rT*tmerate
    adult_g3<-tme[which(cumsum(rT)>=1)[1]]
    adult_g3_doy<-as.numeric(strftime(adult_g3, format = "%j"))
  }
  # Output
  output<-c(emerge_doy,lay_doy,adult_doy,lay_g2_doy,adult_g2_doy,lay_g3_doy,adult_g3_doy)
  names(output)<-c("ips_emerge_doy","ips_lay_doy","ips_adult_doy","ips_g2_lay_doy","ips_g2_adult_doy","ips_g3_lay_doy","ips_g3_adult_doy")

  # Score based on number of complete generations
  if(!all(is.na(output))) score<-as.numeric(which(is.na(output[c("ips_adult_doy","ips_g2_adult_doy","ips_g3_adult_doy")]))[1])-1
  output<-c(output,"generations_complete"=score)
  return(output)
}

########################################################################
### Dendroctonus micans
# Host spruce - outer bark usually low down tree bole often on north sides
# Generations: uni or semivoltine, with typical generation times of UK: 10-18 months; Europe: 1-3 years
# Model approach: multi-year model to estimate voltinism and time for a single generation to develop using linear degree day models
# (Gent et al., 2017) with fitted lower temperature limits beginning with parental generation as overwintering pupa.

# Tair measured from weather
# Assume starting G overwinter as pupa
# Assume Tair thresholds must be met for adult emergence and laying
# Impose diapause when Tback< threshold for 7 consec days (OR when L5 complete??)

# Output
# Completes L5
# Complete Pupate (to adult)

########################################################################
# Version where linear decline between Topt and Tmax - following suggestion of Gent
# Option run yr by yr and simple record which stage completes / proportion of stage completed
# OR run across multiple years - see run_dendrocutonus_multiyear
# T<-values(treetemps1km)[10000,];Tair<-T; tme<-time(airtemps1km)
run_dendroctonus<-function(T,Tair,tme){
  # Setup and parameters
  lay<-NA
  hatch<-NA
  pupate<-NA
  adult<-NA
  dev_complete<-NA
  score<-NA

  if(!all(is.na(T))){
    tstep<-as.numeric(tme[2])-as.numeric(tme[1])
    daysecs<-24*60*60
    tmerate<-tstep/daysecs
    #tmedoys<-yday(tme) # day of year

    # Calculate when latest diapause can occur - >1 October & rolling 7 day mean of daily Tmax  < 17.5C
    doys<-yday(tme)
    dayTmax<-tapply(T,INDEX=doys,FUN=max)
    sdoy<-274-6 # 1 Oct - 7 days
    roll7mean<-as.numeric(zoo::rollmean(dayTmax[sdoy:length((dayTmax))],7,align = 'right', na.pad = TRUE))
    diapause_doy<-which(roll7mean<17.5)[1]+sdoy-1

    # Impose diapause as end of developmental year
    tme<-tme[which((doys<diapause_doy))]
    T<-T[which((doys<diapause_doy))]
    Tair<-Tair[which((doys<diapause_doy))]


    ### Assume overwinter as pupa so must develop to adult - this pushes flight emergence back by 60-100 days (March-April)
    Tbase<-7.23
    ddend<-126
    rT<-ifelse(T>=Tbase, (T-Tbase) * tmerate, 0)
    gdd<-cumsum(rT)
    adult_g0<-tme[which(gdd>ddend)][1] # day when development complete from overwintering as pupa

    ### Point of egg laying based on air temp
    Tair_emg<-14
    days_to_lay<-12
    lay<-tme[which(Tair>Tair_emg & tme>=adult_g0)[1]+(days_to_lay/tmerate)] # first day when adult might emerge

    ### Eggs
    Tbase<-7.39
    ddend<-153
    rT<-ifelse(T>Tbase & tme>lay , (T-Tbase) * tmerate, 0)
    gdd<-cumsum(rT)
    hatch<-tme[which(gdd>ddend)][1]

    ### L1 to L4 - modified version
    Tbase<-6.64
    Tmx<-25
    Topt<-21
    ddend<-527
    # To calculate decline if T above Topt
    m<- (-Topt-Tbase) / (Tmx-Topt)
    n<- -m*Tmx
    # Calculate dev rate
    rT<-ifelse(T>Tbase & T<=Topt & tme>hatch, (T-Tbase) * tmerate, 0)
    rT<-ifelse(T>Topt & T<=Tmx & tme>hatch, (m*T+n) * tmerate, rT)
    # Accumulative dev rate and completion day
    gdd<-cumsum(rT)
    pupate<-tme[which(gdd>ddend)][1]

    ### Pupa to adult
    Tbase<-7.23
    ddend<-126
    rT<-ifelse(T>=Tbase & tme>pupate, (T-Tbase) * tmerate, 0)
    gdd<-cumsum(rT)
    adult<-tme[which(gdd>ddend)][1]

    # Developmental completion 0 to 1 of final phase
    dev_complete<-as.numeric(gdd[length(gdd)]/ddend)
    dev_complete<-ifelse(dev_complete>=1,1.0,dev_complete)
    if(length(dev_complete)==0) dev_complete<-NA

  }

  # Get day of year when growth phases complete
  output<-c(lay,hatch,pupate,adult)
  output<-as.numeric(strftime(output, format = "%j"))
  #if(out=="doy") else output<-as.character(output)
  names(output)<-c("dnd_lay_doy","dnd_hatch_doy","dnd_pupate_doy","dnd_adult_doy")

  # Score based on last phase complete
  complete<-which(!is.na(output))
  if (length(complete)==0) score<-NA else score <-max(complete)-1

  output<-c(output,"dnd_stages_complete"=score,"dnd_incomplete_stage"=dev_complete)

  return(output)
}
#results<-run_dendroctonus(T=values(treetemps1km)[10000,],Tair=values(treetemps1km)[10000,], tme<-time(airtemps1km))

# Convert matrix to spatraster of results
get_results_rast<-function(results,template_r){
  template_r<-ifel(is.na(template_r),NA,0)
  template_r<-rep(template_r,dim(results)[1])
  names(template_r)<-rownames(results)
  results_r<-setValues(template_r,t(results))
  return(results_r)
}

########################################################################
### Agrillus biguttatus - host = oak
# Tobark = outer bark temperature
# Tibark = inner bark temperatures
# Tair = air temperatures
# tme = datetime
# Imposes diapuse at completion of L4 (pupate date)
# Diapause also imposed based on Temp/Day of Year (may be prior to pupate date)
########################################################################
#Tob<-values(outtemps1km)
#Tob<-asplit(Tob,1) # array of lists of timeseries
#Tib<-values(inrtemps1km)+1
#Tib<-asplit(Tib,1) # array of lists of timeseries
#Ta<-values(airtemps1km)
#Ta<-asplit(Ta,1)
#tme<-time(airtemps1km)

#agresults<-mapply(FUN=run_agrillus,Tobark=Tob,Tibark=Tib,Tair=Ta,MoreArgs=list(tme=tme)) # output as day of year
#agresults.r<-get_results_rast(agresults,template_r=treetemps1km[[1]])
#plot(agresults.r)
#Tobark<-as.numeric(unlist(Tob[1000]))+3
#Tibark<-as.numeric(unlist(Tib[1000]))
#Tair<-as.numeric(unlist(Ta[1000]))+4
#tme<-time(airtemps1km)

run_agrillus<-function(Tobark, Tibark, Tair, tme){
  # Setup and parameters
  eclosion<-NA
  lay<-NA
  hatch<-NA
  pupate<-NA
  dev_complete<-NA
  score<-NA
  ddend<-c('eclosion'=76,'hatch'= 157.1,'pupate'=615.9 )
  larval_gdd<-0

  if(!all(is.na(Tobark))){
    tstep<-as.numeric(tme[2])-as.numeric(tme[1])
    daysecs<-24*60*60
    tmerate<-tstep/daysecs

    # Calculate daily max for air and inner bark temperatures
    doys<-yday(tme)
    dayTibmx<-tapply(Tibark,INDEX=doys,FUN=max) # daily inner bark max temp
    dayTairmx<-tapply(Tair,INDEX=doys,FUN=max) # daily air max temp

    # Calculate when latest diapause can occur - 1 Oct or later if temperatures elevated
    Tdiap<-11.9
    sdoy<-274-6 # Calculate rolling 7 day mean of max daily T from 1 October
    roll7mean<-as.numeric(zoo::rollmean(dayTibmx[sdoy:length((dayTibmx))],7,align = 'right', na.pad = TRUE))
    diapause_doy<-which(!is.na(roll7mean) & roll7mean<Tdiap)[1]+sdoy-1

    # Impose diapause as end of developmental year
    tme<-tme[which((doys<diapause_doy))]
    Tair<-Tair[which((doys<diapause_doy))]
    Tibark<-Tibark[which((doys<diapause_doy))]
    Tobark<-Tobark[which((doys<diapause_doy))]

    # 1 Assume overwinter as final larval instar subsequently pupates in outer bark during spring
    Tbase<-15.1
    rT<-ifelse(Tobark>=Tbase, (Tobark-Tbase) * tmerate, 0)
    ecl_gdd<-cumsum(rT)
    eclosion<-tme[which(ecl_gdd>ddend['eclosion'])][1] # day when development complete from overwintering as pupa

    # Flight time and laying
    Tair_emg<-12 # min daily max Tair required for flight
    days_to_lay<-12 # days for egg laying
    lay<-tme[which(Tair>Tair_emg & tme>=eclosion)[1]+(days_to_lay/tmerate)] # first day when adult might emerge

    # Eggs in inner bark
    Tbase<-12.1
    rT<-ifelse(Tibark>Tbase & tme>lay , (Tibark-Tbase) * tmerate, 0)
    hatch_gdd<-cumsum(rT)
    hatch<-tme[which(hatch_gdd>ddend['hatch'])][1]

    # L1 to L4 in inner bark - IMPOSE plateau to dev rate at 40C - see model doc
    if (!is.na(hatch)){
      Tbase<-11.9
      Tmx<-40
      rT<-ifelse(Tibark>Tbase & tme>hatch, (Tibark-Tbase) * tmerate, 0)
      rT<-ifelse(Tibark>Tmx & tme>hatch, (Tmx-Tbase) * tmerate, rT) # impose max rT
      larval_gdd<-cumsum(rT)
      pupate<-tme[which(larval_gdd>ddend['pupate'])][1]
    }
    # Score based on last phase complete 0,1 (larva),2(pupa)
    complete<-which(!is.na(c(eclosion,lay,hatch,pupate)))
    if (length(complete)==0) score<-0 else score <-ifelse(length(complete)<=1,0,max(complete)-2) #

    # Larval developmental completion 0 to 1 (pupates)of final phase
    dev_complete<-as.numeric(larval_gdd[length(larval_gdd)]/ddend['pupate'])
    dev_complete<-ifelse(dev_complete>=1,1.0,dev_complete)
  }
  # Format output
  output<-c(eclosion,lay,hatch,pupate)
  output<-as.numeric(strftime(output, format = "%j"))
  names(output)<-c("eclosion_doy","lay_doy","hatch_doy","pupate_doy")
  output<-c(output,"stages_complete"=score,"larval_development"=dev_complete)
  return(output)
}
