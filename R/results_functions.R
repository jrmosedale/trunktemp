#Adds solar zenith to timeseries database
# df timeseries df with datetime as variable and idvar shared with sites_sf
# sites_sf = id and lat lon data
add_solar_zenith<-function(df,locs_sf){
  # Calc lon/lat of locations/trees
  lat<-st_coordinates(locs_sf)[,2]
  lon<-st_coordinates(locs_sf)[,1]
  locs_df<-data.frame(id=locs_sf$id,lat=lat,lon=lon)
  df<-merge(df,locs_df,by="id",all.x=TRUE)
  solarpos<-micropoint::solarposition(df$lat,df$lon,as.POSIXlt(df$datetime),NA,NA,NA,merid=0,dst=0)
  df$szenith<-solarpos$zen
  # By solar angle group
  minsa<-min(df$szenith)
  df$szgrp<-as.numeric(cut_number(df$szenith,12))
  return(df)
}

# Calculate phase shift in daytime max temperature
calculate_phase_shift<-function(hrlydata,mxfrq,minhr=6,maxhr=20 ){
  hrlydata$date<-as.Date(hrlydata$datetime)
  hrlydata$hour<-hour(hrlydata$datetime)

  ## Calculate phase shift from weather/air to trunk temperatures - daytime only
  sel<-which(hrlydata$hour %in% seq(minhr,maxhr,1))
  dailydata<-hrlydata[sel,] %>%
    group_by(id,date) %>%
    summarise(wtrtemp_max = max(weather_temp), wtrtemp_hr = which.max(weather_temp)/2,
              airtemp_max= max(air_temp), airtemp_hr = which.max(air_temp)/2,
              s_outtemp_max = max(s_outer_temp), s_outtemp_hr = which.max(s_outer_temp)/2,
              s_inrtemp_max = max(s_inner_temp), s_inrtemp_hr = which.max(s_inner_temp)/2,
              n_outtemp_max = max(n_outer_temp), n_outtemp_hr = which.max(n_outer_temp)/2,
              n_inrtemp_max = max(n_inner_temp), n_inrtemp_hr = which.max(n_inner_temp)/2
    )

  phseshift<- dailydata %>%
    group_by(id,date) %>%
    summarise(air_wtr_anom = airtemp_hr-wtrtemp_hr ,
              s_out_anom = s_outtemp_hr-wtrtemp_hr, s_inr_anom = s_inrtemp_hr-wtrtemp_hr,
              n_out_anom = n_outtemp_hr-wtrtemp_hr, n_inr_anom = n_inrtemp_hr-wtrtemp_hr)

  # Plot histograms of phase shift
  par(mfrow=c(3,2))
  hist(phseshift$air_wtr_anom,main="Air temp",ylim=c(0,mxfrq),ylab="",xlab="",xlim=c(-10,10), breaks=seq(-28,28,0.5))
  mn<-mean(phseshift$air_wtr_anom); sd<-sd(phseshift$air_wtr_anom)
  abline(v=mn,col='red',lwd=1.5); #abline(v=mn+sd,col='blue',lwd=1.5); abline(v=mn-sd,col='blue',lwd=1.5)

  hist(phseshift$air_wtr_anom,main="Air temp",ylim=c(0,mxfrq),ylab="",xlab="",xlim=c(-10,10), breaks=seq(-28,28,0.5))
  mn<-mean(phseshift$air_wtr_anom); sd<-sd(phseshift$air_wtr_anom)
  mn<-mean(phseshift$air_wtr_anom); sd<-sd(phseshift$air_wtr_anom)
  abline(v=mn,col='red',lwd=1.5); #abline(v=mn+sd,col='blue',lwd=1.5); abline(v=mn-sd,col='blue',lwd=1.5)

  hist(phseshift$s_out_anom,main="S outer",ylim=c(0,mxfrq),ylab="",xlab="",xlim=c(-10,10), breaks=seq(-28,28,0.5))
  mn<-mean(phseshift$air_wtr_anom); sd<-sd(phseshift$air_wtr_anom)
  mn<-mean(phseshift$s_out_anom); sd<-sd(phseshift$s_out_anom)
  abline(v=mn,col='red',lwd=1.5); #abline(v=mn+sd,col='blue',lwd=1.5); abline(v=mn-sd,col='blue',lwd=1.5)

  hist(phseshift$n_out_anom,main="N outer",ylim=c(0,mxfrq),ylab="",xlim=c(-10,10), breaks=seq(-28,28,0.5))
  mn<-mean(phseshift$air_wtr_anom); sd<-sd(phseshift$air_wtr_anom)
  mn<-mean(phseshift$n_out_anom); sd<-sd(phseshift$n_out_anom)
  abline(v=mn,col='red',lwd=1.5); #abline(v=mn+sd,col='blue',lwd=1.5); abline(v=mn-sd,col='blue',lwd=1.5)

  hist(phseshift$s_inr_anom,main="S inner",xlab="Phase shift in hrs",ylim=c(0,mxfrq),ylab="",xlim=c(-10,10), breaks=seq(-28,28,0.5))
  mn<-mean(phseshift$s_inr_anom); sd<-sd(phseshift$s_inr_anom)
  abline(v=mn,col='red',lwd=1.5); #abline(v=mn+sd,col='blue',lwd=1.5); abline(v=mn-sd,col='blue',lwd=1.5)

  hist(phseshift$n_inr_anom,main="N inner",xlab="Phase shift in hrs",ylim=c(0,mxfrq),ylab="",xlim=c(-10,10), breaks=seq(-28,28,0.5))
  mn<-mean(phseshift$n_inr_anom); sd<-sd(phseshift$n_inr_anom)
  abline(v=mn,col='red',lwd=1.5); #abline(v=mn+sd,col='blue',lwd=1.5); abline(v=mn-sd,col='blue',lwd=1.5)

  # Print mean and sd values
  results_df<-data.frame(Temperature=c("Air","S_outer","S_inner","N_outer","N_inner"),
                         Mean=c(mean(phseshift$air_wtr_anom),mean(phseshift$s_out_anom),mean(phseshift$s_inr_anom),mean(phseshift$n_out_anom),mean(phseshift$n_inr_anom)),
                         sd=c(sd(phseshift$air_wtr_anom),sd(phseshift$s_out_anom),sd(phseshift$s_inr_anom),sd(phseshift$n_out_anom),sd(phseshift$n_inr_anom)) )
  return(results_df)
}

calculate_anomalies<-function(hrlydata,mxfrq,stat="max",minhr=0,maxhr=24){
  hrlydata$date<-as.Date(hrlydata$datetime)
  hrlydata$hour<-hour(hrlydata$datetime)

  sel<-which(hrlydata$hour %in% seq(minhr,maxhr,1))

  ## Calculate daily min/max anomaly
  if(stat=="max"){
    dailydata<-hrlydata[sel,] %>%
      group_by(id,date) %>%
      summarise(wtrtemp_max = max(weather_temp), wtrtemp_hr = which.max(weather_temp)/2,
                airtemp_max= max(air_temp), airtemp_hr = which.max(air_temp)/2,
                s_outtemp_max = max(s_outer_temp), s_outtemp_hr = which.max(s_outer_temp)/2,
                s_inrtemp_max = max(s_inner_temp), s_inrtemp_hr = which.max(s_inner_temp)/2,
                n_outtemp_max = max(n_outer_temp), n_outtemp_hr = which.max(n_outer_temp)/2,
                n_inrtemp_max = max(n_inner_temp), n_inrtemp_hr = which.max(n_inner_temp)/2
      )
  }
  if(stat=="min"){
    dailydata<-hrlydata[sel,] %>%
      group_by(id,date) %>%
      summarise(wtrtemp_max = min(weather_temp), wtrtemp_hr = which.min(weather_temp)/2,
                airtemp_max= min(air_temp), airtemp_hr = which.min(air_temp)/2,
                s_outtemp_max = min(s_outer_temp), s_outtemp_hr = which.min(s_outer_temp)/2,
                s_inrtemp_max = min(s_inner_temp), s_inrtemp_hr = which.min(s_inner_temp)/2,
                n_outtemp_max = min(n_outer_temp), n_outtemp_hr = which.min(n_outer_temp)/2,
                n_inrtemp_max = min(n_inner_temp), n_inrtemp_hr = which.min(n_inner_temp)/2
      )
  }
  anomaly<- dailydata %>%
    group_by(id,date) %>%
    summarise(air_wtr_anom = airtemp_max-wtrtemp_max ,
              s_out_anom = s_outtemp_max-wtrtemp_max, s_inr_anom = s_inrtemp_max-wtrtemp_max,
              n_out_anom = n_outtemp_max-wtrtemp_max, n_inr_anom = n_inrtemp_max-wtrtemp_max)
  # Plot histograms of anomaly
  par(mfrow=c(3,2))
  hist(anomaly$air_wtr_anom,main="Air temp",ylim=c(0,mxfrq),ylab= "",xlab="",xlim=c(-10,10),breaks=seq(-20,20,1))
  mn<-mean(anomaly$air_wtr_anom,na.rm=T); sd<-sd(anomaly$air_wtr_anom,na.rm=T)
  abline(v=mn,col='red',lwd=1.5); #abline(v=mn+sd,col='blue',lwd=1.5); abline(v=mn-sd,col='blue',lwd=1.5)

  hist(anomaly$air_wtr_anom,main="Air temp",ylim=c(0,mxfrq),ylab="",xlab="",xlim=c(-10,10),breaks=seq(-20,20,1))
  mn<-mean(anomaly$air_wtr_anom,na.rm=T); sd<-sd(anomaly$air_wtr_anom,na.rm=T)
  abline(v=mn,col='red',lwd=1.5); #abline(v=mn+sd,col='blue',lwd=1.5); abline(v=mn-sd,col='blue',lwd=1.5)

  hist(anomaly$s_out_anom,main="S outer",ylim=c(0,mxfrq),ylab="",xlab="",xlim=c(-10,10),breaks=seq(-20,40,1))
  mn<-mean(anomaly$s_out_anom,na.rm=T); sd<-sd(anomaly$s_out_anom,na.rm=T)
  abline(v=mn,col='red',lwd=1.5); #abline(v=mn+sd,col='blue',lwd=1.5); abline(v=mn-sd,col='blue',lwd=1.5)

  hist(anomaly$n_out_anom,main="N outer",ylim=c(0,mxfrq),ylab="",xlab="",xlim=c(-10,10),breaks=seq(-20,20,1))
  mn<-mean(anomaly$n_out_anom,na.rm=T); sd<-sd(anomaly$n_out_anom,na.rm=T)
  abline(v=mn,col='red',lwd=1.5); #abline(v=mn+sd,col='blue',lwd=1.5); abline(v=mn-sd,col='blue',lwd=1.5)

  hist(anomaly$s_inr_anom,main="S inner",xlab="Temperature anomaly in deg C",ylim=c(0,mxfrq),ylab="",xlim=c(-10,10),breaks=seq(-20,20,1))
  mn<-mean(anomaly$s_inr_anom,na.rm=T); sd<-sd(anomaly$s_inr_anom,na.rm=T)
  abline(v=mn,col='red',lwd=1.5); #abline(v=mn+sd,col='blue',lwd=1.5); abline(v=mn-sd,col='blue',lwd=1.5)

  hist(anomaly$n_inr_anom,main="N inner",xlab="Temperature anomaly in deg C",ylim=c(0,mxfrq),ylab="",xlim=c(-10,10),breaks=seq(-20,20,1))
  mn<-mean(anomaly$n_inr_anom,na.rm=T); sd<-sd(anomaly$n_inr_anom,na.rm=T)
  abline(v=mn,col='red',lwd=1.5); #abline(v=mn+sd,col='blue',lwd=1.5); abline(v=mn-sd,col='blue',lwd=1.5)

  results_df<-data.frame(Temperature=c("Air","S_outer","S_inner","N_outer","N_inner"),
                         Mean=c(mean(anomaly$air_wtr_anom,na.rm=T),mean(anomaly$s_out_anom,na.rm=T),mean(anomaly$s_inr_anom,na.rm=T),mean(anomaly$n_out_anom,na.rm=T),mean(anomaly$n_inr_anom,na.rm=T)),
                         sd=c(sd(anomaly$air_wtr_anom,na.rm=T),sd(anomaly$s_out_anom,na.rm=T),sd(anomaly$s_inr_anom,na.rm=T),sd(anomaly$n_out_anom,na.rm=T),sd(anomaly$n_inr_anom,na.rm=T)) )
  return(results_df)

}
