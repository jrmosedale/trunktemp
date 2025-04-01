################################################################
# Files to load logger data and get file info from filenames
################################################################
load_surveytag<-function(filepath,addmeta=TRUE, rnd30=TRUE){
  if(file.exists(filepath)==FALSE) stop('File NOT found!!!')
  dat<-read.csv(filepath,skip=3)
  names(dat)<-c('Year', 'Month', 'Day', 'Hour', 'Minute', 'Second', 'air_temp', 'CJ_Reading', 'vBatt', 'Samples','Fault_Code')
  dat$DTstring<-paste0(formatC(dat$Day, width=2, flag="0"),"/",formatC(dat$Month, width=2, flag="0"),"/",  dat$Year," ",formatC(dat$Hour, width=2, flag="0"),":", formatC(dat$Minute, width=2, flag="0"))
  dat$datetime<-as.POSIXct(strptime(dat$DTstring,"%d/%m/%Y %H:%M"), tz='GMT')
  head(dat)
  # Add metadata to df??
  if(addmeta){
    metadata<-as.data.frame(get_surveytag_metadata(filepath))
    metadata<-metadata[rep(seq_len(nrow(metadata)), each = nrow(dat)), ]
    dat<-cbind(metadata,dat)
  }
  # Round times to nearest 30mins??
  if(rnd30){
    dat$datetime<-lubridate::round_date(dat$datetime,"30 mins")
  }
  dat<-dat[,c('id','site','spp','tree','type','datetime','air_temp')]
  return(dat)
}


get_surveytag_metadata<-function(filepath){
  if(file.exists(filepath)==FALSE) stop('File NOT found!!!')
  filename<-basename(filepath)
  site<-tolower(str_split_i(filename,pattern="[_ ]+",i=1))
  spp<-tolower(str_split_i(filename,pattern="[_ ]+",i=2))
  tree<-tolower(str_split_i(filename,pattern="[_ ]+",i=3))
  #logger<-str_split_i(filename,pattern="[_ ]+",i=3))
  type<-'air'
  id<-paste(site,spp,tree,sep='_')
  # Get logger code from within file
  dat<-readLines(filepath,n=2)[2]
  logger<-str_split_i(dat,pattern=",",i=1)
  return(list(id=id,site=site,spp=spp,tree=tree,type=type, logger=logger))
}

get_tinytag_metadata<-function(filename){
  filename<-basename(filename)
  site<-tolower(substr(filename,1,2))
  spp<-tolower(substr(filename,4,6))
  tree<-tolower(substr(filename,8,8))
  asp<-tolower(substr(filename,10,10))
  type<-'trunk'
  id<-paste(site,spp,tree,sep='_')
  fullid<-tolower(paste(site,spp,tree,asp,sep='_'))
  return(list(fullid=fullid,id=id,site=site,spp=spp,tree=tree,asp=asp,type=type))
}

# Uses UTF8 - no unit .csv output from TinyTag software save table
load_tinytag<-function(filepath,addmeta=TRUE, rnd30=TRUE){
  if(file.exists(filepath)==FALSE) stop('File NOT found!!!')
  dat<-read.csv(filepath,skip=4)

  # ASSUMES black probe = outer, white probe = inner
  names(dat)<-c('Reading', 'DTstring', 'inner_temp', 'outer_temp')

  # Convert date to datetime
  #dat$datetime<-as.POSIXct(strptime(dat$DTstring,"%Y-%m-%d %H:%M:%S"))
  dat$datetime<-as.POSIXct(strptime(dat$DTstring,"%Y-%m-%d %H:%M:%S"), tz='GMT')
  head(dat)
  # Add metadata?
  if(addmeta){
    metadata<-as.data.frame(get_tinytag_metadata(basename(filepath)))
    metadata<-metadata[rep(seq_len(nrow(metadata)), each = nrow(dat)), ]
    dat<-cbind(metadata,dat)
  }
    # Round times to nearest 30mins??
  if(rnd30){
    dat$datetime<-lubridate::round_date(dat$datetime,"30 mins")
  }
  # Round times to nearest 30mins
  dat<-dat[,c('fullid','id','site','spp','tree','asp','type','datetime','inner_temp','outer_temp')]
  head(dat)
  return(dat)
}

################################################################
# Filter and merge logger data
################################################################
# Create regular time series to match start and end times of existing datetime
regular_timeseries<-function(dt,interval="30 mins"){
  # Get start and end
  dt_start<-min(dt)
  dt_end<-max(dt)
  # Round to nearest interval
  dt_start<-lubridate::round_date(dt_start,"30 mins")
  dt_end<-lubridate::round_date(dt_end,"30 mins")
  # Create series
  dt_series<-seq(from=dt_start,by=interval,to=dt_end)
  return(dt_series)
}

# Functions for linear interpolating logger values to new datetimes
# temp_in = Temperature values at date times dt_in
# dt_in datetimes of temp_in
# dt_out = required datetimes
# tme_thold = time (secs) within which values are interpolated/assigned
# Returns interpoloated air temperature values to match dt
# If no logger value within tme_thold - returns NA
# If only one value within tme_thold - assigns that value to dt
match_temperatures <- function(temp_in, dt_in, dt_out,tme_thold){
  # Process inputs
  if(length(temp_in)!=length(dt_in)) stop('Input temperature and datetimes to function interp_temp  not of same length!!')
  df<-data.frame("datetime"=as.POSIXct(dt_in,tz='GMT'),"temp"=temp_in)
  df_out<-data.frame(datetime=dt_out)

  # Match logger temp data for dt_out where datetime correspond
  df_out<-left_join(df_out,df,by="datetime")
  nrow(df_out)

  # Where no perfect match interpolate temperatures for dt_out
  missing<-which(is.na(df_out$temp))
  missing_vals<-unlist(sapply(df_out[missing,"datetime"],interp_temperatures,df=df, tme_thold))
  df_out$temp[missing]<-missing_vals
  if(anyNA(df_out$temp)) warning("Failed interpolation of some temperatures in match_temperatures!!")
  return(df_out$temp)
}

interp_temperatures<-function(missing_dt,df,tme_thold){
  # Find if valid temp reading before or after missing datetime
  nxt<-which(as.numeric(df$datetime)-as.numeric(missing_dt)>0 & as.numeric(df$datetime)-as.numeric(missing_dt)<tme_thold)[1]
  prv<-tail(which(as.numeric(missing_dt)-as.numeric(df$datetime)>0 & as.numeric(missing_dt)-as.numeric(df$datetime)<tme_thold))[1]
  # Interpolate / assign or NA according to whether data available
  if(length(prv>0) & length(nxt)>0){
    y0<-df$temp[prv]
    y1<-df$temp[nxt]
    x0<-as.numeric(df$datetime[prv])
    x1<-as.numeric(df$datetime[nxt])
    x<-as.numeric(missing_dt)
    t<-( y0*(x1-x)+y1*(x-x0) ) / ( x1-x0 )
  }
  if(length(prv)>0 & length(nxt)==0) t<-df$temp[prv]
  if(length(prv)==0 & length(nxt)>0) t<-df$temp[nxt]
  if(length(prv)==0 & length(nxt)==0) t<-NA
  return(t)
}

# Filter cols and combine by datetime
filter_data<-function(df){
  df<-subset(df,select=c(datetime,TC_Reading,CJ_Reading))
  #bktemp_df<-subset(bktemp_df,select=c(datetime,bark_temp))
  return(df)
}

# Merge by datetime
merge_data<-function(df1,df2){
  range(df1$datetime)
  range(df2$datetime)
  results<-merge(df1,df2,by=c("datetime"))
  results$day_of_yr<-lubridate::yday(results$datetime)
  head(results); tail(results)
  #TC_reading = air temp
  names(results)<-c("datetime","air_temp", "CJ_Reading", "bark_temp","day_of_yr")
  return(results)
}

# Identify isolated times in timeseries where no other time within limit in seconds - returns index of isolated
# see https://stackoverflow.com/questions/54862344/remove-isolated-elements-of-a-vector
find_isolated_times<-function(x,limit=1800){
  x[order(x)]
  #x[c(FALSE, diff(as.numeric(x)) <= limit) | c(diff(as.numeric(x)) <= limit, FALSE)] # returns non-isolated values
  #sel<-which(c(FALSE, diff(as.numeric(x)) <= limit) | c(diff(as.numeric(x)) <= limit, FALSE) ) # returns non-isolated idx
  sel<-which(c(TRUE, diff(as.numeric(x)) > limit) & c(diff(as.numeric(x))> limit, TRUE) ) # returns idx of isolated
  return(sel)
}


################################################################
# Plot data
################################################################
# Plots all variables ending in "_temp" vs "datetime" variable of a dataframe
plot_tseries<-function(dat,title='Timseries plot',legloc='bottomright',cols=c('red','orange','green','blue','cyan'),
                       ylimits=NA, xlimits=NA){
    varnames<-names(dat)
    tempvars<-varnames[str_sub(varnames, start=-4)=='temp']
    if(length(ylimits)!=2){
      miny<-min(dat[,tempvars],na.rm=TRUE)
      maxy<-max(dat[,tempvars],na.rm=TRUE)
    }
    if(length(ylimits)==2){
      miny<-ylimits[1]
      maxy<-ylimits[2]
    }
    tme<-as.POSIXct(dat$datetime)
    if(length(xlimits)!=2){
      minx<-min(tme,na.rm=TRUE)
      maxx<-max(tme,na.rm=TRUE)
    }
    if(length(xlimits)==2){
      minx<-xlimits[1]
      maxx<-xlimits[2]
    }

    # Plot
    plot(dat[,tempvars[1]]~tme,
         type = "l", cex.axis=0.85,
         col = cols[1],
         ylim = c(miny, maxy),xlim = c(minx,maxx),
         xlab = "Date",
         ylab = 'Temperature C',
         xaxt="n",main=title)
    # Add date time axis
    tstep<-as.numeric(dat$datetime[2]-dat$datetime[1])
    if(tstep!=1){
      day_range <- as.POSIXct(round(range(dat$datetime), "day")) # change if >1 month??
      axis.POSIXct(1, at = seq(day_range[1], day_range[2], by = "day"), format = "%d/%m")
    }
    # Add other temp variables
    if(length(tempvars)>1){
      for(n in 2:length(tempvars)) lines(dat[,tempvars[n]]~tme, type = "l", col = cols[n])
    }
    legend(legloc,tempvars,lty = 1, col=cols,cex=0.75)
    return()
  }


ggplot_tseries<-function(dat){
  ggplot(data=dat, aes(x=datetime, y=air_temp, group = id,
                                              colour = id)) +
    geom_line() +labs(y= "Air temperature loggers", x = "Date")
  #covid_plot + ggtitle("Daily Deaths for European countries in March,2020")+geom_point()

}

# DAILY daily min max
plot_daily<-function(dat){

  # Get daily moin max
  dat$date <- as.Date(dat$datetime)
  air_min<-aggregate(dat$air_temp, by=dat["date"], min)
  air_max<-aggregate(dat$air_temp, by=dat["date"], max)
  bark_min<-aggregate(dat$bark_temp, by=dat["date"], min)
  bark_max<-aggregate(dat$bark_temp, by=dat["date"], max)

  plotdat<-cbind(air_min,air_max$x,bark_min$x,bark_max$x)
  names(plotdat)<-c("date","min_air_temp","max_air_temp","min_bark_temp","max_bark_temp")
  min_t<-min(plotdat[,c(2:5)])
  max_t<-max(plotdat[,c(2:5)])

  # PLlt
  plot(plotdat$date,                              # Draw first time series
       plotdat$min_air_temp,
       type = "l",
       col = 4,
       ylim = c(min_t, max_t),
       xlab = "Date",
       ylab = "Temperature",
       main=paste(a,": air and under-bark daily min/max"),
       xaxt="n")

  lines(plotdat$date,                             # Draw second time series
        plotdat$max_air_temp,
        type = "l",
        col = 3)
  # Bark T
  lines(plotdat$date,                             # Draw third time series
        plotdat$min_bark_temp,
        type = "l",
        col = 1)

  lines(plotdat$date,                             # Draw 4th time series
        plotdat$max_bark_temp,
        type = "l",
        col = 2)

  legend("topright",                           # Add legend to plot
         c("Min air temp", "Max air temp","Min bark temp (paddle)", "Max bark temp (paddle)"),
         lty = 1,
         col = c(4,3,1,2),
         cex = 0.75 # Change legend size
  )
}


