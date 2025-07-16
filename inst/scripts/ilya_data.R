### Variables required
# Tree ID, Species, Sites, Lat, Lon, Veg height
# Datetime, PAI, Tair (survey tags), Quality_flag

# Date range - check
#AH NSP 17/02-19/10
#AH OAK 17/02-22/09
#BE ALL 14/02-01/10
#CW NSP 18/03-24/09
#CW OAK1&2 03/04-11/10
#CW OAK3  27/04-11/10
#KG ALL 01/03-18/09
#WH ALL 22/03-02/10

dir_out<-"data-raw"

## Get survey tag data directories for each site (excluding Kew)
dir_od<-"/Users/jonathanmosedale/Library/CloudStorage/OneDrive-UniversityofExeter"
dir_metoff<-file.path(dir_od,"MET Office Climate and Plant Biosecurity Project","Sites")
site_fullnames<-c("Alice Holt","Bellever","Cornwall","Wakehurst")
site_codes<-c("AH","BE","CW","WH")
dirs_stag<-c()
for(n in 1:length(site_fullnames)){
  dr<-file.path(dir_metoff,site_fullnames[n],paste(site_codes[n],"Data"),paste(site_codes[n],"SurveyTag Data"))
  print(dr)
  dirs_stag<-c(dirs_stag,dr)
}
dir.exists(dirs_stag)


### Survey tag data import
# Get all surveytag files - ignoring raw data dirs
library(stringr)
stfiles<-list.files(dirs_stag, pattern =".csv",recursive=TRUE,full.names=TRUE)
stfiles<-stfiles[which(!str_detect(stfiles,"Raw Data"))]
file.exists(stfiles)

# Extract metadata from file names and create id variable for matching files by tag/location
for(f in stfiles){
  metadata<-as.data.frame(get_surveytag_metadata(f))
  metadata<-cbind(filepath=f,metadata)
  if (f==stfiles[1]) file_metadata<-metadata else file_metadata<-rbind(file_metadata,metadata)
}

unique(file_metadata$id); unique(file_metadata$logger)
head(file_metadata); nrow(file_metadata)
View(table(file_metadata$id) )# =5 to 9 files per tree

# Number of different loggers per tree? BUt some loggers appear named the same?
loggernum_df<-file_metadata %>%
  group_by(id) %>%
  summarize(logger_num = n_distinct(logger))
max(loggernum_df$logger_num)

# Show loggers by id
file_metadata$tag_dir<-str_split_i(file_metadata$filepath,"/",12)
file_metadata$tag_dir_date<-as.Date(unlist(lapply(lapply(str_split(file_metadata$tag_dir," "),rev),"[",1)),format="%d%m%y")
View(file_metadata[order(file_metadata$id, file_metadata$tag_dir_date, decreasing = FALSE), ])


columns <- c("id","name1","name2","name3","name4","dir1","dir2","dir3","dir4")
idloggers_df<- data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(idloggers_df) <- columns

n<-0
for(id in unique(file_metadata$id)){
  n<-n+1
  print(paste(id,file_metadata$logger[which(file_metadata$id==id)]))
  logger_names<-unique(file_metadata$logger[which(file_metadata$id==id)])
  logger_dirs<-unique(file_metadata$tag_dir[which(file_metadata$id==id)])
  rw<-c(id,logger_names)
  idloggers_df[n,1:length(rw)]<-rw
  idloggers_df[n,6:(5+length(logger_dirs))] <-logger_dirs
}
View(idloggers_df)


### Load and merge all files by tree id while performing quality control checks
ids<-unique(file_metadata$id)

columns <- c("id","site","spp","tree","type","datetime", "air_temp","q_code")
stagdata_df<- data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(stagdata_df) <- columns

for(id in ids){ #  For each tree ID...
  alldat<- data.frame(matrix(nrow = 0, ncol = length(columns)))
  colnames(alldat) = columns

  #print(paste("Processing files for",id))
  iddat_df<-file_metadata[which(file_metadata$id==id),]

  # Check if data in file
  finfo<-file.info(iddat_df$filepath)
  ignore<-which(is.na(finfo$size) | finfo$size<500)
  if(length(ignore)>0){
    warning(paste("Ignoring file",iddat_df$filepath[ignore]))
    iddat_df<- iddat_df[-ignore,]
  }

  for(f in iddat_df$filepath){
    tagdat<-load_surveytag(f,rnd30=FALSE)
    tagdat$id<-id
    alldat<-rbind(alldat,tagdat)
  }
  alldat<-alldat[order(alldat$datetime),]

  # Remove readings of 0.00000 as should be NA
  zero_readings<-length(which(alldat$air_temp==0.00000))
  if(zero_readings>0) warning(paste("Removing",zero_readings,"zero readings"))
  alldat$air_temp<-ifelse(alldat$air_temp==0.00000,NA,alldat$air_temp)

  # Remove excessively high/low readings
  minlimit<- -5
  maxlimit<- 40
  extrmT<-which(alldat$air_temp < minlimit | alldat$air_temp > maxlimit )
  if(length(extrmT)>0) {
    warning(paste('Removing extreme reading of',alldat$air_temp[extrmT], 'recorded by logger',alldat$id[extrmT],',   ') )
    alldat<-alldat[-extrmT,]
  }

  # Remove duplicate timed readings - might be better as within xx mins or each other?
  dup<-duplicated(alldat$datetime)
  alldat<-alldat[!dup,]
  print(paste("After removing duplicate",nrow(alldat),"readings remain"))

  # Remove very isolated readings (no other within 1 hour)
  isolated_idx<-find_isolated_times(alldat$datetime,limit=3600)
  if(length(isolated_idx>0)){
    warning(paste('Removing',length(isolated_idx), 'isolated readings from' ,id))
    alldat<-alldat[-isolated_idx,]
  }
  # Join all location data in single df
  stagdata_df<-rbind(stagdata_df,alldat)
}

stagdata_df$q_code<-0
nrow(stagdata_df)
# Save results
saveRDS(stagdata_df,file.path(dir_out,'surveytag_air_logger_results.Rds'))


########## Quality issues ID from manual edit when compared with weather data etc
dir_out<-"data-raw"
stagdata_df<-readRDS(file.path(dir_out,'surveytag_air_logger_results.Rds'))
edits_file<-file.path("/Users/jonathanmosedale/Library/CloudStorage/OneDrive-UniversityofExeter/MetOff_data/field_data","Logger_results","manual_airtemp_edits.csv")
edits_df<-read.csv(edits_file)[,1:3]
edits_df$start<-as.POSIXct(trim(edits_df$start),"%d/%m/%Y %H:%M", tz=tz("GMT"))-(60*60)
edits_df$end<-as.POSIXct(edits_df$end,"%d/%m/%Y %H:%M", tz=tz("GMT"))+(60*60)

for (n in 1:nrow(edits_df)){
  id<-edits_df$id[n]
  sel<-which(stagdata_df$id==id & stagdata_df$datetime>=edits_df$start[n] & stagdata_df$datetime<=edits_df$end[n])
  print(paste(id,length(sel)))
  stagdata_df$q_code[sel]<-1
}


########### Merge with Tree & site data
stagdata_df$subsite<-toupper(substr(stagdata_df$site,1,1))
stagdata_df$subsite<-ifelse(stagdata_df$spp=="nsp",paste0(stagdata_df$subsite,"1"),paste0(stagdata_df$subsite,"2"))

# Site props head(sites_sf)
stagdata_df<-merge(stagdata_df,as.data.frame(sites_sf)[,c("subsite","elev","aspect","slope","soil")],by="subsite",all.x=TRUE)

# Tree lat lon
trees_sf$lat<-st_coordinates(trees_sf)[,2]
trees_sf$lon<-st_coordinates(trees_sf)[,1]

stagdata_df<-merge(stagdata_df,as.data.frame(trees_sf)[,c("id","lat","lon","Tree_height")],by="id",all.x=TRUE)
names(stagdata_df)<-c(names(stagdata_df)[1:length(stagdata_df)-1],"height")


########### Merge with LAI data - enforce min lai for Norway spruce as >= min observed value
lai_models<-readRDS(file.path(dir_lai,'lc_lai_models.Rds'))
# Fitted models file

stagdata_df$lai<-NA
for(t in names(lai_models)){
  model<-lai_models[[t]]
  sel<-which(stagdata_df$id==t)
  # Get day of year from result datetime
  minlai<-min(model$subset$lai)
  maxlai<-max(model$subset$lai)
  doy<-lubridate::yday(stagdata_df$datetime[sel])
  pred_lai<-predict(model, data.frame(date=doy))
  if(substr(t,4,6)=="nsp") pred_lai<-ifelse(pred_lai<minlai,minlai,pred_lai)
  stagdata_df$lai[sel]<-round(pred_lai,2)
}

for(id in unique(stagdata_df$id)){
  sel<-which(stagdata_df$id==id)
  plot(stagdata_df$lai[sel]~lubridate::yday(stagdata_df$datetime[sel]),main=id,ylim=c(0,4.2),xlim=c(80,280))
}


########## Save
dir_out<-"/Users/jonathanmosedale/Library/CloudStorage/OneDrive-UniversityofExeter/Shared folders/treetemps"
write.csv(stagdata_df,file.path(dir_out,'air_logger_data.csv'),row.names=FALSE)

########## timeseries graphs
# Air temp data
stagdata_df<-read.csv(file.path(dir_out,'air_logger_data.csv'))
names(stagdata_df)
# Add time to midnight values
unique(nchar(stagdata_df$datetime))
sel<-which(nchar(stagdata_df$datetime)==10)
stagdata_df$datetime[sel]<-paste0(stagdata_df$datetime[sel]," 00:00:00")
stagdata_df$datetime<-as.POSIXct(stagdata_df$datetime,tz="GMT",format="%Y-%m-%d %H:%M:%OS")


dir_weather<-file.path(dir_od,"MetOff_data","field_data","weather")

ids<-unique(stagdata_df$id)
id<-ids[18]

# Weather data
subsite<-trees_sf$subsite[which(trees_sf$id==id)]
weather<-read.csv(file.path(dir_weather,paste0('weather_',subsite,'_2024.csv')))
weather$obs_time<-as.POSIXct(weather$obs_time,tz="UTC")

sel<-which(stagdata_df$id==id & stagdata_df$q_code==0)
tair<-stagdata_df$air_temp[sel]
tme<-stagdata_df$datetime[sel]
wtemp<-match_temperatures(temp_in=weather$temp,dt_in=weather$obs_time,dt_out=tme,tme_thold=60*60)

figdata<-data.frame(
  "id"=id,
  "datetime"=tme,
  "weather_temp"=wtemp,
  "observed_air"=tair)

tseries<- xts(x = figdata, order.by = tme)
dygraph(tseries, main=id,group='group1') %>% dyAxis("y", label = "Temp (C)") %>%
  dyRangeSelector(height=20) %>%
  dySeries("weather_temp",color = 'black') %>%
  dySeries("observed_air",color = 'green') %>%
  dyLegend(width=600)










# Original version including Kew
df<-readRDS(file.path("/Users/jonathanmosedale/Library/CloudStorage/OneDrive-UniversityofExeter/MetOff_data/field_data/Logger_results",'surveytag_air_logger_results.Rds'))
unique(df$site)
nrow(df[-which(df$site=="kg"),])
# Edited version
df<-readRDS(file.path("/Users/jonathanmosedale/Library/CloudStorage/OneDrive-UniversityofExeter/MetOff_data/field_data/Logger_results",'logger_results_stedit.Rds'))
unique(df$site)
nrow(df[-which(df$site=="kg"|is.na(df$air_temp)),])








