
dir_pests<-"/Users/jonathanmosedale/Data/pest_risks"
dir_outputs<-"/Users/jonathanmosedale/Library/CloudStorage/OneDrive-UniversityofExeter/MetOff_data/pest_raster_share"
ipfiles<-list.files(file.path(dir_pests,"ips"),pattern="*.tif")

r<-rast(file.path(dir_pests,"ips",ipfiles[37]))
plot(r,main='')

# Check ip
dir_ips<-file.path(dir_pests,'ips')
dir_dm<-file.path(dir_pests,'dm')
dir_clim<-file.path(dir_pests,'clim')

for (s in 1:5){
  for (year in 2011:2020){
    # Ips
    f<-file.path(dir_ips,paste0("ip_tt_nsp_s",s,"_",year,".tif"))
    if(file.exists(f)==FALSE) print(paste(f,"does NOT exist!!!")) else{ 
      r<-rast(f)[["generations_complete"]]
      plot(r,main=paste(s,year)) }
  }
}
# ips have NA values - reset to 0???


for (s in 1:5){
  for (year in 2011:2020){
    # DM
    f<-file.path(dir_dm,paste0("dm_tt_nsp_s",s,"_",year,".tif"))
    if(file.exists(f)==FALSE) print(paste(f,"does NOT exist!!!")) else{ 
      r<-rast(f)[["dnd_stages_complete"]]
      plot(r,main=paste(s,year)) }
  }
}


for (s in 1:5){
  for (year in 2011:2020){
    # DM
    f<-file.path(dir_clim,paste0("climate_nsp_s",s,"_",year,".tif"))
    if(file.exists(f)==FALSE) print(paste(f,"does NOT exist!!!")) else{ 
      r<-rast(f)[["trunk_gdd10"]]
      plot(r,main=paste(s,year)) }
  }
}

############ Write Ips and D micans risk rasters ################
template_file<-'/Users/jonathanmosedale/Library/CloudStorage/OneDrive-UniversityofExeter/jasmin/gb_1km.tif'
gb1km<-rast(template_file)
plot(gb1km)

for (s in 1:5){
  for (year in 2011:2020){
    # Dm
    f<-file.path(dir_dm,paste0("dm_tt_nsp_s",s,"_",year,".tif"))
    r<-rast(f)[["dnd_stages_complete"]]
    r<-mask(crop(ifel(is.na(r),0,r),gb1km),gb1km)
    writeRaster(r,file.path(dir_outputs,paste0("Dmican_s",s,"_",year,"_winterstage.tif")),overwrite=TRUE)
    # Ips
    f<-file.path(dir_ips,paste0("ip_tt_nsp_s",s,"_",year,".tif"))
    r<-rast(f)[["generations_complete"]]
    r<-mask(crop(ifel(is.na(r),0,r),gb1km),gb1km)
    writeRaster(r,file.path(dir_outputs,paste0("Itypog_s",s,"_",year,"_gencomplete.tif")),overwrite=TRUE)
    
  }
}
