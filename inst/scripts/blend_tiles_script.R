library(terra)
library(lubridate)
source('/Users/jonathanmosedale/Library/CloudStorage/OneDrive-UniversityofExeter/Rprojects/trunktemp/R/blend_raster_functions.R')


############### Merge Asian Longhorn Risk values ########################
dir_out<-"/Users/jonathanmosedale/Data/pest_risks/alh/riskmaps"
dir_inroot<-"/Users/jonathanmosedale/Data/pest_risks/alh"
modelruns<-c("2011_2015","2012_2016","2013_2017","2014_2018","2015_2019","2016_2020")

scenarios<-c(0,1,2,3,4)
#scenario_name<-"syc_s2"
template_file<-'/Users/jonathanmosedale/Library/CloudStorage/OneDrive-UniversityofExeter/jasmin/gb_1km.tif'
file.exists(template_file)

#year_range<-modelruns[2]
for(s in scenarios){
  scenario_name<-paste0("syc_s",s)
  dir_in<-file.path(dir_inroot,paste0("s",s,"tiles"))
  for(year_range in modelruns){
    print(substring(year_range,1,4))
    ukresults.r<-merge_all_tiles(dir_in,scenario_name,year_range,template_file)
    #ukresults.r<-blend_all_tiles(dir_in,template_file,scenario_name,year_range)
    plot(ukresults.r[[1]],main=paste(s,year_range))
    writeRaster(ukresults.r,file.path(dir_out,paste0("alh_",scenario_name,"1kmrisks_5yr_",year_range,".tif")),overwrite=TRUE)
  }
}




