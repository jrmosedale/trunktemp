library(terra)
library(lubridate)
dir_in<-"/Users/jonathanmosedale/Data/pest_risks/alh/s0tiles"
dir_out<-"/Users/jonathanmosedale/Data/pest_risks/alh"
scenario_name<-"syc_s0"
template_file<-'/Users/jonathanmosedale/Library/CloudStorage/OneDrive-UniversityofExeter/jasmin/gb_1km.tif'
modelruns<-c("2011_2015","2012_2016","2013_2017","2014_2018","2015_2019","2016_2020")

t0<-now()
for(mrun in modelruns){
  r<-blend_all_tiles(dir_in,dir_out,template_file,scenario_name,mrun)
  plot(r)
  print(now()-t0)
}