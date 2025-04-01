library(terra)
library(sf)
library(ggplot2)
library(hrbrthemes)
library(tidyterra)

############### Plot Asian Longhorn Risk values ########################
dir_alh<-"/Users/jonathanmosedale/Data/pest_risks/alh/riskmaps"
template_file<-'/Users/jonathanmosedale/Library/CloudStorage/OneDrive-UniversityofExeter/jasmin/gb_1km.tif'

modelruns<-c("2011_2015","2012_2016","2013_2017","2014_2018","2015_2019","2016_2020")
scenarios<-c(0,1,2,3,4)


coolest<-"2011_2015"
warmest<-"2016_2020"

s0warm.r<-rast(file.path(dir_alh,paste0("alh_syc_","s0","1kmrisks_5yr_",warmest,".tif")))
s0cool.r<-rast(file.path(dir_alh,paste0("alh_syc_","s0","1kmrisks_5yr_",coolest,".tif")))
s4warm.r<-rast(file.path(dir_alh,paste0("alh_syc_","s4","1kmrisks_5yr_",warmest,".tif")))
s4cool.r<-rast(file.path(dir_alh,paste0("alh_syc_","s4","1kmrisks_5yr_",coolest,".tif")))

# Proportion completing lifecyle - classes
lyr<-"pop_completed"
r<-c(s0warm.r[[lyr]],s0cool.r[[lyr]],s4warm.r[[lyr]],s4cool.r[[lyr]])
lyr_names<-c(paste("S0","2016-2020"),
            paste("S0","2011-2015"),
            paste("S4","2016-2020"),
            paste("S4","2011-2015"))
rg<-c(min(global(r,min,na.rm=T)),max(global(r,max,na.rm=T)))


# Plot stage completed as class - convert to % and round up!
r<-ceiling(r*100)
vals<-as.numeric(unlist(unique(r[[1]])))
levels(r[[1]]) <- data.frame(id=vals, Completed=c("None","<1%", "8%", "28%","54%","69%"))
levels(r[[2]]) <- data.frame(id=vals, Completed=c("None","<1%", "8%", "28%","54%","69%"))
levels(r[[3]]) <- data.frame(id=vals, Completed=c("None","<1%", "8%", "28%","54%","69%"))
levels(r[[4]]) <- data.frame(id=vals, Completed=c("None","<1%", "8%", "28%","54%","69%"))
is.factor(r)

#panel(r, col = map.pal("viridis",6),axes=FALSE,main=lyr_names)
#panel(r, col = map.pal("plasma",6),axes=FALSE,main=lyr_names)
# Tailored - reduces space between plot!
panel(r, nc=4,nr=1,col = c("#440154","#6A00A8", "#B12A90", "#E16462", "#FCA636", "#F0F921"),
      plg=list(title.cex=2,size=2,x="top"),axes=FALSE,main=lyr_names)


########### NUmber of days to complete life cycle -
lyr<-"days_to_complete_1"
r<-c(s0warm.r[[lyr]],s0cool.r[[lyr]],s4warm.r[[lyr]],s4cool.r[[lyr]])
lyr_names<-c(paste("S0","2016-2020"),
             paste("S0","2011-2015"),
             paste("S4","2016-2020"),
             paste("S4","2011-2015"))
rg<-c(min(global(r,min,na.rm=T)),max(global(r,max,na.rm=T)))


lsea.r<-rast(template_file)
lsea.r<-ifel(is.na(lsea.r),NA,0)
lsea.v<-as.polygons(lsea.r)

panel(r, nc=4,nr=1,col = rev(map.pal("viridis",8)) ,
      type="continuous",range=c(600,1600),breaks<-c(600,800,1000,1200,1400,1600),
     fun=function()lines(lsea.v),
      plg=list(title.cex=1),
      axes=FALSE,
     main="")

