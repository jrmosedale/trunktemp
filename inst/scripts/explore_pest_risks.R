# Explore pest risk outputs
library(terra)
library(sf)
library(ggplot2)
library(hrbrthemes)
library(tidyterra)

dir_clim<-"/Users/jonathanmosedale/Data/pest_risks/clim"
dir_ips<-"/Users/jonathanmosedale/Data/pest_risks/ips"
dir_ag<-"/Users/jonathanmosedale/Data/pest_risks/ag"
dir_dm<-"/Users/jonathanmosedale/Data/pest_risks/dm"

years<-c(2011:2020)

# Characterise the years - weather gdd10 same for all scenarios
# weather_gdd10, weather_max,weather_min,weather_fdays
r<-rast(file.path(dir_clim,"climate_nsp_s1_2010.tif"))
names(r)
plot(r[['weather_gdd10']])

## Create raster time series for weather variables
wgdd10<-rast()
wtmax<-rast()
wtmin<-rast()
wfdays<-rast()
for(yr in years){
  r<-rast(file.path(dir_clim,paste0("climate_nsp_s1_",yr,".tif")))
  wgdd10<-c(wgdd10,r[['weather_gdd10']])
  wtmax<-c(wtmax,r[['weather_max']])
  wtmin<-c(wtmin,r[['weather_min']])
  wfdays<-c(wfdays,r[['weather_fdays']])
}
names(wgdd10)<-years
names(wtmax)<-years
names(wtmin)<-years
names(wfdays)<-years

# Get mean timeseries
mean_timeseries<-function(r){
  out<-(global(r,mean,na.rm=T))
  out$year<-rownames(out)
  plot(out$mean~out$year)
  return(out)
}
mn_wgdd10<-mean_timeseries(wgdd10)
mn_wtmax<-mean_timeseries(wtmax)
mn_wtmin<-mean_timeseries(wtmin)
mn_wfdays<-mean_timeseries(wfdays)

# Line plots
plot_line<-function(dat,title){
  dat %>%
  tail(10) %>%
  ggplot( aes(x=year, y=mean)) +
  geom_line(group = 1) +
  geom_point() +
  theme_ipsum() +
  ggtitle(title)
}
plot_line(mn_wgdd10,"UK mean GDD10") # Extreme years 2012(low) and 2018(high)
plot_line(mn_wtmax,"UK mean tmax")
plot_line(mn_wtmin,"UK mean tmin")
plot_line(mn_wfdays,"UK mean frostdays")

# Plot rasters for extreme years - common scale
r1<-rast(file.path(dir_clim,paste0("climate_nsp_s1_",2012,".tif")))
r2<-rast(file.path(dir_clim,paste0("climate_nsp_s1_",2018,".tif")))
r<-c(r1[['weather_gdd10']],r2[['weather_gdd10']])
names(r)<-c(2012,2018)
rg<-c(min(global(r,min,na.rm=T)),max(global(r,max,na.rm=T)))
plot( r,range=c(rg[1],rg[2]) )


###### D micans ###################
r12s0<-rast(file.path(dir_dm,"dm_tt_nsp_s0_2012.tif"))
r12s1<-rast(file.path(dir_dm,"dm_tt_nsp_s1_2012.tif"))
r12s4<-rast(file.path(dir_dm,"dm_tt_nsp_s4_2012.tif"))
r18s0<-rast(file.path(dir_dm,"dm_tt_nsp_s0_2018.tif"))
r18s1<-rast(file.path(dir_dm,"dm_tt_nsp_s1_2018.tif"))
r18s4<-rast(file.path(dir_dm,"dm_tt_nsp_s4_2018.tif"))

r18s2<-rast(file.path(dir_dm,"dm_tt_nsp_s2_2018.tif"))
r18s3<-rast(file.path(dir_dm,"dm_tt_nsp_s3_2018.tif"))
r18s5<-rast(file.path(dir_dm,"dm_tt_nsp_s5_2018.tif"))

r12s2<-rast(file.path(dir_dm,"dm_tt_nsp_s2_2012.tif"))
r12s3<-rast(file.path(dir_dm,"dm_tt_nsp_s3_2012.tif"))
r12s5<-rast(file.path(dir_dm,"dm_tt_nsp_s5_2012.tif"))

names(r18s5)

# Egg hatching day of year - coolest
v<-'dnd_hatch_doy'
v<-'dnd_pupate_doy'
v<-"dnd_stages_complete"
v<-"dnd_incomplete_stage"


### Plot extreme scenarios and weather
r<-c(r18s0[[v]],r18s1[[v]],r18s4[[v]],
     r12s0[[v]],r12s1[[v]],r12s4[[v]])

rg<-c(min(global(r,min,na.rm=T)),max(global(r,max,na.rm=T)))

#plot( r, col=rev(map.pal("magma",10)), range=c(rg[1],rg[2]), axes=FALSE, main="", mar=c(0, 0, 1.5, 0), legend=FALSE )

# Plot stage completed
levels(r[[1]]) <- data.frame(id=0:3, Overwintering_stage=c("Egg","Larva", "Pupa", "Adult"))
levels(r[[2]]) <- data.frame(id=0:3, Overwintering_stage=c("Egg","Larva", "Pupa", "Adult"))
levels(r[[3]]) <- data.frame(id=0:3, Overwintering_stage=c("Egg","Larva", "Pupa", "Adult"))
levels(r[[4]]) <- data.frame(id=0:3, Overwintering_stage=c("Egg","Larva", "Pupa", "Adult"))
is.factor(r)
names(r)<-c("S0 2018","S1 2018","S4 2018","S0 2012","S1 2012","S4 2012")

#plot( r, type="classes",axes=FALSE,col=map.pal("viridis",4), main=names(r),mar=c(0, 0, 0, 0),legend=TRUE )

panel(r, nc=2,nr=2,col = c("#440154","#6A00A8", "#B12A90", "#E16462", "#FCA636", "#F0F921"),
      plg=list(title.cex=2,size=2,x="top"),axes=FALSE,main=names(r))


### PLot all scenarios for a year
# hottest year
r<-c(r18s0[[v]],r18s1[[v]],r18s2[[v]],
     r18s3[[v]],r18s4[[v]])
levels(r[[1]]) <- data.frame(id=0:3, Overwintering_stage=c("Egg","Larva", "Pupa", "Adult"))
levels(r[[2]]) <- data.frame(id=0:3, Overwintering_stage=c("Egg","Larva", "Pupa", "Adult"))
levels(r[[3]]) <- data.frame(id=0:3, Overwintering_stage=c("Egg","Larva", "Pupa", "Adult"))
levels(r[[4]]) <- data.frame(id=0:3, Overwintering_stage=c("Egg","Larva", "Pupa", "Adult"))
levels(r[[5]]) <- data.frame(id=0:3, Overwintering_stage=c("Egg","Larva", "Pupa", "Adult"))
#levels(r[[6]]) <- data.frame(id=0:3, Overwintering_stage=c("Egg","Larva", "Pupa", "Adult"))
is.factor(r)
names(r)<-c("S0 2018","S1 2018","S2 2018","S3 2018","S4 2018")

panel(r, nc=5,nr=1,col = c("#440154","#6A00A8", "#B12A90", "#E16462", "#FCA636", "#F0F921"),
      plg=list(title.cex=2,size=2,x="top"),axes=FALSE,main=names(r))

# coolest year
r<-c(r12s0[[v]],r12s1[[v]],r12s2[[v]],
     r12s3[[v]],r12s4[[v]])
levels(r[[1]]) <- data.frame(id=0:3, Overwintering_stage=c("Egg","Larva", "Pupa", "Adult"))
levels(r[[2]]) <- data.frame(id=0:3, Overwintering_stage=c("Egg","Larva", "Pupa", "Adult"))
levels(r[[3]]) <- data.frame(id=0:3, Overwintering_stage=c("Egg","Larva", "Pupa", "Adult"))
levels(r[[4]]) <- data.frame(id=0:3, Overwintering_stage=c("Egg","Larva", "Pupa", "Adult"))
levels(r[[5]]) <- data.frame(id=0:3, Overwintering_stage=c("Egg","Larva", "Pupa", "Adult"))
#levels(r[[6]]) <- data.frame(id=0:3, Overwintering_stage=c("Egg","Larva", "Pupa", "Adult"))
is.factor(r)
names(r)<-c("S0 2012","S1 2012","S2 2012","S3 2012","S4 2012")

panel(r, nc=5,nr=1,col = c("#440154","#6A00A8", "#B12A90", "#E16462", "#FCA636", "#F0F921"),
      plg=list(title.cex=2,size=2,x="top"),axes=FALSE,main=names(r))


#### Plot average development over all years  - extreme years and scenarios
outlist<-list()
for(s in c(0,1,2,3,4)){
  rlist<-list()
  for(y in c(2011:2020)){
    rin<-rast(file.path(dir_dm,paste0("dm_tt_nsp_s",s,"_",y,".tif")))
    rlist<-c(rlist,rin[["dnd_stages_complete"]])
  }
  rmean<-mean(rast(rlist),na.rm=TRUE)
  rmin<-min(rast(rlist),na.rm=TRUE)
  rmax<-max(rast(rlist),na.rm=TRUE)
  outlist<-c(outlist,c(rmin,rmean,rmax))
  panel(c(rmin,rmean,rmax),nc=3,nr=1,axes=FALSE,
        range=c(0,3), breaks=c(0,0.5,1,1.5,2,2.5,3),
        main=paste("Scenario",s,c("min","mean","max")),
        plg=list(labels=c("0 (Egg)","0.5","1.0 (Larva)","1.5","2.0 (Pupa)","2.5","3.0 (Adult)")) )
}

out.r<-rast(outlist)
out.r<-out.r[[c(1,4,7,10,13,2,5,8,11,14,3,6,9,12,15)]]
names(out.r)<-c()
panel(out.r[[c(2,5,8,11,14,1,4,7,10,13,3,6,9,12,15)]],nc=5,nr=3,axes=FALSE,
      range=c(0,3), breaks=c(0,0.5,1,1.5,2,2.5,3),
      plg=list(loc="topright",labels=c("0 (Egg)","0.5","1.0 (Larva)","1.5","2.0 (Pupa)","2.5","3.0 (Adult)")) )





################### ###################I typog ###################  ###################
r12s0<-rast(file.path(dir_ips,"ip_tt_nsp_s0_2012.tif"))
r12s1<-rast(file.path(dir_ips,"ip_tt_nsp_s1_2012.tif"))
r12s4<-rast(file.path(dir_ips,"ip_tt_nsp_s4_2012.tif"))
r18s0<-rast(file.path(dir_ips,"ip_tt_nsp_s0_2018.tif"))
r18s1<-rast(file.path(dir_ips,"ip_tt_nsp_s1_2018.tif"))
r18s4<-rast(file.path(dir_ips,"ip_tt_nsp_s4_2018.tif"))

r18s2<-rast(file.path(dir_ips,"ip_tt_nsp_s2_2018.tif"))
r18s3<-rast(file.path(dir_ips,"ip_tt_nsp_s3_2018.tif"))
r18s5<-rast(file.path(dir_ips,"ip_tt_nsp_s5_2018.tif"))

r12s2<-rast(file.path(dir_ips,"ip_tt_nsp_s2_2012.tif"))
r12s3<-rast(file.path(dir_ips,"ip_tt_nsp_s3_2012.tif"))
r12s5<-rast(file.path(dir_ips,"ip_tt_nsp_s5_2012.tif"))
names(r18s5)

v<-'ips_emerge_doy'
v<-'ips_lay_doy'
v<-"ips_adult_doy" # g1
v<-"ips_g2_lay_doy"
v<-"ips_g2_adult_doy"
v<-"ips_g3_lay_doy"
v<-"ips_g3_adult_doy"
v<-"generations_complete"

### Test plot of extremes
r<-c(r18s0[[v]],r12s0[[v]],r18s5[[v]],r12s5[[v]])
unique(r)
levels(r[[1]]) <- data.frame(id=0:2, Number_generations=c("None","One", "Two"))
levels(r[[2]]) <- data.frame(id=0:2, Number_generations=c("None","One", "Two"))
levels(r[[3]]) <- data.frame(id=0:2, Number_generations=c("None","One", "Two"))
levels(r[[4]]) <- data.frame(id=0:2, Number_generations=c("None","One", "Two"))
is.factor(r)
names(r)<-c("S0 2018","S0 2012","S5 2018","S5 2012")

panel(r, nc=5,nr=1,col = c("#440154","#FCA636",  "#F0F921"),
      plg=list(title.cex=2,size=2,x="top"),axes=FALSE,main=names(r))


# ALL scenarios hottest year
r<-c(r18s0[[v]],r18s1[[v]],r18s2[[v]],
     r18s3[[v]],r18s4[[v]])
levels(r[[1]]) <- data.frame(id=0:2, Number_generations=c("None","One", "Two"))
levels(r[[2]]) <- data.frame(id=0:2, Number_generations=c("None","One", "Two"))
levels(r[[3]]) <- data.frame(id=0:2, Number_generations=c("None","One", "Two"))
levels(r[[4]]) <- data.frame(id=0:2, Number_generations=c("None","One", "Two"))
levels(r[[5]]) <- data.frame(id=0:2, Number_generations=c("None","One", "Two"))
#levels(r[[6]]) <- data.frame(id=0:3, Overwintering_stage=c("Egg","Larva", "Pupa", "Adult"))
is.factor(r)
names(r)<-c("S0 2018","S1 2018","S2 2018","S3 2018","S4 2018")

panel(r, nc=5,nr=1,col = c("#440154","#FCA636",  "#F0F921"),
      plg=list(title.cex=2,size=2,x="top"),axes=FALSE,main=names(r))

# coolest year
r<-c(r12s0[[v]],r12s1[[v]],r12s2[[v]],
     r12s3[[v]],r12s4[[v]]``)
levels(r[[1]]) <- data.frame(id=0:2, Number_generations=c("None","One", "Two"))
levels(r[[2]]) <- data.frame(id=0:2, Number_generations=c("None","One", "Two"))
levels(r[[3]]) <- data.frame(id=0:2, Number_generations=c("None","One", "Two"))
levels(r[[4]]) <- data.frame(id=0:2, Number_generations=c("None","One", "Two"))
levels(r[[5]]) <- data.frame(id=0:2, Number_generations=c("None","One", "Two"))
is.factor(r)
names(r)<-c("S0 2012","S1 2012","S2 2012","S3 2012","S4 2012")

panel(r, nc=5,nr=1,col = c("#440154","#FCA636",  "#F0F921"),
      plg=list(title.cex=2,size=2,x="top"),axes=FALSE,main=names(r))


#### Plot average development over all years  - extreme years and scenarios
outlist<-list()
for(s in c(0,1,2,3,4)){
  rlist<-list()
  for(y in c(2011:2020)){
    rin<-rast(file.path(dir_ips,paste0("ip_tt_nsp_s",s,"_",y,".tif")))
    rlist<-c(rlist,rin[["generations_complete"]])
  }
  rmean<-mean(rast(rlist),na.rm=TRUE)
  rmin<-min(rast(rlist),na.rm=TRUE)
  rmax<-max(rast(rlist),na.rm=TRUE)
  outlist<-c(outlist,c(rmin,rmean,rmax))
  panel(c(rmin,rmean,rmax),nc=3,nr=1,axes=FALSE,col = c("#440154","#FCA636",  "#F0F921"),
        range=c(0,3), breaks=c(0,0.5,1,1.5,2,2.5,3),
        main=paste("Scenario",s,c("min","mean","max")) )
        #plg=list(labels=c("0 (Egg)","0.5","1.0 (Larva)","1.5","2.0 (Pupa)","2.5","3.0 (Adult)")) )
}

out.r<-rast(outlist)
names(out.r)<-c()
panel(out.r[[c(2,5,8,11,14,1,4,7,10,13,3,6,9,12,15)]],nc=5,nr=3,axes=FALSE,
      range=c(0,2), breaks=c(0,0.5,1,1.5,2))
      ,
      plg=list(loc="topright",labels=c("0 (Egg)","0.5","1.0 (Larva)","1.5","2.0 (Pupa)","2.5","3.0 (Adult)")) )

# NUmber of generations
par(mfrow=c(2,2),mar=c(0, 0, 0, 0))
plot(r[[1]], col=c("#440154","#21908C", "#FDE725" ), range=c(rg[1],rg[2]), axes=FALSE,box=FALSE, legend=TRUE)
plot(r[[2]], col=c( "#440154","#21908C"), range=c(rg[1],rg[2]), axes=FALSE,box=FALSE, legend=FALSE)
plot(r[[3]], col=c("#440154","#21908C", "#FDE725" ), range=c(rg[1],rg[2]), axes=FALSE,box=FALSE, legend=FALSE)
plot(r[[4]], col=c( "#440154","#21908C", "#FDE725" ), range=c(rg[1],rg[2]), axes=FALSE,box=FALSE, legend=FALSE)



################### ###################I typog ###################  ###################






##### Climate data
# D micans
r12s1<-rast(file.path(dir_dm,"dm_tt_nsp_s1_2012.tif"))
r18s1<-rast(file.path(dir_dm,"dm_tt_nsp_s1_2018.tif"))

r12w<-rast(file.path(dir_dm,"dm_wt_nsp_s5_2012.tif"))
r18w<-rast(file.path(dir_dm,"dm_wt_nsp_s5_2018.tif"))
v<-'dnd_hatch_doy'
v<-'dnd_pupate_doy'
v<-"dnd_stages_complete"
v<-"dnd_incomplete_stage"

r<-c(r12w[[v]],r12s1[[v]],r18w[[v]],r18s1[[v]])
names(r)<-c("2012 weather","2012 model","2018 weather","2018 model")
rg<-c(min(global(r,min,na.rm=T)),max(global(r,max,na.rm=T)))
par(mar = c(0, 0, 0, 0))

levels(r[[1]]) <- data.frame(id=0:3, Weather_2012=c("Egg","Larva", "Pupa", "Adult"))
levels(r[[2]]) <- data.frame(id=0:3, Model_2012=c("Egg","Larva", "Pupa", "Adult"))
levels(r[[3]]) <- data.frame(id=0:3, Weather_2018=c("Egg","Larva", "Pupa", "Adult"))
levels(r[[4]]) <- data.frame(id=0:3, Model_2018=c("Egg","Larva", "Pupa", "Adult"))
is.factor(r)

plot( r, type="classes",axes=FALSE,col=map.pal("viridis",4),
      main=names(r),mar=c(0, 0, 0, 0),legend=FALSE )
plot( r, type="classes",axes=FALSE,col=map.pal("viridis",4),
      main=names(r),mar=c(0, 0, 2, 0),legend=TRUE )




# I typog
r12s1<-rast(file.path(dir_ips,"ip_tt_nsp_s1_2012.tif"))
r18s1<-rast(file.path(dir_ips,"ip_tt_nsp_s1_2018.tif"))
r12w<-rast(file.path(dir_ips,"ip_tt_nsp_s2_2012.tif"))
r18w<-rast(file.path(dir_ips,"ip_tt_nsp_s2_2018.tif"))
v<-'ips_emerge_doy'
v<-'ips_lay_doy'
v<-"ips_adult_doy" # g1
v<-"ips_g2_lay_doy"
v<-"ips_g2_adult_doy"
v<-"ips_g3_lay_doy"
v<-"ips_g3_adult_doy"
v<-"generations_complete"

r<-c(r12w[[v]],r12s1[[v]],r18w[[v]],r18s1[[v]])
names(r)<-c("2012 weather","2012 model","2018 weather","2018 model")
rg<-c(min(global(r,min,na.rm=T)),max(global(r,max,na.rm=T)))
par(mar = c(0, 0, 0, 0))

# NUmber of generations
names(r)<-c("Weather_2012","Model_2012","Weather_2012","Model_2012")
par(mfrow=c(2,2),mar=c(0, 0, 0, 0))
plot(r[[1]], col=c("#440154","#21908C", "#FDE725" ), range=c(rg[1],rg[2]), axes=FALSE,box=FALSE, legend=FALSE)
plot(r[[2]], col=c("#440154","#21908C", "#FDE725" ), range=c(rg[1],rg[2]), axes=FALSE,box=FALSE, legend=FALSE)
plot(r[[3]], col=c("#440154","#21908C", "#FDE725" ), range=c(rg[1],rg[2]), axes=FALSE,box=FALSE, legend=FALSE)
plot(r[[4]], col=c( "#440154","#21908C", "#FDE725" ), range=c(rg[1],rg[2]), axes=FALSE,box=FALSE, legend=FALSE)


# All scenarios
r18s1<-rast(file.path(dir_ips,"ip_tt_nsp_s1_2018.tif"))
#r18s2-rast(file.path(dir_ips,"ip_tt_nsp_s2_2012.tif"))
r18s3<-rast(file.path(dir_ips,"ip_tt_nsp_s3_2018.tif"))
#r18s4<-rast(file.path(dir_ips,"ip_tt_nsp_s4_2018.tif"))
r18s5<-rast(file.path(dir_ips,"ip_tt_nsp_s5_2018.tif"))
r18w<-rast(file.path(dir_ips,"ip_tt_nsp_s2_2018.tif"))

r<-c(r18s1[[v]],r18s3[[v]],r18s5[[v]],r18w[[v]])
names(r)<-c("2018 S1","2018 S3","2018 S5","2018 weather")
rg<-c(min(global(r,min,na.rm=T)),max(global(r,max,na.rm=T)))
par(mar = c(0, 0, 0, 0))
par(mfrow=c(2,2),mar=c(0, 0, 0, 0))
plot(r[[1]], col=c("#440154","#21908C", "#FDE725" ), range=c(rg[1],rg[2]), axes=FALSE,box=FALSE, legend=FALSE)
plot(r[[2]], col=c("#440154","#21908C", "#FDE725" ), range=c(rg[1],rg[2]), axes=FALSE,box=FALSE, legend=FALSE)
plot(r[[3]], col=c("#440154","#21908C", "#FDE725" ), range=c(rg[1],rg[2]), axes=FALSE,box=FALSE, legend=FALSE)
plot(r[[4]], col=c( "#440154","#21908C", "#FDE725" ), range=c(rg[1],rg[2]), axes=FALSE,box=FALSE, legend=FALSE)


############# all years to get average
# Ips
r_list<-list()
for(s in 1:5){
  ryrs<-rast()
  for(yr in c(2011:2020)){
    r<-rast(file.path(dir_ips,paste0("ip_tt_nsp_s",s,"_",yr,".tif")))[["generations_complete"]]
    ryrs<-c(ryrs,r)
  }
  r_list<-c(r_list,ryrs)
}

names(r_list)<-c('s1','s2','s3','s4','s5')
rmeans<-rast(lapply(r_list,mean,na.rm=TRUE))
rg<-c(0,2)
par(mfrow=c(1,1),mar=c(1, 0, 0, 0))
plot(rmeans[[1]], col=map.pal("viridis",10), range=c(rg[1],rg[2]), axes=FALSE,box=FALSE, legend=FALSE)
plot(rmeans[[2]], col=map.pal("viridis",10), range=c(rg[1],rg[2]), axes=FALSE,box=FALSE, legend=FALSE)
plot(rmeans[[3]], col=map.pal("viridis",10), range=c(rg[1],rg[2]), axes=FALSE,box=FALSE, legend=FALSE)
plot(rmeans[[4]], col=map.pal("viridis",10), range=c(rg[1],rg[2]), axes=FALSE,box=FALSE, legend=FALSE)
plot(rmeans[[5]], col=map.pal("viridis",10), range=c(rg[1],rg[2]), axes=FALSE,box=FALSE, legend=TRUE)


# Dm
r_list<-list()
for(s in 1:5){
  ryrs<-rast()
  for(yr in c(2011:2020)){
    r<-rast(file.path(dir_dm,paste0("dm_tt_nsp_s",s,"_",yr,".tif")))[["dnd_stages_complete"]]
    ryrs<-c(ryrs,r)
  }
  r_list<-c(r_list,ryrs)
}

names(r_list)<-c('s1','s2','s3','s4','s5')
rmeans<-rast(lapply(r_list,mean,na.rm=TRUE))
rg<-c(0,4)
par(mfrow=c(1,1),mar=c(1, 0, 0, 0))
plot(rmeans[[1]], col=map.pal("viridis",10), range=c(rg[1],rg[2]), axes=FALSE,box=FALSE, legend=TRUE)
plot(rmeans[[2]], col=map.pal("viridis",10), range=c(rg[1],rg[2]), axes=FALSE,box=FALSE, legend=FALSE)
plot(rmeans[[3]], col=map.pal("viridis",10), range=c(rg[1],rg[2]), axes=FALSE,box=FALSE, legend=FALSE)
plot(rmeans[[4]], col=map.pal("viridis",10), range=c(rg[1],rg[2]), axes=FALSE,box=FALSE, legend=FALSE)
plot(rmeans[[5]], col=map.pal("viridis",10), range=c(rg[1],rg[2]), axes=FALSE,box=FALSE, legend=FALSE)




# Convert the suitability change raster to a data frame for plotting
r_df <- terra::as.data.frame(c(r[['weather_gdd10']],r2[['weather_gdd10']]), xy = TRUE, na.rm = TRUE)
names(r_df) <- c("x", "y", "vals")
rg<-range(r_df$vals)
# Create the map
p <- ggplot() +
  #geom_sf(data = world, fill = "white", color = "black") +  # Add country borders
  geom_raster(data = r_df, aes(x = x, y = y, fill = vals)) +
  scale_fill_viridis_c()
  #scale_fill_gradient2(low = "blue", mid = "green", high = "yellow",limit = c(rg[1], rg[2]), name = "GDD 10") + # midpoint = 0,
  #coord_sf(xlim = c(-12, 45), ylim = c(33, 72), expand = FALSE) +  # Limit the map to European extent
  labs(title = "something",
       x = "x", y = "y") +
  theme_minimal(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 20),  # Ensure title is centered and bold
        panel.background = element_rect(fill = "white", colour = "white"),  # Set background to white
        plot.background = element_rect(fill = "white", colour = "white"),  # Ensure surrounding background is also white
        legend.background = element_rect(fill = "white"))  # Ensure legend background is white

