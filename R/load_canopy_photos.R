#### Functions concerning LAI nad photo analysis

get_image_metadata<-function(filepath){
  filename<-basename(filepath)
  site<-substr(filename,1,2)
  spp<-substr(filename,4,6)
  tree<-substr(filename,8,8)
  asp<-substr(filename,10,10)
  height<-substr(filename,12,14)
  date<-substr(filename,16,21)
  date<-as.character(as.Date(date,"%d%m%y"))
  treeid<-tolower(paste(site,spp,tree,sep='_'))
  return(list(filepath=filepath,filename=filename,date=date,treeid=treeid,site=site,spp=spp,tree=tree,asp=asp,height=height))
}


read_lai_results<-function(filepath){
  #print(filepath)
  if(file.exists(filepath)==FALSE) stop('File NOT found!!!')

  # Open connection and read lines
  con<-file(filepath,open="r")
  txtlines<-readLines(con,52)
  close(con)

  # Extract relevant info
  ringnum<-as.numeric(str_split_i(txtlines[11],'\t',i=2))
  ringw<-as.numeric(str_sub(str_split_i(txtlines[11],'\t',i=3),end=-2))
  sectornum<-as.numeric(str_split_i(txtlines[12],'\t',i=2))
  s1azi<-as.numeric(str_split_i(txtlines[20],'\t',i=3))
  s2azi<-as.numeric(str_split_i(txtlines[21],'\t',i=3))

  # Types of LAI
  s1lai_licor<-as.numeric(str_split_i(txtlines[38],'\t',i=3))
  s2lai_licor<-as.numeric(str_split_i(txtlines[39],'\t',i=3))
  s1lai_g<-as.numeric(str_split_i(txtlines[41],'\t',i=3))
  s2lai_g<-as.numeric(str_split_i(txtlines[42],'\t',i=3))
  s1lai_t<-as.numeric(str_split_i(txtlines[44],'\t',i=3))
  s2lai_t<-as.numeric(str_split_i(txtlines[45],'\t',i=3))

  # gaps - remove %
  transm<-as.numeric(str_sub(str_split_i(txtlines[29],'\t',i=2),end=-2))
  skyview<-as.numeric(str_sub(str_split_i(txtlines[31],'\t',i=2),end=-2))
  lggaps<-as.numeric(str_sub(str_split_i(txtlines[32],'\t',i=2),end=-2))

  Fmv<-as.numeric(str_sub(str_split_i(txtlines[50],'\t',i=3),end=-2))
  #Frv<-as.numeric(str_split_i(txtlines[50,'\t',i=9))

  # Clumping omega by ring
  clump1<-as.numeric(str_split_i(txtlines[15],'\t',i=9))
  clump2<-as.numeric(str_split_i(txtlines[16],'\t',i=9))
  clump3<-as.numeric(str_split_i(txtlines[17],'\t',i=9))
  clump<-mean(c(clump1,clump2,clump3),na.rm=TRUE)

  # Leaf angles - use T. & al. (2010)
  s1angle<-as.numeric(str_split_i(txtlines[44],'\t',i=4))
  s2angle<-as.numeric(str_split_i(txtlines[45],'\t',i=4))

  # Metadata from filename
  meta<-get_image_metadata(filepath)

  # Which sector and lai and angle to use
  if(meta['asp']=='n') azi<-0
  if(meta['asp']=='s') azi<-180
  lai_lc<-c(s1lai_licor,s2lai_licor)[which(c(s1azi,s2azi)==azi)]
  lai_g<-c(s1lai_g,s2lai_g)[which(c(s1azi,s2azi)==azi)]
  lai_t<-c(s1lai_t,s2lai_t)[which(c(s1azi,s2azi)==azi)]

  angle<-c(s1angle,s2angle)[which(c(s1azi,s2azi)==azi)]

  return(c(meta,lai_lc=lai_lc,lai_g=lai_g,lai_t=lai_t,angle=angle,transm=transm,skyview=skyview,lggaps=lggaps,
           Fmv=Fmv,clump1=clump1,clump2=clump2,clump3=clump3,clump=clump))
}

plot_param<-function(df,site,height,param=c('angle','clump')){
  v<-paste0(param,'_',height)
  minv<-min(df[,v])
  maxv<-max(df[,v])

  v<-sym(v)
  plotdata<-df[which(df$site==site),]
  plotdata<-plotdata[order(plotdata$treeid),]
  plotdata$label<-paste0(plotdata$spp,'_',plotdata$tree)

  fig<-ggplot(data = plotdata, aes(x = factor(date), y = !!v, color = spp)) +
    ggtitle(paste(param,"estimate from site",site)) + ylim(minv,maxv) +
    xlab("Date") + ylab(paste(param,'at',height,'metres')) +
    geom_line(aes(group = treeid)) + geom_point() +
    geom_dl(aes(label=label),method=list(dl.combine("first.points","last.points"),cex=0.7))
  return(fig)
}
