#### Function to create overlapping tileset
# Return tile extens as list and also a y/n character vector of whether tile includes non NA
# template.r = land mask to inform whether land present
#tileset<-create_overlapping_tiles(gb)
#elist<-tileset$tile_extents
#etype<-tileset$tile_land
create_overlapping_tiles<-function(template.r,overlap=5000,sz=75000){
  xmax<-ext(template.r)[2]
  xmin<-ext(template.r)[1]
  ymax<-ext(template.r)[4]
  ymin<-ext(template.r)[3]
  xstart<-seq(xmin,xmax,75000)
  xend<-c(xstart[2:length(xstart)]+overlap,as.numeric(xmax))
  ystart<-seq(ymin,ymax,75000);
  yend<-c(ystart[2:length(ystart)]+overlap,as.numeric(ymax))

  elist<-list()
  etype<-c()
  for(x in 1:length(xstart)){
    for(y in 1:length(ystart)){
      e<-ext(xstart[x],xend[x],ystart[y],yend[y])
      r<-crop(template.r,e)
      elist<-c(elist,e)
      if(!all(is.na(values(r)))){
        etype<-c(etype,'y')
      } else etype<-c(etype,'n')
    }
  }
  #length(which(etype=='y'))
  #missing_list<-list()
  #for(e in elist[which(etype=='n')]){
  #  missing_list<-c(missing_list,crop(template.r,e))
  #}
  tileset<-list("tile_extents"=elist,"tile_land"=etype)
  return(tileset)
}


### FUNCTIONS from: https://github.com/ilyamaclean/microclimf/tree/main
#' @title Mosaics a list of overlapping SpatRasters blending overlap areas
#' @description Mosaics a list of overlapping SpatRasters blending
#' the areas of overlap using a distance weighting to eliminate tiling effects
#' @param rlist a list of SpatRasters
#' @details
#' If rlist contains SpatRasters that are not
#' overlapping the conventional terra::moasic function is used.
#' If rlist contains SpatRasters that do overlap, they should comprise
#' a list of adjacent rasters in a single row or column.
#' @import terra
#' @export
mosaicblend <- function(rlist) {
  # order by row and then by column
  xmn<-0
  ymn<-0
  for (i in 1:length(rlist)) {
    e<-ext(rlist[[i]])
    xmn[i]<-e$xmin
    ymn[i]<-e$ymin
  }
  xmn1<-unique(xmn)
  ymn1<-unique(ymn)
  le<-min(length(xmn1),length(ymn1))
  if (le > 1) warning("rlist not a row or column. Blended mosaicing may not work")
  if (length(xmn1) > length(ymn1)) {
    o<-order(xmn)
  } else o<-order(ymn)
  rlist2<-list()
  for (i in 1:length(o)) rlist2[[i]]<-rlist[[o[i]]]
  rlist<-NULL
  rma<-rlist2[[1]]
  for (i in 2:length(rlist2)) {
    r<-rlist2[[i]]
    it<-terra::intersect(ext(rma),ext(r))
    a<-as.numeric((it$xmax-it$xmin)*(it$ymax-it$ymin))
    if (a>0) {
      rma<-.blendmosaic(rma, r)
    } else rma<-mosaic(rma,r)
  }
  return(rma)
}
#' blend two adjacent rasters that have overlap
.blendmosaic<-function(r1, r2) {
  # run checks
  reso1<-res(r1)
  reso2<-res(r2)
  if (reso1[1] != reso2[1]) stop("resolutions must match")
  if (reso1[2] != reso2[2]) stop("resolutions must match")
  # Find whether r2 is TT, TR, RR, BR, BB, BL, LL or TL
  e1<-ext(r1)
  e2<-ext(r2)
  corner<-"ID"
  if (e2$ymax > e1$ymax & e2$xmax == e1$xmax) corner<-"TT"
  if (e2$ymax > e1$ymax & e2$xmax > e1$xmax) corner<-"TR"
  if (e2$ymax == e1$ymax & e2$xmax > e1$xmax) corner<-"RR"
  if (e2$ymax < e1$ymax & e2$xmax > e1$xmax) corner<-"BR"
  if (e2$ymax < e1$ymax & e2$xmax == e1$xmax) corner<-"BB"
  if (e2$ymax < e1$ymax & e2$xmax < e1$xmax) corner<-"BL"
  if (e2$ymax == e1$ymax & e2$xmax < e1$xmax) corner<-"LL"
  if (e2$ymax > e1$ymax & e2$xmax < e1$xmax) corner<-"TL"
  if (corner == "ID") {
    ro<-mosaic(r1,r2,fun="mean")
  } else {
    # Calculate overlap area
    if (corner == "TR" || corner == "TT") eo<-ext(e2$xmin,e1$xmax,e2$ymin,e1$ymax)
    if (corner == "BR" || corner == "RR") eo<-ext(e2$xmin,e1$xmax,e1$ymin,e2$ymax)
    if (corner == "BL" || corner == "BB") eo<-ext(e1$xmin,e2$xmax,e1$ymin,e2$ymax)
    if (corner == "TL" || corner == "LL") eo<-ext(e2$xmin,e1$xmax,e2$ymin,e1$ymax)
    # Calculate weights
    nx<-as.numeric((eo$xmax-eo$xmin)/res(r1)[1])
    ny<-as.numeric((eo$ymax-eo$ymin)/res(r1)[2])
    wx<-matrix(rep(seq(0,1,length.out=nx),each=ny),ncol=nx,nrow=ny)
    wy<-matrix(rep(seq(0,1,length.out=ny),nx),ncol=nx,nrow=ny)
    # Create a SpatRast of the blended area
    if (corner == "TT") w1<-wy
    if (corner == "TR") w1<-sqrt((1-wx)^2+wy^2)/sqrt(2)
    if (corner == "RR") w1<-1-wx
    if (corner == "BR") w1<-sqrt((1-wx)^2+(1-wy)^2)/sqrt(2)
    if (corner == "BB") w1<-1-wy
    if (corner == "BL") w1<-sqrt(wx^2+(1-wy)^2)/sqrt(2)
    if (corner == "LL") w1<-wx
    if (corner == "TL") w1<-sqrt(wx^2+wy^2)/sqrt(2)
    # Apply weights to raster
    nn<-dim(r1)[3]
    r1c<-crop(r1,eo)
    r2c<-crop(r2,eo)
    w1<-.rast(.rta(w1,nn),r1c[[1]])
    rb<-r1c*w1+r2c*(1-w1)
    # Clip out the overlap area
    ro<-mosaic(r1,r2)
    re<-w1*0-9999
    ro<-mosaic(ro,re,fun="min")
    ro[ro == -9999]<-NA
    # Mosaic with belnded data
    ro<-mosaic(ro,rb,fun="mean")
  }
  return(ro)
}
#' Convert matrix to array
.rta <- function(r,n) {
  m<-.is(r)
  a<-array(rep(m,n),dim=c(dim(r)[1:2],n))
  a
}

### Wrapper to blend all lad tile files in dir and add allsea tiles
blend_all_tiles<-function(dir_in,dir_out,template_file,scenario_name="syc_s0",year_range="2011_2015"){
  gb1km<-rast(template_file)
  gb<-trim(aggregate(gb1km,25,fun="mean",na.rm=TRUE))
  gb1km<-extend(gb1km,gb)

  scenario_params<-read.csv(scenario_file)
  spp_params<-read.csv(spp_file)
  sp<-substr(scenario_name,1,3)

  # Create overlapping tileset and record which are sea only which with land
  tileset<-create_overlapping_tiles(gb1km,overlap=5000,sz=75000)
  elist<-tileset$tile_extents
  etype<-tileset$tile_land
  landtiles<-which(elist[[which(etype=="y")]])

  tile_filelist<-file.path(dir_in,paste0("alh_",scenario_name,"_1kmrisks_t",landtiles,"_5yr_",year_range,".tif"))
  if(!all(file.exists(tile_filelist))) stop(paste("Missing input tile files:",tile_filelist[which(!file.exists(tile_filelist))]))

  # Load tiles into list of raster tiles, adding sea tiles
  num_layers<-nlyr(rast(tile_filelist[[1]]))
  lyrnames<-names(rast(tile_filelist[[1]]))

  tile_rlist<-list()
  for(n in 1:length(elist)){
    print(n)
    if(etype[n]=="n"){
      empty.r<-crop(gb1km,crop(gb,elist[[n]],snap="out"))
      r<-rep(empty.r,num_layers)
    } else{
      r<-rast(file.path(dir_in,paste0("alh_",scenario_name,"_1kmrisks_t",n,"_5yr_",year_range,".tif")))
    }
    tile_rlist[[n]]<-r
  } # end for

  ### Blend tile list
  numrows<-17 # number of tiles in a column
  startseq<-seq(1,length(elist),numrows)
  endseq<-startseq+(numrows-1)
  # blend rows to columns
  col_list<-list()
  for(c in 1:length(startseq)){
    s<-startseq[c]; e<-endseq[c]
    blend.r<-mosaicblend(rlist=inrtemps_list[s:e])
    col_list<-c(col_list,blend.r)
  }
  # blend columns
  ukresults.r<-mosaicblend(col_list)
  names(ukresults.r)<-lyrnames
  writeRaster(ukresults.r,file.path(dir_out,paste0("alh_",scenario_name,"1kmrisks_5yr_",year_range,".tif")),overwrite=TRUE)
  return(ukresults.r)
}
