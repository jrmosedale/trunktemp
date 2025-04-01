library(devtools)
load_all()

library(reticulate)
library(rgee)
library(terra)
library(sf)
library(microclimdata)
# Set up python and GEE
pathtopython<-"/usr/local/bin/python3.12"
pathtopython<-"/Library/Frameworks/Python.framework/Versions/3.12/Python"

reticulate::use_python(pathtopython)
ee_Authenticate()

GoogleDrivefolder<-"ee-treeheight"
projectname<-"ee-treeheight"

# Get site info
dir_od<-"/Users/jonathanmosedale/Library/CloudStorage/OneDrive-UniversityofExeter"
dir_field<-file.path(dir_od,"MetOff_data/field_data")
dir_dtm<-file.path(dir_od,"Data/Terrain50")
dtm_file<-file.path(dir_dtm,"uk_dtm.tif")


# Get veg heights for each site
# local v
vegheight_download<-function(r, GoogleDrivefolder, pathtopython, projectname = NA, silent = FALSE) {
  #use_python(paste0(pathtopython,"python.exe"), required = TRUE)
  if (is.na(projectname) == FALSE) ee$Initialize(project=projectname)
  # Get bounding box
  e<-ext(r)
  r2<-rast(e)
  crs(r2)<-crs(r)
  r2<-project(r2,"EPSG:4326")
  e<-ext(r2)
  # work out crs
  proj_string <- crs(r, describe=T)
  epsg_code <- paste0("EPSG:",proj_string$code)
  # get aoi_coordinates
  aoi <- ee$Geometry$Rectangle(c(e$xmin, e$ymin, e$xmax, e$ymax))
  aoi_bounds <- aoi$bounds()$getInfo()
  aoi_coordinates <- aoi_bounds$coordinates[[1]]
  canopy_height <- ee$Image('users/nlang/ETH_GlobalCanopyHeight_2020_10m_v1')
  task <- ee$batch$Export$image$toDrive(
    image = canopy_height,  # The image to export
    description = 'canopy_height_export',  # Description for the export task
    folder = GoogleDrivefolder,  # Folder in Google Drive where the file will be saved
    fileNamePrefix = 'canopy_height_2020',  # Prefix for the filename
    region =   aoi_coordinates,  # Define the export region
    scale = 10,  # Scale in meters (resolution)
    crs = epsg_code  # Coordinate reference system
  )
  task$start()
  if (silent == FALSE) .monitor_task(task$id)
}
dtm<-rast(dtm_file)


# Get veg heights
site_file<-file.path(dir_field,"site_properties.csv")
tree_file<-file.path(dir_field,"Field_data_locations.csv")

tree_locs<-read.csv(tree_file)
trees_sf$id<-paste(tolower(trees_sf$Site),trees_sf$Tree_species,trees_sf$Tree,sep='_')
trees_sf<-st_as_sf(tree_locs,coords = c("lon", "lat"), crs = 4326)
trees_27700<-st_transform(trees_sf,27700)

site_props<-read.csv(site_file,strip.white = TRUE )
sites_sf<-st_as_sf(site_props,coords = c("lon", "lat"), crs = 4326)
sites_27700<-st_transform(sites_sf,27700)

plot(crop(dtm,e))


for(s in sites_27700$site){
  trees<-trees_27700[which(trees_27700$Site==s),]
  e<-ext(st_bbox(st_buffer(trees,100)))
  plot(crop(dtm,e),main=s); plot(vect(trees),add=T,col='red')
  r <- rast(e)
  crs(r) <- "EPSG:27700"
  vegheight_download(r, GoogleDrivefolder, pathtopython, projectname)
}


vhgt <- rast(file.path(dir_field,'treeheight','cw_canopy_height_2020.tif'))
plot(vhgt)

# Extract tree heights
tree_locs<-read.csv(tree_file)
trees_sf<-st_as_sf(tree_locs,coords = c("lon", "lat"), crs = 4326)
trees_sf$id<-paste(tolower(trees_sf$Site),trees_sf$Tree_species,trees_sf$Tree,sep='_')
trees_27700<-st_transform(trees_sf,27700)

height_list<-list()
for (site in c('AH','BE','CW','KG','WH')){
  vhgt <- rast(file.path(dir_field,'treeheight',paste0(tolower(site),'_canopy_height_2020.tif')))
  names(vhgt)<-'Tree_Height'
  sitetrees<-trees_27700[which(trees_27700$Site==site),]
  hght<-extract(vhgt,sitetrees)
  hght$ID<-sitetrees$id
  height_list[[site]]<-hght
}
df<-do.call(rbind,height_list)
