## code to prepare `DATASET` dataset goes here
library(terra)
dir_od<-"/Users/jonathanmosedale/Library/CloudStorage/OneDrive-UniversityofExeter"
dir_loggers<-file.path(dir_od,"MetOff_data/field_data/Logger_results")
dir_weather<-file.path(dir_od,"MetOff_data/field_data/weather")
dir_lai<-file.path(dir_od,"MetOff_data/field_data/LAI_results")


# Air logger data
stagdata_df<-readRDS(file.path(dir_loggers,'surveytag_air_logger_results.Rds'))

# Tree logger data
ttagdata_df<-readRDS(file.path(dir_loggers,'tinytag_bark_logger_results.Rds'))

# Weather data
for(site in c('AH','BE','CW','KG','WH'))
weather<-(file.path(dir_weather,paste0("weather_",site,"_2024.Rds")))

# Also weather byu sub-site - CW2 etc .csv


# LAI results
lai_models<-readRDS(file.path(dir_lai,'lc_lai_models.Rds'))
lairesults_df<-read.csv(file.path(dir_lai,'lai_allsites_results.csv'),sep=',')

usethis::use_data(lai_models, overwrite = TRUE)

# Merge loggers and weather temp to same hourly timesteps
# colnames<-c("id", "site","spp","tree","datetime","air_temp","n_outer_temp","n_inner_temp","s_outer_temp","s_inner_temp")
#allresults_df<-readRDS(file.path(dir_loggers,'logger_results.Rds'))
loggers_df<-readRDS(file.path(dir_loggers,'logger_results_stedit.Rds'))
usethis::use_data(loggers_df, overwrite = TRUE)


# Field tree information combined
usethis::use_data(trees_sf, overwrite = TRUE)
usethis::use_data(sites_sf, overwrite = TRUE)
