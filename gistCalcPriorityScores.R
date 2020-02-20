
#########################################################################################################
#
# Code to accompany "Green Infrastructure Space and Traits (GIST) Model: Integrating green infrastructure 
# spatial placement and plant traits to maximize multifunctionality"
# Urban Forestry & Urban Greening
#
# Tyler J. Tran, Matthew R. Helmus, and Jocelyn E. Behm
# 2020
#
# Contact: TJT (tylerjtran [at] gmail [dot] com); MRH (mrhelmus [at] temple [dot] edu);
# JEB, corresponding author (jebehm [at] temple [dot] edu)
#
# Note: This code does not apply to all of the analyses in the Tran et al. 2020 paper. Because of our 
# data sharing agreement with partners, we are not allowed to share some code. The code below concerns
# calculating condition scores and priority scores for the spatial component of the GIST model. The 
# outputs of each code chunk are aggregates for each benefit at the census tract level. To calculate
# condition scores and priority scores, these aggregates were scaled from 0 to 1 to be relative to each
# other.
#
#########################################################################################################


# To run the code below, download the data folder from URL TO REPOSITORY.
setwd('') # Set path pointing to folder where you downloaded data


library(tidyverse); library(sf)
# library(raster); library(rgdal) # These packages will only be necessary if you choose to do image processing on Landsat imagery for heat island mediation condition scores and priority scores


# Load census tract shapefile
tracts <- st_read('./censusTracts.shp') %>%
  mutate(NAME10 = as.character(NAME10),
         GEOID10 = as.character(GEOID10))


###########################################################
# Stormwater diversion
# The proportional area of impervious surface in census tracts


# Below is land cover data for each census tract, with sums of # pixels of each class for each tract
# Column percentImp is (buildings + roads + otherPaved)/totalPixels
impervious <- read_csv('./landCoverTracts.csv') %>%
  mutate(tract = as.character(tract),
         percentGreen = (tree + grassShrub)/totalPixels)


# Join census tract shapefile to impervious surface data
tracts <- left_join(tracts, impervious, by = c('NAME10' = 'tract'))

tracts <- tracts %>%
  mutate(tractArea = as.numeric(st_area(tracts)))

# Now do max-min scaling to scale values between 0 and 1
tracts <- tracts %>%
  mutate(stormwater01 = (percentImp - min(percentImp, na.rm = T))/
           (max(percentImp, na.rm = T) - min(percentImp, na.rm = T)))


###########################################################
# Biodiversity, 2 scenarios: 
  # Biodiversity-complementation:
  # Biodiversity-equity:


greenspace <- st_read('./greenspace.shp')
# greenspace <- greenspace %>%
#   mutate(area2 = st_area(greenspace))

greenspaceCentroids <- st_centroid(greenspace)

tractGreen <- st_intersection(greenspaceCentroids, tracts) %>%
  group_by(GEOID10) %>%
  summarise(areaGreenspace = sum(area, na.rm = T),
            tractArea = unique(tractArea),
            nPatches = n()) %>%
  mutate(meanPatchSize = areaGreenspace/nPatches,
         logMPS = log(meanPatchSize)) %>%
  as.tibble() %>%
  dplyr::select(GEOID10, meanPatchSize, logMPS)


tracts <- left_join(tracts, tractGreen, by = 'GEOID10') 

# Now do max-min scaling to scale values between 0 and 1
tracts <- tracts %>%
  mutate(biodiversity1_01 = (logMPS - min(logMPS, na.rm = T))/
           (max(logMPS, na.rm = T) - min(logMPS, na.rm = T)),
         biodiversity2_01 = (percentGreen - min(percentGreen, na.rm = T))/
           (max(percentGreen, na.rm = T) - min(percentGreen, na.rm = T)),
         biodiversity2_01 = 1 - biodiversity2_01)


################################################################
# Air quality improvements
# See https://www.epa.gov/air-research/downscaler-model-predicting-daily-air-pollution


pm25 <- read_csv('./airQualityTracts.csv') %>% # daily average PM 2.5 from summer 2013 by census tract
  mutate(FIPS = as.character(FIPS),
         pm25dailyAvg = as.numeric(pm25dailyAvg)) %>%
  group_by(Date, FIPS) %>%
  summarise(meanPM25 = mean(pm25dailyAvg, na.rm = T)) %>%
  # ungroup() %>%
  group_by(FIPS) %>%
  summarise(meanMeanPM25 = mean(meanPM25, na.rm = T))

tracts <- left_join(tracts, pm25, by = c('GEOID10' = 'FIPS'))

# Now do max-min scaling to scale values between 0 and 1
tracts <- tracts %>%
  mutate(airQuality01 = (meanMeanPM25 - min(meanMeanPM25, na.rm = T))/
           (max(meanMeanPM25, na.rm = T) - min(meanMeanPM25, na.rm = T)))



################################################################
# Improved human well-being
# See https://www.cdc.gov/500cities/index.htm

health <- read_csv('./wellBeingTracts.csv') %>%
  filter(StateAbbr == 'PA' & PlaceName == 'Philadelphia') %>%
  mutate(mentalHealth = MHLTH_CrudePrev, # the data are percentages
         physHealth = PHLTH_CrudePrev,
         wellBeing = MHLTH_CrudePrev + PHLTH_CrudePrev) %>%
  dplyr::select(TractFIPS, wellBeing, Population2010)

tracts <- left_join(tracts, health, by = c('GEOID10' = 'TractFIPS'))

# Now do max-min scaling to scale values between 0 and 1
tracts <- tracts %>%
  mutate(wellBeing01 = (wellBeing - min(wellBeing, na.rm = T))/
           (max(wellBeing, na.rm = T) - min(wellBeing, na.rm = T)))


################################################################
# Crime
# Public crime incident data, see https://www.opendataphilly.org/dataset/crime-incidents 

crime <- st_read('./crimeTracts.shp') %>%
  group_by(GEOID10) %>%
  summarise(crime = n()) %>%
  as.tibble() %>%
  dplyr::select(GEOID10, crime)

tracts <- left_join(tracts, crime, by = 'GEOID10')

# Now do max-min scaling to scale values between 0 and 1
tracts <- tracts %>%
  mutate(crime01 = (crime - min(crime, na.rm = T))/
           (max(crime, na.rm = T) - min(crime, na.rm = T)))


################################################################
# Property values: Data from Philadelphia public OPA data
# https://www.opendataphilly.org/dataset/opa-property-assessments

propVal <- read_csv('./propValuesTracts.csv') %>%
  mutate(GEOID10 = as.character(GEOID10)) %>%
  filter(category_1 %in% c('Single Family', 'Multi Family'),
         total_liva > 0) %>%
  mutate(pricePerSqft = market_val/total_liva) %>%
  group_by(GEOID10) %>%
  summarise(meanPropVal = mean(market_val, na.rm = T),
            medianPropVal = median(market_val, na.rm = T),
            meanPricePerSqft = mean(pricePerSqft, na.rm = T))

tracts <- left_join(tracts, propVal, by = 'GEOID10')

# Now do max-min scaling to scale values between 0 and 1
tracts <- tracts %>%
  mutate(propValues01 = (meanPricePerSqft - min(meanPricePerSqft, na.rm = T))/
           (max(meanPricePerSqft, na.rm = T) - min(meanPricePerSqft, na.rm = T)))

################################################################
# Heat island mediation
#
# Much of the code below is adapted from this blog post: https://www.gis-blog.com/calculation-of-land-surface-temperature-lst-from-landsat-8-using-r/
# We have not included large file-size Landsat imagery in our included data. However, the data can be downloaded for free at https://www.usgs.gov/land-resources/nli/landsat/landsat-data-access?qt-science_support_page_related_con=0#qt-science_support_page_related_con


airTemp <- read_csv('./airTempTracts.csv') %>%
  dplyr::select(GEOID10, jjaTemp2016) %>%
  mutate(GEOID10 = as.character(GEOID10))

tracts <- left_join(tracts, airTemp, by = 'GEOID10')

# Now do max-min scaling to scale values between 0 and 1
tracts <- tracts %>%
  mutate(heatIsland01 = (jjaTemp2016 - min(jjaTemp2016, na.rm = T))/
           (max(jjaTemp2016, na.rm = T) - min(jjaTemp2016, na.rm = T)))


# # This function will be used below to pre-process the imagery
# calcLST <- function(RADIANCE_MULT_BAND_10, RADIANCE_MULT_BAND_11, RADIANCE_ADD_BAND_10, RADIANCE_ADD_BAND_11, K1_B10, K1_B11, K2_B10, K2_B11, sinSunAngle,
#                     band4, band5, band10, band11, censusTracts){
#   
#   tract <- spTransform(censusTracts, CRS(proj4string(band_10)))
#   band_10 <- crop(band10, tract)
#   band_11 <- crop(band11, tract)
#   band_4 <- crop(band4, tract)
#   band_5 <- crop(band5, tract)
#   
#   #Calculate TOA from DN:
#   toa_band10 <- calc(band_10, fun=function(x){RADIANCE_MULT_BAND_10 * x + RADIANCE_ADD_BAND_10})
#   toa_band11 <- calc(band_11, fun=function(x){RADIANCE_MULT_BAND_11 * x + RADIANCE_ADD_BAND_11})
#   
#   
#   #Calculate LST in Kelvin for Band 10 and Band 11
#   temp10_kelvin <- calc(toa_band10, fun=function(x){K2_B10/log(K1_B10/x + 1)})
#   temp11_kelvin <- calc(toa_band11, fun=function(x){K2_B11/log(K1_B11/x + 1)})
#   
#   #Convert Kelvin to Celsius for Band 10 and 11
#   temp10_celsius <- calc(temp10_kelvin, fun=function(x){x - 273.15})
#   temp11_celsius <- calc(temp11_kelvin, fun=function(x){x - 273.15})
#   
#   atSat_kelvin <- mean(temp10_kelvin, temp11_kelvin)
#   atSat_celsius <- mean(temp10_celsius, temp11_celsius)
#   
#   # Calculate NDVI
#   red_ref <- ((0.00002*band_4) - 0.1)/sinSunAngle
#   nir_ref <- ((0.00002*band_5) - 0.1)/sinSunAngle
#   ndvi <- (nir_ref - red_ref)/(nir_ref + red_ref)
#   
#   # Proportion of vegetation
#   pVeg <- ((ndvi - minValue(ndvi))/(maxValue(ndvi) - minValue(ndvi)))^2 # Sobrino et al. 2004: Land surface temperature retrieval from LANDSAT TM 5
#   
#   # Emissivity
#   emiss <- 0.004*pVeg + 0.986
#   
#   # Now calculate land surface temperature (LST)
#   lst_10_k <- temp10_kelvin/(1 + (11.5 * (temp10_kelvin/14380))*log(emiss)) # Weng et al. 2004: Estimation of land surface temperature-vegetation abundance relationship for urban heat island studies
#   lst_11_k <- temp11_kelvin/(1 + (11.5 * (temp11_kelvin/14380))*log(emiss))
#   
#   lst_10_c <- lst_10_k - 273.15
#   lst_11_c <- lst_11_k - 273.15
#   
#   lst_avg <- mean(lst_10_c, lst_11_c)
#   return(lst_avg)
# }
# 
# 
# #Values from Metafile
# RADIANCE_MULT_BAND_10 <- 3.3420E-04
# RADIANCE_MULT_BAND_11 <- 3.3420E-04
# 
# RADIANCE_ADD_BAND_10 <- 0.10000
# RADIANCE_ADD_BAND_11 <- 0.10000
# 
# #Values from Metafile
# K1_CONSTANT_BAND_10 <- 774.8853
# K1_CONSTANT_BAND_11 <- 480.8883
# K2_CONSTANT_BAND_10 <- 1321.0789
# K2_CONSTANT_BAND_11 <- 1201.1442
# 
# junSinSunAngle <- 0.91154
# julSinSunAngle <- 0.88129982
# augSinSunAngle <- 0.81191
# 
# # Landsat imagery can be downloaded for free from the USGS (see link above)
# jun.band_10 <- raster("DOWNLOADED RASTER IMAGE") #change image name accordingly
# jun.band_11 <- raster("DOWNLOADED RASTER IMAGE") #change image name accordingly
# jun.band_4 <- raster("DOWNLOADED RASTER IMAGE") #change image name accordingly
# jun.band_5 <- raster("DOWNLOADED RASTER IMAGE") #change image name accordingly
# 
# jul.band_10 <- raster("DOWNLOADED RASTER IMAGE") #change image name accordingly
# jul.band_11 <- raster("DOWNLOADED RASTER IMAGE") #change image name accordingly
# jul.band_4 <- raster("DOWNLOADED RASTER IMAGE") #change image name accordingly
# jul.band_5 <- raster("DOWNLOADED RASTER IMAGE") #change image name accordingly
# 
# aug.band_10 <- raster("DOWNLOADED RASTER IMAGE") #change image name accordingly
# aug.band_11 <- raster("DOWNLOADED RASTER IMAGE") #change image name accordingly
# aug.band_4 <- raster("DOWNLOADED RASTER IMAGE") #change image name accordingly
# aug.band_5 <- raster("DOWNLOADED RASTER IMAGE") #change image name accordingly
# 
# jun.lst <- calcLST(RADIANCE_MULT_BAND_10, RADIANCE_MULT_BAND_11, RADIANCE_ADD_BAND_10, RADIANCE_ADD_BAND_11, K1_B10=K1_CONSTANT_BAND_10, 
#                    K1_B11=K1_CONSTANT_BAND_11, K2_B10 = K2_CONSTANT_BAND_10, K2_B11 = K2_CONSTANT_BAND_11, junSinSunAngle, band4=jun.band_4, band5=jun.band_5, band10=jun.band_10, band11=jun.band_11)
# jul.lst <- calcLST(RADIANCE_MULT_BAND_10, RADIANCE_MULT_BAND_11, RADIANCE_ADD_BAND_10, RADIANCE_ADD_BAND_11, K1_B10=K1_CONSTANT_BAND_10, 
#                    K1_B11=K1_CONSTANT_BAND_11, K2_B10 = K2_CONSTANT_BAND_10, K2_B11 = K2_CONSTANT_BAND_11, julSinSunAngle, band4=jul.band_4, band5=jul.band_5, band10=jul.band_10, band11=jul.band_11)
# aug.lst <- calcLST(RADIANCE_MULT_BAND_10, RADIANCE_MULT_BAND_11, RADIANCE_ADD_BAND_10, RADIANCE_ADD_BAND_11, K1_B10=K1_CONSTANT_BAND_10, 
#                    K1_B11=K1_CONSTANT_BAND_11, K2_B10 = K2_CONSTANT_BAND_10, K2_B11 = K2_CONSTANT_BAND_11, augSinSunAngle, band4=aug.band_4, band5=aug.band_5, band10=aug.band_10, band11=aug.band_11)
# 
# jja.lst <- mean(jun.lst, jul.lst, aug.lst)
# 
#   
# tracts <- spTransform(tracts, CRS(proj4string(jja.lst)))
# tempTracts <- extract(jja.lst, tracts, fun=mean)
# tempTracts <- data.frame(NAME10=tracts@data$NAME10, jjaTemp=tempTracts)



###############################################
# Clean up tracts shapefile with priority scores for each benefit

priorityScoresTracts <- tracts %>%
  dplyr::select(GEOID10, stormwater01, biodiversity1_01, biodiversity2_01,
                airQuality01, wellBeing01, crime01, propValues01, heatIsland01,
                geometry) %>%
  mutate(overallPriority = stormwater01 + biodiversity1_01 +
         airQuality01 + wellBeing01 + crime01 + propValues01 + heatIsland01)

# # If you choose to load the shapefile in the data folder instead of doing the calculations above
# priorityScoresTracts <- st_read('./priorityScoresTracts.shp')

# Below is code to make a replicate of Figure 2 in Tran et al. 2020 using the calculated priority scores for Philadelphia
myTheme <- theme(legend.position = 'none',
                 plot.background = element_blank(),
                 panel.background = element_blank(),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 axis.text = element_blank(),
                 axis.ticks = element_blank(),
                 axis.line = element_blank())

stormwaterPlot <- ggplot() +
  geom_sf(data = priorityScoresTracts, aes(fill = stormwater01), col = 'white') +
  myTheme +
  labs(title = 'Stormwater')
heatIslandPlot <- ggplot() +
  geom_sf(data = priorityScoresTracts, aes(fill = heatIsland01), col = 'white') +
  myTheme +
  labs(title = 'Heat Island')
airQualityPlot <- ggplot() +
  geom_sf(data = priorityScoresTracts, aes(fill = airQuality01), col = 'white') +
  myTheme +
  labs(title = 'Air Quality')
wellBeingPlot <- ggplot() +
  geom_sf(data = priorityScoresTracts, aes(fill = wellBeing01), col = 'white') +
  myTheme +
  labs(title = 'Well-being')
multifunctionalityPlot <- ggplot() +
  geom_sf(data = priorityScoresTracts, aes(fill = overallPriority), col = 'white') +
  myTheme +
  labs(title = 'Multifunctionality')
propValuesPlot <- ggplot() +
  geom_sf(data = priorityScoresTracts, aes(fill = propValues01), col = 'white') +
  myTheme +
  labs(title = 'Property Values')
crimePlot <- ggplot() +
  geom_sf(data = priorityScoresTracts, aes(fill = crime01), col = 'white') +
  myTheme +
  labs(title = 'Crime')
biodCompPlot <- ggplot() +
  geom_sf(data = priorityScoresTracts, aes(fill = biodiversity1_01), col = 'white') +
  myTheme +
  labs(title = 'Biodiversity-\nComplementation')
biodEquityPlot <- ggplot() +
  geom_sf(data = priorityScoresTracts, aes(fill = biodiversity1_01), col = 'white') +
  myTheme +
  labs(title = 'Biodiversity-\nEquity')


library(gridExtra)
grid.arrange(stormwaterPlot, heatIslandPlot, airQualityPlot,
             wellBeingPlot, multifunctionalityPlot, propValuesPlot,
             crimePlot, biodCompPlot, biodEquityPlot, nrow = 3)




