
#########################################################################################################
#
# Code to accompany "Green Infrastructure Space and Traits (GIST) Model: Integrating green infrastructure 
# spatial placement and plant traits to maximize multifunctionality"
# Urban Forestry & Urban Greening
#
# Tyler J. Tran, Matthew R. Helmus, and Jocelyn E. Behm
# 2020
#
# Note: This code does not apply to all of the analyses in the Tran et al. 2020 paper. Because of our 
# data sharing agreement with partners, we are not allowed to share some code. The code below concerns
# calculating condition scores and priority scores for the spatial component of the GIST model. The 
# outputs of each code chunk are aggregates for each benefit at the census tract level. To calculate
# condition scores and priority scores, these aggregates were scaled from 0 to 1 to be relative to each
# other.
#
#########################################################################################################


# You can download the data from the iEcoLab github. If you put it all in one folder and point to the 
# folder's directory here, the data should load without issue. Some of the larger files are not provided
# on our github because of file size. In those cases, links to download are given below.
setwd('')


library(tidyverse); library(sf); library(lwgeom); library(raster); library(rgdal)


# Load census tract shapefile
tracts <- st_read('./Census_Tracts_2010.shp')
tracts <- tracts %>%
  mutate(NAME10 = as.character(NAME10))


# Below is land cover data for each census tract, with sums of # pixels of each class for each tract
# Column percentImp is (buildings + roads + otherPaved)/totalPixels
impervious <- read_csv('./landCover_tract.csv') %>%
  mutate(tract = as.character(tract))


# Join census tract shapefile to impervious surface data
tracts <- left_join(tracts, impervious, by = c('NAME10' = 'tract'))

tracts <- tracts %>%
  mutate(tractArea = as.numeric(st_area(tracts)))




###########################################################
# Biodiversity
# First show using sf methods, then show using sp methods

greenSpace <- st_read('./greenspace.shp') %>%
  mutate(valid = st_is_valid(greenSpace))

corruptPolys <- greenSpace %>%
  filter(valid == FALSE) %>%
  st_make_valid()

greenSpace <- greenSpace %>%
  filter(valid == TRUE) %>%
  rbind(corruptPolys) %>%
  st_transform(4326)

greenSpace <- greenSpace %>%
  mutate(area = as.numeric(st_area(greenSpace)))

greenspace_centroids <- st_centroid(greenSpace)

tractGreen <- st_intersection(greenspace_centroids, tracts) %>%
  group_by(GEOID10) %>%
  summarise(areaGreenspace = sum(area, na.rm = T),
            tractArea = unique(tractArea),
            nPatches = n()) %>%
  mutate(meanPatchSize = areaGreenspace/nPatches)

# now sp methods
# vectorized version of raster land cover data of green space (vegetated land cover types)
green <- readOGR(dsn='.', layer='greenspace')
green <- spTransform(green, CRS(proj4string(tracts)))
green@data$area <- areaPolygon(green, a=6378137, f=1/298.257223563)
green.cent <- getSpPPolygonsLabptSlots(green)
green.cent <- as.data.frame(green.cent)
colnames(green.cent) <- c('long','lat')
coordinates(green.cent) <- ~long+lat
proj4string(green.cent) <- proj4string(green)
tractGreen <- over(green.cent, tracts)
tractGreen <- data.frame(NAME10=tractGreen$NAME10, area=green@data$area)
tractGreen <- group_by(tractGreen, NAME10)
tractGreen <- summarise(tractGreen, nGreen=n(), area=sum(area))
tractGreen$meanPatchSize <- tractGreen$area/tractGreen$nGI
tractGreen$NAME10 <- as.character(tractGreen$NAME10)
tractGreen$NAME10 <- as.numeric(tractGreen$NAME10)
tractGreen <- data.frame(NAME10=tractGreen$NAME10, meanPatchSize=tractGreen$meanPatchSize)


################################################################
# Crime
# Public crime incident data, see https://www.opendataphilly.org/dataset/crime-incidents 

crime <- st_read('./crime2016tracts.shp') %>%
  group_by(GEOID10) %>%
  summarise(nCrimes = n())


################################################################
# Air quality improvements

# Because this is a large file, you should download.
# See https://www.epa.gov/air-research/downscaler-model-predicting-daily-air-pollution

pm25 <- read_csv('./2013_pm25_daily_average.txt') %>%
  rename(dailyAvgPM25 = `pm25_daily_average(ug/m3)`) %>%
  mutate(FIPS = as.character(FIPS),
         dailyAvgPM25 = as.numeric(dailyAvgPM25)) %>%
  filter(Date >= lubridate::ymd('2013-06-01') & Date <= lubridate::ymd('2013-06-30'),
         startsWith(FIPS, prefix = '42101')) %>% # this was national data, narrow down to Philadelphia County
  group_by(FIPS) %>%
  summarise(meanPM25 = mean(dailyAvgPM25, na.rm = T),
            medianPM25 = median(dailyAvgPM25, na.rm = T))



################################################################
# Improved human well-being
# See https://www.cdc.gov/500cities/index.htm

health <- read_csv('./500cities.csv') %>%
  filter(StateAbbr == 'PA' & PlaceName == 'Philadelphia') %>%
  mutate(mentalHealth = MHLTH_CrudePrev/Population2010, physHealth = PHLTH_CrudePrev/Population2010) %>%
  dplyr::select(TractFIPS, mentalHealth, physHealth)


################################################################
# Property values: Data from Philadelphia public OPA data
# https://www.opendataphilly.org/dataset/opa-property-assessments

propVal <- read_csv('./propValuesTracts.csv') %>%
  filter(category_1 %in% c('Single Family', 'Multi Family')) %>%
  mutate(pricePerSqft = market_val/total_liva) %>%
  group_by(GEOID10) %>%
  summarise(meanPropVal = mean(market_val, na.rm = T),
            medianPropVal = median(market_val, na.rm = T))


################################################################
# Heat island mediation
#
# Much of the code below is adapted from this blog post: https://www.gis-blog.com/calculation-of-land-surface-temperature-lst-from-landsat-8-using-r/
# We have not included large file-size Landsat imagery in our included data. However, the data can be downloaded for free at https://www.usgs.gov/land-resources/nli/landsat/landsat-data-access?qt-science_support_page_related_con=0#qt-science_support_page_related_con

# This function will be used below to pre-process the imagery
calcLST <- function(RADIANCE_MULT_BAND_10, RADIANCE_MULT_BAND_11, RADIANCE_ADD_BAND_10, RADIANCE_ADD_BAND_11, K1_B10, K1_B11, K2_B10, K2_B11, sinSunAngle,
                    band4, band5, band10, band11, censusTracts){
  
  tract <- spTransform(censusTracts, CRS(proj4string(band_10)))
  band_10 <- crop(band10, tract)
  band_11 <- crop(band11, tract)
  band_4 <- crop(band4, tract)
  band_5 <- crop(band5, tract)
  
  #Calculate TOA from DN:
  toa_band10 <- calc(band_10, fun=function(x){RADIANCE_MULT_BAND_10 * x + RADIANCE_ADD_BAND_10})
  toa_band11 <- calc(band_11, fun=function(x){RADIANCE_MULT_BAND_11 * x + RADIANCE_ADD_BAND_11})
  
  
  #Calculate LST in Kelvin for Band 10 and Band 11
  temp10_kelvin <- calc(toa_band10, fun=function(x){K2_B10/log(K1_B10/x + 1)})
  temp11_kelvin <- calc(toa_band11, fun=function(x){K2_B11/log(K1_B11/x + 1)})
  
  #Convert Kelvin to Celsius for Band 10 and 11
  temp10_celsius <- calc(temp10_kelvin, fun=function(x){x - 273.15})
  temp11_celsius <- calc(temp11_kelvin, fun=function(x){x - 273.15})
  
  atSat_kelvin <- mean(temp10_kelvin, temp11_kelvin)
  atSat_celsius <- mean(temp10_celsius, temp11_celsius)
  
  # Calculate NDVI
  red_ref <- ((0.00002*band_4) - 0.1)/sinSunAngle
  nir_ref <- ((0.00002*band_5) - 0.1)/sinSunAngle
  ndvi <- (nir_ref - red_ref)/(nir_ref + red_ref)
  
  # Proportion of vegetation
  pVeg <- ((ndvi - minValue(ndvi))/(maxValue(ndvi) - minValue(ndvi)))^2 # Sobrino et al. 2004: Land surface temperature retrieval from LANDSAT TM 5
  
  # Emissivity
  emiss <- 0.004*pVeg + 0.986
  
  # Now calculate land surface temperature (LST)
  lst_10_k <- temp10_kelvin/(1 + (11.5 * (temp10_kelvin/14380))*log(emiss)) # Weng et al. 2004: Estimation of land surface temperature-vegetation abundance relationship for urban heat island studies
  lst_11_k <- temp11_kelvin/(1 + (11.5 * (temp11_kelvin/14380))*log(emiss))
  
  lst_10_c <- lst_10_k - 273.15
  lst_11_c <- lst_11_k - 273.15
  
  lst_avg <- mean(lst_10_c, lst_11_c)
  return(lst_avg)
}


#Values from Metafile
RADIANCE_MULT_BAND_10 <- 3.3420E-04
RADIANCE_MULT_BAND_11 <- 3.3420E-04

RADIANCE_ADD_BAND_10 <- 0.10000
RADIANCE_ADD_BAND_11 <- 0.10000

#Values from Metafile
K1_CONSTANT_BAND_10 <- 774.8853
K1_CONSTANT_BAND_11 <- 480.8883
K2_CONSTANT_BAND_10 <- 1321.0789
K2_CONSTANT_BAND_11 <- 1201.1442

junSinSunAngle <- 0.91154
julSinSunAngle <- 0.88129982
augSinSunAngle <- 0.81191

# Landsat imagery can be downloaded for free from the USGS (see link above)
jun.band_10 <- raster("DOWNLOADED RASTER IMAGE") #change image name accordingly
jun.band_11 <- raster("DOWNLOADED RASTER IMAGE") #change image name accordingly
jun.band_4 <- raster("DOWNLOADED RASTER IMAGE") #change image name accordingly
jun.band_5 <- raster("DOWNLOADED RASTER IMAGE") #change image name accordingly

jul.band_10 <- raster("DOWNLOADED RASTER IMAGE") #change image name accordingly
jul.band_11 <- raster("DOWNLOADED RASTER IMAGE") #change image name accordingly
jul.band_4 <- raster("DOWNLOADED RASTER IMAGE") #change image name accordingly
jul.band_5 <- raster("DOWNLOADED RASTER IMAGE") #change image name accordingly

aug.band_10 <- raster("DOWNLOADED RASTER IMAGE") #change image name accordingly
aug.band_11 <- raster("DOWNLOADED RASTER IMAGE") #change image name accordingly
aug.band_4 <- raster("DOWNLOADED RASTER IMAGE") #change image name accordingly
aug.band_5 <- raster("DOWNLOADED RASTER IMAGE") #change image name accordingly

jun.lst <- calcLST(RADIANCE_MULT_BAND_10, RADIANCE_MULT_BAND_11, RADIANCE_ADD_BAND_10, RADIANCE_ADD_BAND_11, K1_B10=K1_CONSTANT_BAND_10, 
                   K1_B11=K1_CONSTANT_BAND_11, K2_B10 = K2_CONSTANT_BAND_10, K2_B11 = K2_CONSTANT_BAND_11, junSinSunAngle, band4=jun.band_4, band5=jun.band_5, band10=jun.band_10, band11=jun.band_11)
jul.lst <- calcLST(RADIANCE_MULT_BAND_10, RADIANCE_MULT_BAND_11, RADIANCE_ADD_BAND_10, RADIANCE_ADD_BAND_11, K1_B10=K1_CONSTANT_BAND_10, 
                   K1_B11=K1_CONSTANT_BAND_11, K2_B10 = K2_CONSTANT_BAND_10, K2_B11 = K2_CONSTANT_BAND_11, julSinSunAngle, band4=jul.band_4, band5=jul.band_5, band10=jul.band_10, band11=jul.band_11)
aug.lst <- calcLST(RADIANCE_MULT_BAND_10, RADIANCE_MULT_BAND_11, RADIANCE_ADD_BAND_10, RADIANCE_ADD_BAND_11, K1_B10=K1_CONSTANT_BAND_10, 
                   K1_B11=K1_CONSTANT_BAND_11, K2_B10 = K2_CONSTANT_BAND_10, K2_B11 = K2_CONSTANT_BAND_11, augSinSunAngle, band4=aug.band_4, band5=aug.band_5, band10=aug.band_10, band11=aug.band_11)

jja.lst <- mean(jun.lst, jul.lst, aug.lst)

  
tracts <- spTransform(tracts, CRS(proj4string(jja.lst)))
tempTracts <- extract(jja.lst, tracts, fun=mean)
tempTracts <- data.frame(NAME10=tracts@data$NAME10, jjaTemp=tempTracts)




