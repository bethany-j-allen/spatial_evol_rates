#Erin Saupe   downloaded 27th November 2018
#[Code to palaeo-rotate modern coordinates into palaeocoordinates via Gplates]
#1. Code to convert modern coordinates into a shapefile for Gplates
#2. Code to convert Gplates shapefile output into palaeocoordinates

#Rotated as follows:
#Artinskian <- Map #54 "Early Permian (Artinskian, 280Ma)"
#Kungurian <- Map #53 "Early Permian (Kungurian, 275Ma)"
#Roadian & Wordian -> Map #52 "Middle Permian (Roadian & Wordian, 270Ma)"
#Capitanian -> Map #51 "Late Permian (Capitanian, 260Ma)"
#Wuchiapingian & Changhsingian -> Map #50 "Late Permian (Lopingian, 255Ma)"
#Induan -> Map #49 "Permo-Triassic boundary (250Ma)"
#Olenekian -> Map #48 "Early Triassic (Induan & Olenekian, 245Ma)"
#Anisian -> Map #47 "Middle Triassic (Anisian, 240Ma)"
#Ladinian -> Map #46 "Middle Triassic (Ladinian, 230Ma)"
#Carnian <- Map #45 "Late Triassic (Carnian, 220Ma)"
#Norian <- Map #44 "Late Triassic (Norian, 210Ma)"

#setwd("#####")

#Load packages
library(tidyverse)
library(sp)
library(rgdal)
library(rgeos)
library(maptools)


###Part 1. Convert modern lat-longs to a shape file
#Read in csv file with the present-day lat/long coordinates
dat <- read_csv("data/PT_marine_inverts.csv")


#Artinskian
#Create data frame with collection numbers, longs, lats
ArtinskianFilter <- filter(dat, stage_bin == "Artinskian")
ArtinskianPoints <- ArtinskianFilter[c("collection_no", "lng", "lat")]

#Filter to unique points (one per collection)
ArtinskianPoints <- distinct(ArtinskianPoints, collection_no, .keep_all = T)
Artinskian_s2xy <- ArtinskianPoints[c("lng", "lat")]

#Label as coordinates
coordinates(Artinskian_s2xy) <- ~ lng + lat
proj4string(Artinskian_s2xy) <- CRS("+proj=longlat +datum=WGS84")

#Convert to 'Spatial Points Data Frame'
Artinskian_shp <- SpatialPointsDataFrame(Artinskian_s2xy, ArtinskianPoints)

#Write shape file
writeOGR(Artinskian_shp, "data/Modern shapefiles", "Artinskian", driver="ESRI Shapefile")


#Kungurian
KungurianFilter <- filter(dat, stage_bin == "Kungurian")
KungurianPoints <- KungurianFilter[c("collection_no", "lng", "lat")]

KungurianPoints <- distinct(KungurianPoints, collection_no, .keep_all = T)
Kungurian_s2xy <- KungurianPoints[c("lng", "lat")]

coordinates(Kungurian_s2xy) <- ~ lng + lat
proj4string(Kungurian_s2xy) <- CRS("+proj=longlat +datum=WGS84")

Kungurian_shp <- SpatialPointsDataFrame(Kungurian_s2xy, KungurianPoints)

writeOGR(Kungurian_shp, "data/Modern shapefiles", "Kungurian", driver="ESRI Shapefile")


#Roadian & Wordian
GuadalupianFilter <- filter(dat, stage_bin == "Roadian" | stage_bin == "Wordian")
GuadalupianPoints <- GuadalupianFilter[c("collection_no", "lng", "lat")]

GuadalupianPoints <- distinct(GuadalupianPoints, collection_no, .keep_all = T)
Guadalupian_s2xy <- GuadalupianPoints[c("lng", "lat")]

coordinates(Guadalupian_s2xy) <- ~ lng + lat
proj4string(Guadalupian_s2xy) <- CRS("+proj=longlat +datum=WGS84")

Guadalupian_shp <- SpatialPointsDataFrame(Guadalupian_s2xy, GuadalupianPoints)

writeOGR(Guadalupian_shp, "data/Modern shapefiles", "Guadalupian", driver="ESRI Shapefile")


#Capitanian
CapitanianFilter <- filter(dat, stage_bin == "Capitanian")
CapitanianPoints <- CapitanianFilter[c("collection_no", "lng", "lat")]

CapitanianPoints <- distinct(CapitanianPoints, collection_no, .keep_all = T)
Capitanian_s2xy <- CapitanianPoints[c("lng", "lat")]

coordinates(Capitanian_s2xy) <- ~ lng + lat
proj4string(Capitanian_s2xy) <- CRS("+proj=longlat +datum=WGS84")

Capitanian_shp <- SpatialPointsDataFrame(Capitanian_s2xy, CapitanianPoints)

writeOGR(Capitanian_shp, "data/Modern shapefiles", "Capitanian", driver="ESRI Shapefile")


#Lopingian
LopingianFilter <- filter(dat, stage_bin == "Wuchiapingian" | stage_bin == "Changhsingian")
LopingianPoints <- LopingianFilter[c("collection_no", "lng", "lat")]

LopingianPoints <- distinct(LopingianPoints, collection_no, .keep_all = T)
Lopingian_s2xy <- LopingianPoints[c("lng", "lat")]

coordinates(Lopingian_s2xy) <- ~ lng + lat
proj4string(Lopingian_s2xy) <- CRS("+proj=longlat +datum=WGS84")

Lopingian_shp <- SpatialPointsDataFrame(Lopingian_s2xy, LopingianPoints)

writeOGR(Lopingian_shp, "data/Modern shapefiles", "Lopingian", driver="ESRI Shapefile")


#Induan
InduanFilter <- filter(dat, stage_bin == "Induan")
InduanPoints <- InduanFilter[c("collection_no", "lng", "lat")]

InduanPoints <- distinct(InduanPoints, collection_no, .keep_all = T)
Induan_s2xy <- InduanPoints[c("lng", "lat")]

coordinates(Induan_s2xy) <- ~ lng + lat
proj4string(Induan_s2xy) <- CRS("+proj=longlat +datum=WGS84")

Induan_shp <- SpatialPointsDataFrame(Induan_s2xy, InduanPoints)

writeOGR(Induan_shp, "data/Modern shapefiles", "Induan", driver="ESRI Shapefile")


#Olenekian
OlenekianFilter <- filter(dat, stage_bin == "Olenekian")
OlenekianPoints <- OlenekianFilter[c("collection_no", "lng", "lat")]

OlenekianPoints <- distinct(OlenekianPoints, collection_no, .keep_all = T)
Olenekian_s2xy <- OlenekianPoints[c("lng", "lat")]

coordinates(Olenekian_s2xy) <- ~ lng + lat
proj4string(Olenekian_s2xy) <- CRS("+proj=longlat +datum=WGS84")

Olenekian_shp <- SpatialPointsDataFrame(Olenekian_s2xy, OlenekianPoints)

writeOGR(Olenekian_shp, "data/Modern shapefiles", "Olenekian", driver="ESRI Shapefile")


#Anisian
AnisianFilter <- filter(dat, stage_bin == "Anisian")
AnisianPoints <- AnisianFilter[c("collection_no", "lng", "lat")]

AnisianPoints <- distinct(AnisianPoints, collection_no, .keep_all = T)
Anisian_s2xy <- AnisianPoints[c("lng", "lat")]

coordinates(Anisian_s2xy) <- ~ lng + lat
proj4string(Anisian_s2xy) <- CRS("+proj=longlat +datum=WGS84")

Anisian_shp <- SpatialPointsDataFrame(Anisian_s2xy, AnisianPoints)

writeOGR(Anisian_shp, "data/Modern shapefiles", "Anisian", driver="ESRI Shapefile")


#Ladinian
LadinianFilter <- filter(dat, stage_bin == "Ladinian")
LadinianPoints <- LadinianFilter[c("collection_no", "lng", "lat")]

LadinianPoints <- distinct(LadinianPoints, collection_no, .keep_all = T)
Ladinian_s2xy <- LadinianPoints[c("lng", "lat")]

coordinates(Ladinian_s2xy) <- ~ lng + lat
proj4string(Ladinian_s2xy) <- CRS("+proj=longlat +datum=WGS84")

Ladinian_shp <- SpatialPointsDataFrame(Ladinian_s2xy, LadinianPoints)

writeOGR(Ladinian_shp, "data/Modern shapefiles", "Ladinian", driver="ESRI Shapefile")


#Carnian
CarnianFilter <- filter(dat, stage_bin == "Carnian")
CarnianPoints <- CarnianFilter[c("collection_no", "lng", "lat")]

CarnianPoints <- distinct(CarnianPoints, collection_no, .keep_all = T)
Carnian_s2xy <- CarnianPoints[c("lng", "lat")]

coordinates(Carnian_s2xy) <- ~ lng + lat
proj4string(Carnian_s2xy) <- CRS("+proj=longlat +datum=WGS84")

Carnian_shp <- SpatialPointsDataFrame(Carnian_s2xy, CarnianPoints)

writeOGR(Carnian_shp, "data/Modern shapefiles", "Carnian", driver="ESRI Shapefile")


#Norian
NorianFilter <- filter(dat, stage_bin == "Norian")
NorianPoints <- NorianFilter[c("collection_no", "lng", "lat")]

NorianPoints <- distinct(NorianPoints, collection_no, .keep_all = T)
Norian_s2xy <- NorianPoints[c("lng", "lat")]

coordinates(Norian_s2xy) <- ~ lng + lat
proj4string(Norian_s2xy) <- CRS("+proj=longlat +datum=WGS84")

Norian_shp <- SpatialPointsDataFrame(Norian_s2xy, NorianPoints)

writeOGR(Norian_shp, "data/Modern shapefiles", "Norian", driver="ESRI Shapefile")


###Part 2. Convert shape file back to palaeo lat-longs
#Read shape file, convert to data frame, write as csv
dat <- read_csv("data/PT_brach_biv.csv")


#Artinskian
Artinskian_shp <- readOGR("data/Rotated shapefiles/Artinskian/reconstructed_280.00Ma.shp")
Artinskian_dat <- data.frame(Artinskian_shp)
Artinskian_s2xy <- dplyr::select(Artinskian_dat, collection_no = cllctn_, coords.x1, coords.x2)


#Kungurian
Kungurian_shp <- readOGR("data/Rotated shapefiles/Kungurian/reconstructed_275.00Ma.shp")
Kungurian_dat <- data.frame(Kungurian_shp)
Kungurian_s2xy <- dplyr::select(Kungurian_dat, collection_no = cllctn_, coords.x1, coords.x2)


#Roadian & Wordian
Guadalupian_shp <- readOGR("data/Rotated shapefiles/Guadalupian/reconstructed_270.00Ma.shp")
Guadalupian_dat <- data.frame(Guadalupian_shp)
Guadalupian_s2xy <- dplyr::select(Guadalupian_dat, collection_no = cllctn_, coords.x1, coords.x2)


#Capitanian
Capitanian_shp <- readOGR("data/Rotated shapefiles/Capitanian/reconstructed_260.00Ma.shp")
Capitanian_dat <- data.frame(Capitanian_shp)
Capitanian_s2xy <- dplyr::select(Capitanian_dat, collection_no = cllctn_, coords.x1, coords.x2)


#Lopingian
Lopingian_shp <- readOGR("data/Rotated shapefiles/Lopingian/reconstructed_255.00Ma.shp")
Lopingian_dat <- data.frame(Lopingian_shp)
Lopingian_s2xy <- dplyr::select(Lopingian_dat, collection_no = cllctn_, coords.x1, coords.x2)


#Induan
Induan_shp <- readOGR("data/Rotated shapefiles/Induan/reconstructed_250.00Ma.shp")
Induan_dat <- data.frame(Induan_shp)
Induan_s2xy <- dplyr::select(Induan_dat, collection_no = cllctn_, coords.x1, coords.x2)


#Olenekian
Olenekian_shp <- readOGR("data/Rotated shapefiles/Olenekian/reconstructed_245.00Ma.shp")
Olenekian_dat <- data.frame(Olenekian_shp)
Olenekian_s2xy <- dplyr::select(Olenekian_dat, collection_no = cllctn_, coords.x1, coords.x2)


#Anisian
Anisian_shp <- readOGR("data/Rotated shapefiles/Anisian/reconstructed_240.00Ma.shp")
Anisian_dat <- data.frame(Anisian_shp)
Anisian_s2xy <- dplyr::select(Anisian_dat, collection_no = cllctn_, coords.x1, coords.x2)


#Ladinian
Ladinian_shp <- readOGR("data/Rotated shapefiles/Ladinian/reconstructed_230.00Ma.shp")
Ladinian_dat <- data.frame(Ladinian_shp)
Ladinian_s2xy <- dplyr::select(Ladinian_dat, collection_no = cllctn_, coords.x1, coords.x2)


#Carnian
Carnian_shp <- readOGR("data/Rotated shapefiles/Carnian/reconstructed_220.00Ma.shp")
Carnian_dat <- data.frame(Carnian_shp)
Carnian_s2xy <- dplyr::select(Carnian_dat, collection_no = cllctn_, coords.x1, coords.x2)


#Norian
Norian_shp <- readOGR("data/Rotated shapefiles/Norian/reconstructed_210.00Ma.shp")
Norian_dat <- data.frame(Norian_shp)
Norian_s2xy <- dplyr::select(Norian_dat, collection_no = cllctn_, coords.x1, coords.x2)


#Compile
AllPoints <- rbind(Artinskian_s2xy, Kungurian_s2xy, Guadalupian_s2xy, Capitanian_s2xy, Lopingian_s2xy,
                   Induan_s2xy, Olenekian_s2xy, Anisian_s2xy, Ladinian_s2xy, Carnian_s2xy, Norian_s2xy)
AllPoints <- rename(AllPoints, rot_lng = coords.x1, rot_lat = coords.x2)

#Write a csv of points
write.csv(AllPoints, "data/rotated_points.csv")

#Match points to collection numbers
new_db <- left_join(dat, AllPoints, by = "collection_no")

#Write a new csv with the paleolat & long added
write_csv(new_db, path = "data/PT_brach_biv.csv")
