# Author: Maria A. Hurtado-Materon
# Date: May 1, 2025
# Description: This script automates IUCN species reassessments by generating
# updated range polygons, cleaning occurrence data, and calculating key metrics.
# This script:
# - Loads and cleans occurrence data
# - Generates a new range polygon
# - Saves an IUCN shapefile
# - Produces an interactive map
# - Calculates EOO and AOO

# ========================
# Load necessary packages
# ========================
library(rgbif)
library(leafem)
library(leaflet.esri)
library(plyr)
library(htmlwidgets)
library(sf)
library(raster)
library(dplyr)
library(ggplot2)
library(mapedit)
library(mapview)
library(redlistr)
library(leaflet)
library(leaflet.extras)
library(tidyverse)
library(smoothr)
library(red)
library(devtools)

# Install the SMSGRedList package from GitHub
devtools::install_github("mariahm1995/SMSGRedList", force = TRUE)
library(SMSGRedList)

# If using within a package, you can skip this and use the package namespace
# source("scripts/functions.R")  # <- only needed if not loading as a package

# ========================
# Configuration
# ========================
# Download the files needed for the assessment (only one time)
options(timeout = 600)
download.file("https://ndownloader.figshare.com/files/54137138", destfile = "data.zip", mode = "wb")
unzip("data.zip")
unlink("data.zip")

# ========================
# Species selection
# ========================
spname <- c("Chinchilla_lanigera")
scientific_name <- c("Chinchilla lanigera")

# ========================
# Load previous IUCN polygon
# ========================
previous_map <- st_read("data/mammal_IUCN_2024/MAMMALS.shp") %>%
  filter(sci_name == scientific_name)

plot(previous_map$geometry)

# ========================
# Country overlap
# ========================
countries_list <- countries_intersection(map = previous_map)
countries_list

# ========================
# Occurrence data
# ========================

# GBIF data
occsGBIF <- sRL_SimpleGBIF(co_EXT = c(-170,-50,-90,90))

# Literature records
occsLit <- read.csv(paste0("species/", spname, "/lit_", spname, "_Occurrences.csv")) %>%
  convert_to_decimal() %>%
  dplyr::select(c("decimalLatitude", "decimalLongitude",
                  "coordinateUncertaintyInMeters", "year", "month", "day",
                  "country", "occurrenceID", "Source_type", "Link", "gbifID"))

# ========================
# Combine and clean occurrences
# ========================
occs <- combineOccs(occsGBIF = occsGBIF, occsLit = occsLit, combine = FALSE)
occs <- extractElevation(occs = occs)

occs_cleaned <- sRL_cleanDataGBIF(
  flags = occs,
  year_GBIF = 1800,
  uncertainty_GBIF = 0.5,
  GBIF_xmin = -180,
  GBIF_xmax = 180,
  GBIF_ymin = -90,
  GBIF_ymax = 90,
  AltMIN = 0,
  AltMAX = 700,
  land = TRUE
)

occs_cleaned <- sRL_PopRecords(occs_cleaned)
occs_filter <- occs_cleaned %>% filter(is.na(Reason))

# ========================
# Quick map to check points
# ========================
ggplot() +
  geom_sf(data = previous_map) +
  geom_point(data = occs_cleaned,
             aes(decimalLongitude, decimalLatitude, color = "blue")) +
  geom_point(data = occs_filter,
             aes(decimalLongitude, decimalLatitude, color = "red"))

# ========================
# Map generation from occurrences
# ========================
occs_filter_spatial <- st_as_sf(
  occs_filter[c("decimalLongitude", "decimalLatitude")],
  coords = c("decimalLongitude", "decimalLatitude"),
  crs = crs(previous_map)
)

modification <- sRL_MapDistributionGBIF(
  dat = occs_filter_spatial,
  First_step = "mcp",
  Buffer_km = 15,
  GBIF_crop = "cropsea",
  Gbif_Param = NULL
)

processed_maps <- checkMap(
  map = modification$geometry,
  method = "ksmooth",
  smoothness = 20,
  smooth = TRUE,
  elevation_range = TRUE,
  AltMIN = 0,
  AltMAX = 800,
  keep = 0.05,
  continental = TRUE
)

# ========================
# Manual editing (if needed)
# ========================
new_area <- manualEdition(
  modification_type = "add",
  maps = processed_maps,
  previous_map = previous_map,
  occs = occs_filter_spatial
)

# ========================
# Export IUCN shapefile
# ========================
final_map <- exportSpatialobject(
  previous_map = previous_map,
  new_map = processed_maps$mcc_occs,
  currentyear = 2025
)

final_map <- st_read(paste0("species/", spname, "/final_map/final_map.shp"))

# ========================
# Create interactive map
# ========================
leaflet_map <- sRL_LeafletFlags(
  flags = occs_cleaned,
  previous_map = previous_map,
  final_map = final_map,
  elevation_map = processed_maps$crop_alt_simplified,
  occs_filter_spatial = occs_filter_spatial
)

# ========================
# EOO and AOO calculations
# ========================
eoo_km2 <- EOOcalculation(map = previous_map)
aoo_points <- aooCalculation(spatial = occs_filter_spatial, type = "occs")
aoo_poly <- aooCalculation(spatial = final_map, type = "poly")

