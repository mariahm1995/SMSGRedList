#' Extract country names from overlapping map polygons
#'
#' Identifies which country polygons intersect with a given species distribution map.
#'
#' @param map An `sf` object of the species' spatial range.
#' @return A character vector of country names that intersect the input map.
#' @export
countries_intersection <- function(map){
  countries_pol <- st_read("data/geodata/polygons/Countries/") %>% st_make_valid()
  map <- st_make_valid(map)
  overlap_countries <- st_intersection(countries_pol, map)
  countries_list <- overlap_countries$ADM0_NAME
  return(countries_list)
}

#' Convert degrees, minutes, and seconds to decimal coordinates
#'
#' Converts occurrence coordinates stored in DMS format to decimal degrees.
#'
#' @param data A data frame with columns: `degreeLat`, `minLat`, `secLat`, `degreeLon`, `minLon`, `secLon`, `cardinalDirectionLat`, `cardinalDirectionLon`.
#' @return A data frame with added `decimalLatitude` and `decimalLongitude` columns.
#' @export
convert_to_decimal <- function(data) {
  for (i in seq_len(nrow(data))) {
    if (is.na(data$degreeLat[i]) || is.na(data$minLat[i]) ||
        is.na(data$degreeLon[i]) || is.na(data$minLon[i])) {
      message("Skipping row ", i, ": missing coordinate components")
      next
    }

    secLat <- suppressWarnings(as.numeric(data$secLat[i]))
    secLon <- suppressWarnings(as.numeric(data$secLon[i]))
    if (is.na(secLat)) secLat <- 0
    if (is.na(secLon)) secLon <- 0

    data$decimalLatitude[i] <- data$degreeLat[i] + data$minLat[i]/60 + secLat/3600
    data$decimalLongitude[i] <- data$degreeLon[i] + data$minLon[i]/60 + secLon/3600

    if (data$cardinalDirectionLat[i] == "S") data$decimalLatitude[i] <- -data$decimalLatitude[i]
    if (data$cardinalDirectionLon[i] == "W") data$decimalLongitude[i] <- -data$decimalLongitude[i]
  }
  return(data)
}

#' Combine GBIF and literature-based occurrences
#'
#' Merges occurrence datasets, ensuring column consistency and optional source control.
#'
#' @param occsGBIF A data frame of GBIF records.
#' @param occsLit A data frame of literature-sourced records.
#' @param combine Logical. If TRUE (default), merges both datasets. If FALSE, only GBIF is used.
#' @return A unified occurrence data frame with harmonized columns.
#' @export
combineOccs <- function(occsGBIF, occsLit, combine = TRUE){
  if (!combine) {
    occs <- occsGBIF
    occs$Link <- NA
  } else {
    occsLit[setdiff(names(occsGBIF), names(occsLit))] <- NA
    occsGBIF[setdiff(names(occsLit), names(occsGBIF))] <- NA
    occs <- rbind(occsLit, occsGBIF)
  }
  if (!"Link" %in% colnames(occs)) occs$Link <- NA
  return(occs)
}

#' Extract elevation for each occurrence record
#'
#' Adds elevation values (in meters) to each record based on a global DEM raster.
#'
#' @param occs A data frame with decimal coordinates and occurrence information.
#' @return The same data frame with a new column `elevationCalculated`.
#' @export
extractElevation <- function(occs){
  raster <- raster::raster("data/geodata/rasters/DEM_globe_raster/dem_globe1.tif")
  occsSpatial <- st_as_sf(occs[c("decimalLongitude", "decimalLatitude")],
                          coords = c("decimalLongitude", "decimalLatitude"), crs = crs(previous_map))
  occs$elevationCalculated <- terra::extract(raster, occsSpatial)
  return(occs)
}

#' Export a spatial object for IUCN Red List submission
#'
#' Combines metadata and new geometry into an IUCN-compatible shapefile.
#'
#' @param previous_map `sf` object from previous IUCN assessment.
#' @param new_map New `sf` geometry (processed or manually edited).
#' @param currentyear Numeric. The current assessment year.
#' @return The exported geometry object.
#' @export
exportSpatialobject <- function(previous_map, new_map, currentyear){
  template <- st_read("data/geodata/empty_polygon/")
  final_map <- template
  final_map <- plyr::rbind.fill(final_map, data.frame(OBJECTID = NA))

  empty_matrix <- as.data.frame(matrix(data = NA, nrow = length(new_map)-1, ncol = ncol(final_map)))
  colnames(empty_matrix) <- colnames(final_map)
  final_map <- rbind(final_map, empty_matrix)

  final_map$OBJECTID <- seq(1, nrow(final_map), 1)
  for (i in seq_along(new_map)) {
    final_map$geometry[i] <- new_map[i]
  }

  final_map$id_no <- previous_map$id_no
  final_map$sci_name <- previous_map$sci_name
  final_map$presence <- previous_map$presence
  final_map$origin <- previous_map$origin
  final_map$seasonal <- previous_map$seasonal
  final_map$compiler <- previous_map$compiler
  final_map$yrcompiled <- currentyear
  final_map$citation <- previous_map$citation

  shape_areas <- round(as.numeric(st_area(final_map$geometry))/1000, 0)
  final_map$Shape_Area <- ifelse(any(shape_areas > 1000000), NA, shape_areas)
  final_map <- final_map[, colSums(is.na(final_map)) < nrow(final_map)]

  final_map <- st_as_sf(final_map)

  dir_path <- paste0("species/", spname, "/", spname, "_final_map")
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
    message("The folder has been created.")
  } else {
    message("The folder already exists.")
  }

  st_write(final_map, paste0(dir_path, "/", spname, "_final_map.shp"), append = FALSE)
  return(final_map$geometry)
}

#' Clean and simplify distribution map geometry
#'
#' Smooths and filters a map using optional elevation and land boundaries.
#'
#' @param map An `sf` geometry object to be processed.
#' @param method Character. Smoothing method: `"ksmooth"`, `"chaikin"`, or `"spline"`.
#' @param smoothness Numeric. Degree of smoothing (only for `"ksmooth"`).
#' @param smooth Logical. Whether to apply smoothing.
#' @param elevation_range Logical. Whether to clip to elevation.
#' @param AltMIN,AltMAX Minimum and maximum elevation thresholds.
#' @param keep Numeric between 0 and 1. Proportion of detail to keep in simplification.
#' @param continental Logical. If TRUE, limits to continental land only.
#'
#' @return A named list of cleaned maps: \code{mcc_occs} and \code{crop_alt_simplified}.
#' @export
checkMap <- function(map, method, smoothness = NULL, smooth, elevation_range,
                     AltMIN, AltMAX, keep = .05, continental) {

  elevation <- raster::raster("data/geodata/rasters/DEM_globe_raster/dem_globe1.tif")
  land_full <- st_read("data/geodata/polygons/Land_Masses_and_Ocean_Islands/") %>%
    st_make_valid()

  if (continental == TRUE){
    land_full <- land_full %>% slice(1:7)
  }

  sf::sf_use_s2(FALSE)
  land <- st_crop(land_full, st_bbox(map) + c(-2, -2, 2, 2)) %>% st_union()
  sf::sf_use_s2(TRUE)

  if (smooth == TRUE) {
    if (method == "ksmooth") {
      map_smoothed <- smoothr::smooth(map, method = method, smoothness = smoothness)
      new <- st_intersection(st_as_sf(map_smoothed), land)
    }
    if (method %in% c("chaikin", "spline")) {
      map_smoothed <- smoothr::smooth(map, method = method)
      new <- st_intersection(st_as_sf(map_smoothed), land)
    }
  } else {
    new <- st_intersection(st_as_sf(map), land)
  }

  if (elevation_range == TRUE){
    m <- c(AltMIN, AltMAX, 1,
           AltMAX, Inf, NA,
           -Inf, AltMIN, NA)
    rclmat <- matrix(m, ncol = 3, byrow = TRUE)

    elevation_crop <- crop(elevation, st_bbox(map))
    raster_masked <- terra::classify(terra::rast(elevation_crop), rclmat)
    raster_polygons <- terra::as.polygons(raster_masked, na.rm = TRUE)
    raster_polygons_sf <- st_as_sf(raster_polygons)

    polygon_intersect <- st_intersection(map, raster_polygons_sf)
    simplified_polygon <- st_union(polygon_intersect)

    oversimplified <- rmapshaper::ms_simplify(simplified_polygon, keep = keep,
                                              keep_shapes = TRUE)
  } else {
    oversimplified <- NULL
  }

  list_maps <- list(new$x, oversimplified)
  names(list_maps) <- c("mcc_occs", "crop_alt_simplified")
  return(list_maps)
}

#' Manually edit a distribution map via interactive interface
#'
#' Allows users to add, remove, or redraw polygon areas using `mapview::editMap()`.
#'
#' @param modification_type Character. `"add"`, `"remove"`, or `"keep_as_drawn"`.
#' @param maps A list output from `checkMap()`, including two geometry layers.
#' @param previous_map The existing `sf` map to edit.
#' @param occs An `sf` object of valid occurrence points for reference.
#'
#' @return An `sf` object representing the manually modified geometry.
#' @export
manualEdition <- function(modification_type, maps, previous_map, occs){
  sf::sf_use_s2(FALSE)
  on.exit(sf::sf_use_s2(TRUE))

  map1 <- previous_map
  map1$geometry <- st_make_valid(st_union(maps[[1]]))

  elevation_map <- previous_map
  elevation_map$geometry <- st_make_valid(st_union(maps[[2]]))

  maps_all <- list(map1, elevation_map, previous_map)
  names(maps_all) <- c("new_range", "elevation_limits", "previous_map")

  modification <- (mapview(maps_all, col.regions = c("#335c67", "#e09f3e", "#9e2a2b")) + mapview(occs)) %>% editMap()

  if (modification_type == "add") {
    new <- st_union(modification$drawn$geometry, previous_map$geometry)
  } else if (modification_type == "remove") {
    new <- st_difference(previous_map$geometry, modification$drawn$geometry)
  } else if (modification_type == "keep_as_drawn") {
    new <- st_sf(geometry = modification$drawn$geometry)
  }

  land <- st_read("data/geodata/polygons/Land_Masses_and_Ocean_Islands/") %>% st_make_valid()
  map_smoothed <- smoothr::smooth(new, method = "ksmooth", smoothness = 0.5)
  new_cropped <- st_intersection(map_smoothed, land)

  return(new_cropped)
}

#' Calculate Extent of Occurrence (EOO)
#'
#' Uses a minimum convex polygon (MCP) to compute EOO area in km².
#'
#' @param map An `sf` object representing a species' range.
#'
#' @return Numeric. EOO in km².
#' @export
EOOcalculation <- function(map){
  sf::sf_use_s2(FALSE)
  on.exit(sf::sf_use_s2(TRUE))

  EOO <- st_as_sf(st_convex_hull(st_union(map)))
  EOO_km2 <- round(as.numeric(st_area(st_transform(EOO, st_crs(4326)))) / 1000000) %>% max(c(., 4), na.rm = TRUE)

  ggsave(paste0("species/", spname, "/plot_eoo.png"),
         ggplot() +
           geom_sf(data = map, fill = "#fcbba1", col = NA) +
           geom_sf(data = EOO, col = "#ef3b2c", fill = NA, lwd = 2) +
           ggtitle("EOO map") +
           theme_minimal(),
         width = 6, height = 6)
  return(EOO_km2)
}

#' Calculate Area of Occupancy (AOO)
#'
#' Computes the Area of Occupancy (AOO) using either point-based or polygon-based methods,
#' with 2x2 km grid cells as defined by IUCN.
#'
#' @param spatial An `sf` object. Either cleaned occurrence points (type = "occs") or species range polygon (type = "poly").
#' @param type Character. Method to use: `"occs"` for point-based, `"poly"` for polygon-based.
#'
#' @return Numeric. AOO in km².
#' @export
aooCalculation <- function(spatial, type){
  if (type == "occs") {
    grid22 <- raster::raster("data/geodata/sRedList/Empty.grid.2x2.Mollweide.tif")
    occs_proj <- st_transform(spatial, crs(grid22))
    grid_crop <- crop(grid22, (extent(occs_proj) + c(-100000, 100000, -100000, 100000)), snap = "out")

    pts <- occs_proj %>% as_Spatial() %>% as(., 'SpatialPoints')
    AOO_pts <- terra::rasterize(pts, grid_crop, fun = 'count') >= 1
    aooValue <- sum(as.vector(AOO_pts), na.rm = TRUE) * 4
  }

  if (type == "poly") {
    multipolygon <- spatial %>%
      dplyr::group_by(sci_name) %>%
      dplyr::summarise(N = dim(spatial)[1])

    spat_raster <- rasterize(multipolygon,
                             crop(raster::raster("data/geodata/rasters/DEM_globe_raster/dem_globe1.tif"),
                                  st_bbox(spatial))) >= 1
    spat_raster_SR <- terra::rast(spat_raster)
    aooValue <- red::aoo(spat_raster_SR)
  }

  return(aooValue)
}

#' Download GBIF occurrences for a species
#'
#' A simplified wrapper for downloading GBIF records within custom coordinate bounds.
#'
#' @param co_EXT A numeric vector of length 4 indicating longitude and latitude bounds (xmin, xmax, ymin, ymax).
#'
#' @return A data frame of GBIF records.
#' @export
sRL_SimpleGBIF <- function(co_EXT){
  dat_gbif <- rgbif::occ_data(
    scientificName = scientific_name,
    hasCoordinate = TRUE,
    decimalLongitude = paste(co_EXT[1], co_EXT[2], sep = ","),
    decimalLatitude = paste(co_EXT[3], co_EXT[4], sep = ",")
  )$data

  dat_gbif$Source_type <- "GBIF"
  return(dat_gbif)
}

#' Clean GBIF occurrence records based on quality filters
#'
#' Applies temporal, spatial, uncertainty, elevation, and land/sea filters to flag invalid records.
#'
#' @param flags A data frame of occurrence records to filter.
#' @param year_GBIF Minimum acceptable year (e.g., 1800).
#' @param uncertainty_GBIF Maximum acceptable uncertainty in km.
#' @param GBIF_xmin,GBIF_xmax,GBIF_ymin,GBIF_ymax Coordinate bounds for filtering.
#' @param AltMIN,AltMAX Elevation filtering limits (in meters).
#' @param land Logical. If TRUE, restricts to continental areas.
#'
#' @return A data frame with quality flags and a `Reason` column explaining rejections.
#' @export
sRL_cleanDataGBIF <- function(flags, year_GBIF, uncertainty_GBIF, GBIF_xmin, GBIF_xmax,
                              GBIF_ymin, GBIF_ymax, AltMIN, AltMAX, land) {

  sea_GBIF <- st_read("data/geodata/polygons/Countries/")
  occstemp <- st_as_sf(flags[c("decimalLongitude", "decimalLatitude")],
                       coords = c("decimalLongitude", "decimalLatitude"), crs = crs(sea_GBIF))

  if (land == TRUE) {
    list_sea <- st_within(occstemp, st_make_valid(st_as_sf(sea_GBIF)))
    flags$sea <- ifelse(lengths(list_sea) > 0, "inside", NA)
    flags$.sea <- !is.na(flags$sea)
  } else {
    flags$.sea <- NA
  }

  flags$year <- replace(flags$year, flags$year == 0, NA)
  flags$.year <- flags$year > year_GBIF

  if (!"coordinateUncertaintyInMeters" %in% names(flags)) {
    flags$coordinateUncertaintyInMeters <- NA
  }
  flags$.uncertainty <- as.numeric(flags$coordinateUncertaintyInMeters) < uncertainty_GBIF * 1000

  flags$.limits <- !(flags$decimalLongitude < GBIF_xmin |
                       flags$decimalLongitude > GBIF_xmax |
                       flags$decimalLatitude < GBIF_ymin |
                       flags$decimalLatitude > GBIF_ymax)

  flags$elevationCalculated <- replace(flags$elevationCalculated, flags$elevationCalculated == 0, NA)
  flags$.elevation <- (flags$elevationCalculated < AltMAX & flags$elevationCalculated > AltMIN)

  if ("presence" %in% names(flags)) {
    flags$.pres <- as.character(!flags$presence %in% c("4", "5", "6"))
  }

  flags$Reason <- apply(
    flags[, grepl("^\\.", names(flags))],
    1,
    function(x) {
      reasons <- names(x)[!x & !is.na(x)]
      if (length(reasons) == 0) return(NA)
      plyr::revalue(reasons, c(
        ".val" = "Validity", ".equ" = "Equal_LonLat", ".zer" = "Zero_Coordinates",
        ".cap" = "Capitals", ".cen" = "Country_centroids", ".gbf" = "GBIF_headquarters",
        ".inst" = "Institutions", ".sea" = "Sea", ".year" = "Year",
        ".uncertainty" = "Coordinates_uncertainty", ".limits" = "Outside_extent",
        ".pres" = "Presence_456", ".elevation" = "Elevation"
      ), warn_missing = FALSE) %>% paste(collapse = "; ")
    }
  )

  for (i in seq_len(nrow(flags))) {
    if (!is.na(flags$Link[i])) next
    flags$Link[i] <- paste0("https://www.gbif.org/occurrence/", flags$gbifID[i])
  }

  return(flags)
}

#' Create HTML popup content for map visualization
#'
#' Generates HTML-formatted popup text for each occurrence record to display in interactive maps.
#'
#' @param flags A data frame of cleaned and flagged occurrence records.
#'
#' @return A modified data frame with an added `PopText` column for leaflet popups.
#' @export
sRL_PopRecords <- function(flags){
  flags$Only_for_syn <- ""
  if ("TRUE" %in% grepl("Synonyms_", flags$Source_type)) {
    flags$Only_for_syn <- paste0("<b>Species: </b>", flags$species_download, "<br>")
  }

  flags$PopText <- paste0(
    "<b>", plyr::revalue(as.character(is.na(flags$Reason)),
                         c("TRUE" = "VALID OBSERVATION", "FALSE" = "NOT VALID OBSERVATION")), "</b><br><br>",
    "<b>Source: </b>", flags$Source_type, "<br>",
    flags$Only_for_syn,
    "<b>Observation ID: </b>",
    ifelse(!is.na(flags$Link),
           paste0("<a href='", flags$Link, "' target='_blank'>", flags$gbifID, "</a>"),
           flags$gbifID), "<br>",
    "<b>Year: </b>", flags$year, "<br>",
    "<b>Uncertainty (km): </b>", as.numeric(as.character(flags$coordinateUncertaintyInMeters)) / 1000, "<br>",
    "<b>Elevation: </b>", as.numeric(as.character(flags$elevationCalculated)), "<br>"
  )

  flags$PopText[!is.na(flags$Reason)] <- paste0(flags$PopText[!is.na(flags$Reason)],
                                                "<b>Reason flagged: </b>", flags$Reason[!is.na(flags$Reason)], "<br>")
  return(flags)
}

#' Generate distribution map from occurrence points
#'
#' Builds a species range map from cleaned occurrence records using MCP, kernel, or alpha hull methods.
#'
#' @param dat An `sf` object of occurrence points.
#' @param First_step Character. Mapping method: `"mcp"`, `"kernel"`, `"alpha"`, or `"indivsites"`.
#' @param Buffer_km Numeric. Buffer (in km) to apply to the resulting shape.
#' @param GBIF_crop Character. Cropping method: `"cropsea"` or `"cropland"`.
#' @param Gbif_Param Optional numeric parameter for alpha or kernel smoothing.
#'
#' @return An `sf` object of the resulting buffered and cropped species range.
#' @export
sRL_MapDistributionGBIF <- function(dat, First_step, Buffer_km, GBIF_crop, Gbif_Param) {
  distCountries_mapping <- st_read("data/geodata/sRedList/Red_List_countries_msSimplif0.05_MOLL.shp")
  realms_raw <- st_read("data/geodata/sRedList/RL_Realms.shp")
  CRSMOLL <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"
  realms_mcp <- st_union(realms_raw) %>% st_as_sf() %>% st_transform(CRSMOLL) %>% st_convex_hull() %>% st_buffer(1000)

  if (First_step == "mcp" || First_step == "") {
    distGBIF <- st_as_sf(st_convex_hull(st_union(dat)))
    st_geometry(distGBIF) <- "geometry"
  }

  if (First_step == "kernel") {
    dat_subsample <- distinct(dat, as.character(geometry), .keep_all = TRUE)
    kernel.ref <- adehabitatHR::kernelUD(as_Spatial(dat), h = "href")
    distGBIF <- adehabitatHR::getverticeshr(kernel.ref, percent = 100 * Gbif_Param) %>% st_as_sf()
  }

  if (First_step == "alpha") {
    Round_Fact <- ifelse(as.numeric(max(st_distance(dat))) < 10000, -1, -2)
    Coords_simplif <- st_coordinates(dat) %>% round(digits = Round_Fact)
    dat$Coord_simplif <- paste(Coords_simplif[, 1], Coords_simplif[, 2], sep = ":")
    dat_subsample <- distinct(dat, Coord_simplif, .keep_all = TRUE)

    EX <- raster::extent(dat_subsample)
    Alpha_scaled <- (0.5 * Gbif_Param)^2 * sqrt((EX@xmin - EX@xmax)^2 + (EX@ymin - EX@ymax)^2) %>% as.numeric()

    distGBIF <- spatialEco::convexHull(dat_subsample, alpha = Alpha_scaled)
    st_crs(distGBIF) <- st_crs(dat_subsample)
  }

  if (First_step == "indivsites") {
    distGBIF <- st_buffer(dat, 1)
  }

  distGBIF <- st_buffer(distGBIF, Buffer_km * 1000) %>% st_as_sf()
  return(distGBIF)
}

#' Create interactive map for expert validation
#'
#' Builds an interactive leaflet map combining occurrences, range polygons, and AOO cells.
#'
#' @param flags A data frame of flagged occurrence records (with popup text).
#' @param previous_map An `sf` object of the IUCN 2024 map.
#' @param final_map An `sf` object of the newly proposed map.
#' @param elevation_map Optional `sf` elevation-filtered map.
#' @param occs_filter_spatial `sf` points of valid occurrences.
#'
#' @return A leaflet map widget.
#' @export
sRL_LeafletFlags <- function(flags, previous_map, final_map, elevation_map,
                             occs_filter_spatial){

  grid22 <- raster::raster("data/geodata/sRedList/Empty.grid.2x2.Mollweide.tif")

  occs_proj <- st_transform(occs_filter_spatial, crs(grid22))
  grid_crop <- crop(grid22, (extent(occs_proj) + c(-100000, 100000, -100000, 100000)), snap = "out")

  pts <- occs_proj %>% as_Spatial() %>% as('SpatialPoints')
  AOO_pts <- terra::rasterize(pts, grid_crop, fun = 'count') >= 1
  aooValue <- sum(as.vector(AOO_pts), na.rm = TRUE) * 4

  AOO_bin <- AOO_pts == 1
  AOO_poly <- rasterToPolygons(AOO_bin, fun = function(x){x == 1}, dissolve = FALSE)
  cells <- st_as_sf(AOO_poly)
  cells <- st_transform(cells, crs = st_crs(occs_filter_spatial))

  dist_transf <- st_transform(final_map, crs(grid22))
  grid_crop_map <- crop(grid22, extent(dist_transf), snap = "out")
  grid_crop_map_2 <- grid_crop_map == 1
  grid_crop_poly <- rasterToPolygons(grid_crop_map_2, fun = function(x){x == 1}, dissolve = FALSE)
  cells_total <- st_as_sf(grid_crop_poly)
  cells_total <- st_transform(cells_total, crs = st_crs(final_map))
  cells_cropped <- st_intersection(cells_total, final_map$geometry)

  Leaf <- leaflet(flags) %>%
    addTiles(group = "OpenStreetMap") %>%
    addEsriBasemapLayer(esriBasemapLayers$Imagery, group = "Satellite") %>%
    addEsriBasemapLayer(esriBasemapLayers$Topographic, group = "Topography") %>%
    addPolygons(data = as_Spatial(cells_cropped), fillColor = "#d3d3d3", fillOpacity = 0.01,
                color = "#d3d3d3") %>%
    addPolygons(data = as_Spatial(cells), fillColor = "#38b000", fillOpacity = 0.1,
                color = "#38b000") %>%
    addPolygons(data = as_Spatial(previous_map), fillColor = "#E84A5F", fillOpacity = 0.1,
                color = "#E84A5F") %>%
    addPolygons(data = as_Spatial(final_map), fillColor = "transparent",
                color = "black") %>%
    addCircleMarkers(
      lng = flags$decimalLongitude,
      lat = flags$decimalLatitude,
      color = ifelse(is.na(flags$Reason) == TRUE, "#fdcb25ff", "#440154ff"),
      fillOpacity = 0.8,
      stroke = FALSE,
      popup = flags$PopText,
      radius = 8,
      group = "Occurrence records"
    ) %>%
    addLegend(position = "bottomleft", colors = c('#fdcb25ff', '#440154ff'),
              labels = c("Valid", "Not valid")) %>%
    addLayersControl(baseGroups = c("OpenStreetMap", "Satellite", "Topography"),
                     overlayGroups = "Occurrence records",
                     position = "topleft") %>%
    addMouseCoordinates() %>%
    addScaleBar(position = "bottomright") %>%
    addLegend(
      colors = if (!is.null(elevation_map)) c("#F7AAA6", "#686868", "#34a0a4", "#38b000")
      else c("#F7AAA6", "#686868", "#38b000"),
      labels = if (!is.null(elevation_map)) c("Previous range", "New range", "Elevation limit", "Occupied cells")
      else c("Previous range", "New range", "Occupied cells"),
      title = "Ranges",
      opacity = 1,
      position = "bottomleft"
    )

  saveWidget(Leaf, file = paste0("species/", spname, "/", spname, "_map_experts.html"))
  return(Leaf)
}

