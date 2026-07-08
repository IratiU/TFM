################################################################################
# 00_load_data.R
#
# Carga los datos observados (ERA5), la simulación histórica del RCM y la
# proyección futura RCP8.5 para las simulaciones REMO usadas en el TFM
# (NorESM, HadGEM, MPI). Qué partes cargar se controla con los parámetros
# cargar_obs / cargar_rcm / cargar_rcp.
#
# Salidas:
#   Data/tas_obs.rda, Data/hurs_obs.rda, Data/grid_obs.rda
#   Data/<model>/tas_rcm.rda, Data/<model>/hurs_rcm.rda, Data/<model>/grid_rcm.rda
#   Data/<model>/tas_rcp.rda, Data/<model>/hurs_rcp.rda, Data/<model>/grid_rcp.rda
################################################################################


################################################################################
# 1. PARÁMETROS
################################################################################

# Qué partes cargar
cargar_obs <- TRUE
cargar_rcm <- TRUE
cargar_rcp <- TRUE

# Modelo a procesar (solo hace falta si cargar_rcm y/o cargar_rcp = TRUE)
model <- "NorESM"  # "NorESM", "HadGEM" o "MPI"

# Dominio espacial
lonLim <- c(-10, 5)
latLim <- c(35, 45)

# Periodo histórico (observado y RCM)
years.hist <- 1986:2005

# Periodo a conservar de la proyección RCP8.5 (NULL = todo el periodo disponible)
years.rcp <- NULL

# Directorio principal del proyecto
dir.project <- "/lustre/gmeteo/WORK/uribei/TFM"

# Shapefile para recortar la malla al territorio de interés
shp.path <- file.path(dir.project, "Scripts", "España.shp", "gadm41_ESP_0.shp")


################################################################################
# 2. PAQUETES
################################################################################

library(rJava)
.jaddClassPath("/nfs/home/gmeteo/uribei/miniconda3/envs/c4r/lib/R/library/loadeR.java/java/netcdfAll-4.6.0-SNAPSHOT.jar")
.jinit(classpath = "/nfs/home/gmeteo/uribei/miniconda3/envs/c4r/lib/R/library/loadeR.java/java/netcdfAll-4.6.0-SNAPSHOT.jar")

library(loadeR)
library(sf)
library(transformeR)
library(visualizeR)
library(convertR)


################################################################################
# 3. DIRECTORIOS
################################################################################

dir  <- file.path(dir.project, "Data")
diro <- file.path(dir, model)

dir.create(diro, recursive = TRUE, showWarnings = FALSE)


################################################################################
# 4. SHAPEFILE Y FUNCIONES AUXILIARES
################################################################################

shp <- st_read(shp.path)
shp <- st_make_valid(shp)
shp <- st_transform(shp, crs = 4326)

# Enmascara con NA las celdas de la malla que caen fuera del shapefile.
trimGrid <- function(grid, shp) {

  lon <- grid$xyCoords$x
  lat <- grid$xyCoords$y

  pts <- expand.grid(lon = lon, lat = lat)
  pts_sf <- st_as_sf(pts, coords = c("lon", "lat"), crs = 4326)

  ins <- st_intersects(pts_sf, shp, sparse = FALSE)
  ins <- apply(ins, 1, any)

  mask <- matrix(ins, nrow = length(lat), ncol = length(lon), byrow = TRUE)

  for (t in seq_len(dim(grid$Data)[1])) {
    aux <- grid$Data[t, , ]
    aux[!mask] <- NA
    grid$Data[t, , ] <- aux
  }

  grid
}

# Descarga y ensambla una variable CORDEX a partir de un catálogo THREDDS:
# lista los .nc disponibles, los carga uno a uno y los une en el tiempo.
load_cordex_var <- function(base_catalog, base_dods, var, version,
                            lonLim, latLim, shp, label = "") {

  catalog_url <- paste0(base_catalog, var, "/", version, "/catalog.html")
  cat_lines <- readLines(catalog_url, warn = FALSE)

  nc_files <- sub(
    ".*<code>([^<]+\\.nc)</code>.*",
    "\\1",
    grep("\\.nc", cat_lines, value = TRUE)
  )
  nc_files <- unique(nc_files)
  nc_files <- nc_files[grepl("\\.nc$", nc_files)]

  datasets <- paste0(base_dods, var, "/", version, "/", nc_files)

  grid_list <- vector("list", length(datasets))

  for (i in seq_along(datasets)) {
    message("--------------------------------------------------")
    message("Cargando ", toupper(var), " ", label, ": ", nc_files[i])
    message("URL: ", datasets[i])

    tmp <- loadGridData(datasets[i], var = var, lonLim = lonLim, latLim = latLim)
    tmp <- trimGrid(tmp, shp)

    grid_list[[i]] <- tmp

    rm(tmp)
    gc()
  }

  do.call("bindGrid", c(grid_list, list(dimension = "time")))
}


################################################################################
# 5. DATOS OBSERVADOS (ERA5)
################################################################################

if (cargar_obs) {

  dataset.era5 <- "https://hub.climate4r.ifca.es/thredds/dodsC/fao/observations/ERA5/0.25/ERA5_025.ncml"

  hurs.era5 <- loadGridData(
    dataset.era5,
    var = "hurs", years = years.hist,
    lonLim = lonLim, latLim = latLim
  )
  hurs.era5 <- trimGrid(hurs.era5, shp)
  save(hurs.era5, file = file.path(dir, "hurs_obs.rda"))

  tas.era5 <- loadGridData(
    dataset.era5,
    var = "tas", years = years.hist,
    lonLim = lonLim, latLim = latLim
  )
  tas.era5 <- trimGrid(tas.era5, shp)
  tas.era5 <- udConvertGrid(tas.era5, new.units = "celsius")
  tas.era5$Variable$varName <- "tas"
  save(tas.era5, file = file.path(dir, "tas_obs.rda"))

  grid.obs <- makeMultiGrid(tas = tas.era5, hurs = hurs.era5)
  save(grid.obs, file = file.path(dir, "grid_obs.rda"))

  rm(hurs.era5, tas.era5, grid.obs)
  gc()
}


################################################################################
# 6. AJUSTES DEL CATÁLOGO DEL MODELO (necesario para RCM y/o RCP)
################################################################################

if (cargar_rcm || cargar_rcp) {

  gcm <- switch(
    model,
    "NorESM" = "NCC-NorESM1-M",
    "HadGEM" = "MOHC-HadGEM2-ES",
    "MPI"    = "MPI-M-MPI-ESM-LR"
  )

  version.hist <- switch(
    model,
    "NorESM" = "v20191029",
    "HadGEM" = "v20191029",
    "MPI"    = "v20191015"
  )

  version.rcp <- switch(
    model,
    "NorESM" = "v20191029",
    "HadGEM" = "v20191029",
    "MPI"    = "v20191029"
  )
}


################################################################################
# 7. SIMULACIÓN HISTÓRICA DEL RCM
################################################################################

if (cargar_rcm) {

  base.catalog.hist <- paste0(
    "https://hub.climate4r.ifca.es/thredds/catalog/fao/interp025/raw/",
    "CORDEX/output/EUR-22/GERICS/", gcm, "/historical/r1i1p1/REMO2015/v1/day/"
  )
  base.dods.hist <- paste0(
    "https://hub.climate4r.ifca.es/thredds/dodsC/fao/interp025/raw/",
    "CORDEX/output/EUR-22/GERICS/", gcm, "/historical/r1i1p1/REMO2015/v1/day/"
  )

  hurs.rcm <- load_cordex_var(
    base.catalog.hist, base.dods.hist, "hurs", version.hist,
    lonLim, latLim, shp, label = "historical"
  )
  hurs.rcm <- subsetGrid(hurs.rcm, years = years.hist)
  save(hurs.rcm, file = file.path(diro, "hurs_rcm.rda"))

  tas.rcm <- load_cordex_var(
    base.catalog.hist, base.dods.hist, "tas", version.hist,
    lonLim, latLim, shp, label = "historical"
  )
  tas.rcm <- subsetGrid(tas.rcm, years = years.hist)
  tas.rcm <- udConvertGrid(tas.rcm, new.units = "celsius")
  save(tas.rcm, file = file.path(diro, "tas_rcm.rda"))

  grid.rcm <- makeMultiGrid(tas = tas.rcm, hurs = hurs.rcm)
  save(grid.rcm, file = file.path(diro, "grid_rcm.rda"))

  rm(hurs.rcm, tas.rcm, grid.rcm)
  gc()
}


################################################################################
# 8. PROYECCIÓN FUTURA RCP8.5
################################################################################

if (cargar_rcp) {

  base.catalog.rcp <- paste0(
    "https://hub.climate4r.ifca.es/thredds/catalog/fao/interp025/raw/",
    "CORDEX/output/EUR-22/GERICS/", gcm, "/rcp85/r1i1p1/REMO2015/v1/day/"
  )
  base.dods.rcp <- paste0(
    "https://hub.climate4r.ifca.es/thredds/dodsC/fao/interp025/raw/",
    "CORDEX/output/EUR-22/GERICS/", gcm, "/rcp85/r1i1p1/REMO2015/v1/day/"
  )

  hurs.rcp <- load_cordex_var(
    base.catalog.rcp, base.dods.rcp, "hurs", version.rcp,
    lonLim, latLim, shp, label = "RCP8.5"
  )
  if (!is.null(years.rcp)) hurs.rcp <- subsetGrid(hurs.rcp, years = years.rcp)
  save(hurs.rcp, file = file.path(diro, "hurs_rcp.rda"))

  tas.rcp <- load_cordex_var(
    base.catalog.rcp, base.dods.rcp, "tas", version.rcp,
    lonLim, latLim, shp, label = "RCP8.5"
  )
  if (!is.null(years.rcp)) tas.rcp <- subsetGrid(tas.rcp, years = years.rcp)
  tas.rcp <- udConvertGrid(tas.rcp, new.units = "celsius")
  save(tas.rcp, file = file.path(diro, "tas_rcp.rda"))

  grid.rcp <- makeMultiGrid(tas = tas.rcp, hurs = hurs.rcp)
  save(grid.rcp, file = file.path(diro, "grid_rcp.rda"))

  rm(hurs.rcp, tas.rcp, grid.rcp)
  gc()
}


################################################################################
# 9. COMPROBACIÓN DE SALIDAS
################################################################################

list.files(dir, pattern = "_obs\\.rda$", full.names = TRUE)
list.files(diro, pattern = "\\.rda$", full.names = TRUE)
