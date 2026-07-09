# =============================================================================
# 07_estadisticas.R
# Calcula valores medios observados y valores futuros por especie, categoria
# y regiones ganaderas relevantes.
#
# IMPORTANTE:
# Los valores futuros se calculan como:
#
#   futuro_estimado_GCM = observado_ERA5_historico + delta_GCM
#
# donde:
#
#   delta_GCM = futuro_corregido_GCM - historico_corregido_GCM
#
# La incertidumbre se representa con la horquilla:
#
#   [min(futuro_MPI, futuro_HadGEM, futuro_NorESM),
#    max(futuro_MPI, futuro_HadGEM, futuro_NorESM)]
# =============================================================================

# ---- Librerias ----------------------------------------------------------------
library(sf)
library(transformeR)

# ---- Funciones auxiliares -------------------------------------------------------
# (ordenadas de mas basica a mas compuesta: cada una usa las anteriores)

select_ccaa <- function(ccaa_name) {

  out <- ccaa[ccaa[[name_col]] == ccaa_name, ]

  if (nrow(out) == 0) {
    stop(
      "No se ha encontrado la CCAA: ",
      ccaa_name,
      ". Nombres disponibles: ",
      paste(unique(ccaa[[name_col]]), collapse = ", ")
    )
  }

  out
}

get_lon_lat <- function(grid) {

  lon <- grid$xyCoords$x
  lat <- grid$xyCoords$y

  if (is.null(lon)) lon <- grid$xyCoords$lon
  if (is.null(lat)) lat <- grid$xyCoords$lat

  list(
    lon = as.numeric(lon),
    lat = as.numeric(lat)
  )
}

build_mask <- function(grid, ccaa_names) {

  coords <- get_lon_lat(grid)
  lon <- coords$lon
  lat <- coords$lat

  xy <- expand.grid(lon = lon, lat = lat)

  pts <- st_as_sf(xy, coords = c("lon", "lat"), crs = 4326, remove = FALSE)

  selected <- do.call(rbind, lapply(ccaa_names, select_ccaa))
  selected_union <- st_union(selected)

  inside <- st_intersects(pts, selected_union, sparse = FALSE)[, 1]

  matrix(
    inside,
    nrow = length(lat),
    ncol = length(lon),
    byrow = TRUE
  )
}

extract_grid_values <- function(grid, mask) {

  grid_clim <- climatology(grid)
  dat <- grid_clim$Data

  dims <- dim(dat)
  idx <- which(mask, arr.ind = TRUE)

  if (length(dims) == 2) {

    vals <- dat[mask]

  } else if (length(dims) == 3) {

    vals <- unlist(
      lapply(seq_len(nrow(idx)), function(k) {
        dat[, idx[k, 1], idx[k, 2]]
      })
    )

  } else if (length(dims) == 4) {

    vals <- unlist(
      lapply(seq_len(nrow(idx)), function(k) {
        dat[, , idx[k, 1], idx[k, 2]]
      })
    )

  } else {

    stop(
      "Estructura de grid$Data no esperada despues de climatology(). Dimensiones: ",
      paste(dims, collapse = " x ")
    )
  }

  vals <- as.numeric(vals)
  vals[!is.na(vals)]
}

safe_mean <- function(x) {
  if (length(x) == 0) return(NA_real_)
  mean(x, na.rm = TRUE)
}

safe_min <- function(x) {
  if (length(x) == 0) return(NA_real_)
  min(x, na.rm = TRUE)
}

safe_max <- function(x) {
  if (length(x) == 0) return(NA_real_)
  max(x, na.rm = TRUE)
}

get_delta_grid <- function(delta_obj, sp, method, category) {

  if (is.null(delta_obj$list.maps.ndays.Y[[sp]][[method]][[category]])) {
    stop(
      "No se ha encontrado la delta para especie = ", sp,
      ", metodo = ", method,
      ", categoria = ", category
    )
  }

  delta_obj$list.maps.ndays.Y[[sp]][[method]][[category]]
}

summarise_projected_from_deltas <- function(delta_grids_by_model, obs_grid, mask) {

  obs_vals <- extract_grid_values(obs_grid, mask)
  obs_mean <- safe_mean(obs_vals)

  delta_means <- sapply(names(delta_grids_by_model), function(model) {
    vals <- extract_grid_values(delta_grids_by_model[[model]], mask)
    safe_mean(vals)
  })

  # Valor futuro estimado = observado historico + delta proyectada
  projected_values <- obs_mean + delta_means

  data.frame(
    obs_mean = obs_mean,
    mean_days = mean(projected_values, na.rm = TRUE),
    min_days = min(projected_values, na.rm = TRUE),
    max_days = max(projected_values, na.rm = TRUE),
    MPI = projected_values["MPI"],
    HadGEM = projected_values["HadGEM"],
    NorESM = projected_values["NorESM"],
    delta_MPI = delta_means["MPI"],
    delta_HadGEM = delta_means["HadGEM"],
    delta_NorESM = delta_means["NorESM"]
  )
}

# ---- Configuracion ---------------------------------------------------------------

season_to_analyse <- "jja"        # "dfj", "mam", "jja", "son"
gwl_to_analyse <- "GWL1.5"        # "GWL1.5", "GWL2", "GWL3"

models <- c("MPI", "HadGEM", "NorESM")

methods_future <- c("mbcn")

species_to_use <- c("dairy", "beef", "sheep", "goats", "swine")
categories <- c("mild", "moderate", "severe")

species_regions <- list(
  dairy = c("Galicia", "Principado de Asturias", "Cantabria"),
  swine = c("Extremadura", "Andalucía"),
  beef = c("Extremadura", "Castilla y León"),
  sheep = c("Extremadura", "Castilla y León", "Castilla-La Mancha"),
  goats = c("Extremadura", "Andalucía")
)

dir_project <- "~/lustre/gmeteo/WORK/uribei/TFM"

dir_obs <- file.path(dir_project, "Data/obs")
obs_file <- file.path(dir_obs, "THI_obs_ndays_by_season_1986-2005.rds")

delta_file_model <- function(model) {
  file.path(
    dir_project, "Data",
    model,
    gwl_to_analyse,
    paste0("THI_deltas_results_", gwl_to_analyse, "_", toupper(season_to_analyse), ".rds")
  )
}

ccaa_shp <- file.path(dir_project, "Scripts/España.shp/gadm41_ESP_1.shp")

out_dir <- file.path(dir_project, "Data/REGIONAL_VALUES")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

name_col <- "NAME_1"

# ---- Cargar datos -----------------------------------------------------------------

if (!file.exists(obs_file)) {
  stop("No existe el archivo de observaciones: ", obs_file)
}

obs_THI <- readRDS(obs_file)

delta_THI_models <- list()

for (model in models) {

  f <- delta_file_model(model)

  if (!file.exists(f)) {
    stop("No existe el archivo de deltas para ", model, ": ", f)
  }

  delta_THI_models[[model]] <- readRDS(f)
}

ccaa <- st_read(ccaa_shp, quiet = TRUE)
ccaa <- st_transform(ccaa, 4326)
ccaa <- st_make_valid(ccaa)

# -----------------------------------------------------------------------------
# OBSERVACIONES POR REGION COMBINADA DE CADA ESPECIE
# -----------------------------------------------------------------------------

obs_combined_rows <- list()
k <- 1

for (sp in species_to_use) {

  ccaa_names <- species_regions[[sp]]

  for (category in categories) {

    grid_obs <- obs_THI$maps.by.cat[[sp]][[category]][[season_to_analyse]]

    mask <- build_mask(grid_obs, ccaa_names)

    vals <- extract_grid_values(grid_obs, mask)
    mean_vals <- safe_mean(vals)

    obs_combined_rows[[k]] <- data.frame(
      period = "historical",
      source = "ERA5",
      gwl = NA,
      season = toupper(season_to_analyse),
      species = sp,
      category = category,
      region = paste(ccaa_names, collapse = " + "),
      obs_mean = mean_vals,
      mean_days = mean_vals,
      min_days = safe_min(vals),
      max_days = safe_max(vals),
      n_values = length(vals)
    )

    k <- k + 1
  }
}

obs_summary_species_region <- do.call(rbind, obs_combined_rows)

# -----------------------------------------------------------------------------
# FUTURO POR REGION COMBINADA DE CADA ESPECIE
# Se calcula:
#
#   futuro_GCM = observado_ERA5 + delta_GCM
#
# Luego:
#
#   mean_days = media de futuro_MPI, futuro_HadGEM y futuro_NorESM
#   min_days  = minimo entre futuro_MPI, futuro_HadGEM y futuro_NorESM
#   max_days  = maximo entre futuro_MPI, futuro_HadGEM y futuro_NorESM
# -----------------------------------------------------------------------------

future_combined_rows <- list()
k <- 1

for (sp in species_to_use) {

  ccaa_names <- species_regions[[sp]]

  for (method in methods_future) {

    for (category in categories) {

      delta_grids_by_model <- lapply(models, function(model) {
        get_delta_grid(
          delta_obj = delta_THI_models[[model]],
          sp = sp,
          method = method,
          category = category
        )
      })

      names(delta_grids_by_model) <- models

      grid_obs <- obs_THI$maps.by.cat[[sp]][[category]][[season_to_analyse]]

      mask <- build_mask(grid_obs, ccaa_names)

      gcm_summary <- summarise_projected_from_deltas(
        delta_grids_by_model = delta_grids_by_model,
        obs_grid = grid_obs,
        mask = mask
      )

      future_combined_rows[[k]] <- data.frame(
        period = "future",
        source = method,
        gwl = gwl_to_analyse,
        season = toupper(season_to_analyse),
        species = sp,
        category = category,
        region = paste(ccaa_names, collapse = " + "),
        obs_mean = gcm_summary$obs_mean,
        mean_days = gcm_summary$mean_days,
        min_days = gcm_summary$min_days,
        max_days = gcm_summary$max_days,
        n_values = NA_integer_,
        MPI = gcm_summary$MPI,
        HadGEM = gcm_summary$HadGEM,
        NorESM = gcm_summary$NorESM,
        delta_MPI = gcm_summary$delta_MPI,
        delta_HadGEM = gcm_summary$delta_HadGEM,
        delta_NorESM = gcm_summary$delta_NorESM
      )

      k <- k + 1
    }
  }
}

future_summary_species_region <- do.call(rbind, future_combined_rows)

# -----------------------------------------------------------------------------
# UNIFICAR COLUMNAS
# -----------------------------------------------------------------------------

obs_summary_species_region$MPI <- NA_real_
obs_summary_species_region$HadGEM <- NA_real_
obs_summary_species_region$NorESM <- NA_real_
obs_summary_species_region$delta_MPI <- NA_real_
obs_summary_species_region$delta_HadGEM <- NA_real_
obs_summary_species_region$delta_NorESM <- NA_real_

common_cols <- c(
  "period",
  "source",
  "gwl",
  "season",
  "species",
  "category",
  "region",
  "obs_mean",
  "mean_days",
  "min_days",
  "max_days",
  "n_values",
  "MPI",
  "HadGEM",
  "NorESM",
  "delta_MPI",
  "delta_HadGEM",
  "delta_NorESM"
)

obs_summary_species_region <- obs_summary_species_region[, common_cols]
future_summary_species_region <- future_summary_species_region[, common_cols]

summary_species_region <- rbind(
  obs_summary_species_region,
  future_summary_species_region
)

# -----------------------------------------------------------------------------
# GUARDAR RESULTADOS
# -----------------------------------------------------------------------------

out_file <- file.path(
  out_dir,
  paste0(
    "THI_obs_future_summary_region_",
    gwl_to_analyse,
    "_",
    toupper(season_to_analyse),
    "_OBS_plus_delta.csv"
  )
)

write.csv2(
  summary_species_region,
  file = out_file,
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# IMPRIMIR RESUMEN
# -----------------------------------------------------------------------------

print("Resumen observaciones por region combinada:")
print(obs_summary_species_region)

print("Resumen futuro por region combinada:")
print(future_summary_species_region)

print("Resumen total guardado:")
print(summary_species_region)

print(paste("Archivo guardado en:", out_file))