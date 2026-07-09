# =============================================================================
# 03_compute_THI_seasons.R
#
# Calcula, a partir de datos observacionales (ERA5), el nº de días bajo
# estrés termico (THI) por estacion del anio, especie ganadera y categoria
# de severidad (mild/moderate/severe). Guarda el resultado en un RDS que
# consume plot_THI_seasons.R
# =============================================================================

# ---- Librerias --------------------------------------------------------------
library(convertR)
library(transformeR)
library(loadeR)
library(sp)

# ---- Funciones propias -------------------------------------------------------
source(file.path(dir.project <- "/lustre/gmeteo/WORK/uribei/TFM",
                  "Scripts/THI/functions_THI.R"))

computeTrend <- function(ts) {  # compute trends
  df <- data.frame(x = 1:length(ts), y = ts)
  if (sum(is.na(ts)) < round(0.75 * length(ts))) {  # ask for a minimum of 75% of non-missing data to compute the trend
    return(lm(y ~ x, df)$coefficients[2])
  } else {
    return(NA)
  }
}

computeSigTrend <- function(ts) {  # compute trends
  df <- data.frame(x = 1:length(ts), y = ts)
  if (sum(is.na(ts)) < round(0.75 * length(ts))) {  # ask for a minimum of 75% of non-missing data to compute the trend
    return(summary(lm(y ~ x, df))$coefficients[, 4][2])
  } else {
    return(NA)
  }
}

# ---- Parametros ---------------------------------------------------------------
dir.data <- file.path(dir.project, "Data")
dir.obs  <- file.path(dir.project, "Data/obs")

f.tas  <- file.path(dir.data, "tas_obs.rda")
f.hurs <- file.path(dir.data, "hurs_obs.rda")

per_hist <- 1986:2005

seasons <- list(dfj = c(12, 1, 2), mam = 3:5, jja = 6:8, son = 9:11)

species     <- c("cattle", "cattle", "ruminants", "ruminants", "poultry", "swine")
species.lab <- c("dairy", "beef", "sheep", "goats", "poultry", "swine")

mild     <- list(c(72, 78),  c(75, 81),  c(25, 30), c(27, 32), c(27.8, 28.8),  c(23.33, 26.11))
moderate <- list(c(78, 88),  c(81, 91),  c(30, 35), c(32, 37), c(28.8, 30),    c(26.11, 28.88))
severe   <- list(c(88, Inf), c(91, Inf), c(35, Inf), c(37, Inf), c(30, Inf),   c(28.88, Inf))

thresholds.df <- data.frame(rbind(mild, moderate, severe))
colnames(thresholds.df) <- species.lab
rownames(thresholds.df) <- c("mild", "moderate", "severe")

categories <- rownames(thresholds.df)

f.out <- file.path(dir.obs, sprintf("THI_obs_ndays_by_season_%d-%d.rds",
                                     min(per_hist), max(per_hist)))

# ---- Comprobaciones -----------------------------------------------------------
stopifnot(
  file.exists(f.tas),
  file.exists(f.hurs),
  dir.exists(dir.obs)
)

# ---- Carga de datos observacionales --------------------------------------------
load(f.tas)   # tas.era5
load(f.hurs)  # hurs.era5

tas.era5$Variable$varName <- "tas"
grid_obs <- makeMultiGrid("tas" = tas.era5, "hurs" = hurs.era5)
grid_obs <- subsetGrid(grid_obs, years = per_hist)

rm(tas.era5, hurs.era5)
gc()

# ---- Calculo de THI y nº de dias por categoria/especie/estacion ----------------
maps.by.cat <- list()
for (sp in species.lab) {
  maps.by.cat[[sp]] <- list()
  for (category in categories) {
    maps.by.cat[[sp]][[category]] <- list()
  }
}

for (i in seq_along(species)) {

  sp_lab <- species.lab[i]
  sp_fun <- species[i]

  for (s in names(seasons)) {

    thi <- list()

    for (M in seasons[[s]]) {

      # tas
      data <- subsetGrid(grid_obs, season = M, years = per_hist, var = "tas")
      attributes(data$Variable)$units <- "degC"
      tas <- data
      rm(data)

      # hurs
      data <- subsetGrid(grid_obs, season = M, years = per_hist, var = "hurs")
      hurs <- data
      rm(data)

      thi[[as.character(M)]] <- computeTHI(tas, hurs, species = sp_fun)
    }

    thi.D <- do.call("bindGrid", c(thi, list(dimension = "time")))

    for (j in seq_along(categories)) {

      category  <- categories[j]
      threshold <- thresholds.df[[sp_lab]][[j]]

      thi.Y <- analyzeTHI_Y(
        thi           = thi.D,
        threshold.min = threshold[1],
        threshold.max = threshold[2]
      )

      maps.by.cat[[sp_lab]][[category]][[s]] <- thi.Y$ndays
    }

    rm(thi, thi.D)
    gc()
  }
}

# ---- Guardado -------------------------------------------------------------------
obs_THI <- list(
  maps.by.cat = maps.by.cat,
  species.lab = species.lab,
  categories  = categories,
  seasons     = seasons,
  per_hist    = per_hist
)

saveRDS(obs_THI, file = f.out)

message("Guardado: ", f.out)