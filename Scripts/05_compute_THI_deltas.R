# =============================================================================
# 05_compute_THI_deltas.R
# Calcula el THI (hist y pred) y las deltas de dias de estres por especie,
# metodo y categoria. Itera sobre todos los modelos y guarda un RDS por modelo.
# =============================================================================

library(gridExtra)
library(grid)
library(visualizeR)
library(RColorBrewer)
library(lattice)
library(convertR)
library(transformeR)
library(loadeR)

dir.project <- "/lustre/gmeteo/WORK/uribei/TFM"
source(file.path(dir.project, "Scripts", "THI", "functions_THI.R"))

# Nota: definidas aqui pero no se usan en este script; si ya no hacen falta,
# se pueden quitar (o mover a functions_THI.R si se usan en otro sitio).
computeTrend <- function(ts) {
  df <- data.frame(x = 1:length(ts), y = ts)
  if (sum(is.na(ts)) < round(0.75 * length(ts))) return(lm(y ~ x, df)$coefficients[2])
  return(NA)
}

computeSigTrend <- function(ts) {
  df <- data.frame(x = 1:length(ts), y = ts)
  if (sum(is.na(ts)) < round(0.75 * length(ts))) return(summary(lm(y ~ x, df))$coefficients[, 4][2])
  return(NA)
}


# -----------------------------------------------------------------------------
# CONFIGURACION GLOBAL
# -----------------------------------------------------------------------------

models <- c("MPI", "HadGEM", "NorESM")
met    <- c("eqm", "qdm", "mbcn")

seasons <- list(ANN = 1:12, DJF = c(1, 2, 12), MAM = 3:5, JJA = 6:8, SON = 9:11)
season_to_compute <- "MAM"
stopifnot(season_to_compute %in% names(seasons))
months_to_compute <- seasons[[season_to_compute]]

# GWL a calcular
gwl_to_compute <- "GWL3"

gwl_periods <- list(
  HadGEM = list("GWL1.5" = 2014:2033, "GWL2" = 2026:2045, "GWL3" = 2045:2064),
  MPI    = list("GWL1.5" = 2008:2027, "GWL2" = 2028:2047, "GWL3" = 2052:2071),
  NorESM = list("GWL1.5" = 2023:2042, "GWL2" = 2039:2058, "GWL3" = 2063:2082)
)

stopifnot(all(models %in% names(gwl_periods)))
stopifnot(gwl_to_compute %in% c("GWL1.5", "GWL2", "GWL3"))

per_hist <- 1986:2005

# Rutas: solo cambia el nombre del modelo
dir_base <- function(model) file.path(dir.project, "Data", model)
dir_gwl  <- function(model) file.path(dir_base(model), gwl_to_compute)
dir_hist <- function(model) file.path(dir_base(model), "hist")

# Especies y umbrales de estres
species.lab <- c("dairy", "beef", "sheep", "goats", "poultry", "swine")
mild        <- list(c(72, 78),  c(75, 81),  c(25, 30),  c(27, 32),  c(27.8, 28.8), c(23.33, 26.11))
moderate    <- list(c(78, 88),  c(81, 91),  c(30, 35),  c(32, 37),  c(28.8, 30),   c(26.11, 28.88))
severe      <- list(c(88, Inf), c(91, Inf), c(35, Inf), c(37, Inf), c(30, Inf),    c(28.88, Inf))

df <- data.frame(rbind(mild, moderate, severe))
colnames(df) <- species.lab
rownames(df) <- c("mild", "moderate", "severe")


# -----------------------------------------------------------------------------
# FUNCION AUXILIAR: computa THI mensual y lo une en un solo grid temporal
# -----------------------------------------------------------------------------
compute_thi_list <- function(grid, years, months = 1:12) {

  thi.cattle    <- vector("list", length(months))
  thi.ruminants <- vector("list", length(months))
  thi.swine     <- vector("list", length(months))
  thi.poultry   <- vector("list", length(months))

  for (k in seq_along(months)) {

    M <- months[k]

    tas.m  <- subsetGrid(grid, season = M, years = years, var = "tas")
    attributes(tas.m$Variable)$units <- "degC"

    hurs.m <- subsetGrid(grid, season = M, years = years, var = "hurs")

    thi.cattle[[k]]    <- computeTHI(tas.m, hurs.m, species = "cattle")
    thi.ruminants[[k]] <- computeTHI(tas.m, hurs.m, species = "ruminants")
    thi.swine[[k]]     <- computeTHI(tas.m, hurs.m, species = "swine")
    thi.poultry[[k]]   <- computeTHI(tas.m, hurs.m, species = "poultry")
  }

  # dairy/beef comparten la formula "cattle" y sheep/goats la formula
  # "ruminants": se une una vez por grupo y se reutiliza en ambas especies.
  thi.cattle.grid    <- bindGrid(thi.cattle,    dimension = "time")
  thi.ruminants.grid <- bindGrid(thi.ruminants, dimension = "time")

  list(
    dairy   = thi.cattle.grid,
    beef    = thi.cattle.grid,
    sheep   = thi.ruminants.grid,
    goats   = thi.ruminants.grid,
    poultry = bindGrid(thi.poultry, dimension = "time"),
    swine   = bindGrid(thi.swine,   dimension = "time")
  )
}

# -----------------------------------------------------------------------------
# FUNCION AUXILIAR: diferencia de dias de estres entre metodos (mbcn - X)
# -----------------------------------------------------------------------------
compute_ndays_diff <- function(ndays_list, comparator, cats) {

  diffs  <- setNames(vector("list", length(cats)), cats)
  titles <- setNames(vector("list", length(cats)), cats)

  for (cat in cats) {
    diffs[[cat]] <- gridArithmetics(
      ndays_list[["mbcn"]][[cat]],
      ndays_list[[comparator]][[cat]],
      operator = "-"
    )
    titles[[cat]] <- sprintf("MBCn - %s - %s", toupper(comparator), toupper(cat))
  }

  list(diffs = diffs, titles = titles)
}


# -----------------------------------------------------------------------------
# BUCLE PRINCIPAL POR MODELO
# -----------------------------------------------------------------------------
for (model in models) {
  per <- gwl_periods[[model]][[gwl_to_compute]]
  message("\n============================")
  message(sprintf("  MODELO: %s", model))
  message("============================")

  # --- Carga de datos --------------------------------------------------------
  message(sprintf("  Cargando datos futuros RCP (%s)...", model))
  load(file.path(dir_base(model), "tas_rcp.rda"))
  load(file.path(dir_base(model), "hurs_rcp.rda"))
  grid_rcp <- makeMultiGrid("tas" = tas.rcp, "hurs" = hurs.rcp)
  grid_rcp <- subsetGrid(grid_rcp, years = per)

  message(sprintf("  Cargando datos historicos RCM (%s)...", model))
  load(file.path(dir_base(model), "tas_rcm.rda"))
  load(file.path(dir_base(model), "hurs_rcm.rda"))
  grid_rcm <- makeMultiGrid("tas" = tas.rcm, "hurs" = hurs.rcm)
  grid_rcm <- subsetGrid(grid_rcm, years = per_hist)

  # Listas de grids por metodo
  methods.pred <- list(rcm = grid_rcp)
  methods.hist <- list(rcm = grid_rcm)

  for (m in met) {
    methods.pred[[m]] <- readRDS(file.path(dir_gwl(model),  paste0(m, "_pred.rds")))
    methods.hist[[m]] <- readRDS(file.path(dir_hist(model), paste0(m, "_hist.rds")))
  }


  # --- Calculo ---------------------------------------------------------------
  maps.dc                  <- list()
  list.maps.ndays.Y        <- list()
  tit.days                 <- list()
  list.maps.ndays.diff     <- list()
  tit.days.diff            <- list()
  list.maps.ndays.diff.eqm <- list()
  tit.days.diff.eqm        <- list()


  for (m in names(methods.pred)) {
    message(sprintf("  [%s] Metodo: %s", model, toupper(m)))

    thi.pred <- compute_thi_list(methods.pred[[m]], per, months = months_to_compute)
    thi.hist <- compute_thi_list(methods.hist[[m]], per_hist, months = months_to_compute)

    # Delta climatologico del THI
    maps.dc[[m]] <- lapply(names(thi.pred), function(s) {
      gridArithmetics(climatology(thi.pred[[s]]),
                      climatology(thi.hist[[s]]),
                      operator = "-")
    })
    names(maps.dc[[m]]) <- names(thi.pred)

    # Delta de dias de estres por especie y categoria
    for (i in seq_along(species.lab)) {
      s <- species.lab[i]
      message(sprintf("    Especie: %s", toupper(s)))

      for (j in seq_len(nrow(df))) {
        cat <- rownames(df)[j]

        thi.pred.Y <- analyzeTHI_Y(thi           = thi.pred[[s]],
                                   threshold.min = df[[j, i]][1],
                                   threshold.max = df[[j, i]][2])

        thi.hist.Y <- analyzeTHI_Y(thi           = thi.hist[[s]],
                                   threshold.min = df[[j, i]][1],
                                   threshold.max = df[[j, i]][2])

        list.maps.ndays.Y[[s]][[m]][[cat]] <- gridArithmetics(climatology(thi.pred.Y$ndays), climatology(thi.hist.Y$ndays),
                                                              operator = "-")
        tit.days[[s]][[m]][[cat]] <- sprintf("%s - %s", toupper(m), toupper(cat))
      }

      # Diferencias entre metodos (solo cuando ya estan disponibles)
      if (all(c("qdm", "mbcn") %in% names(list.maps.ndays.Y[[s]]))) {
        res <- compute_ndays_diff(list.maps.ndays.Y[[s]], "qdm", rownames(df))
        list.maps.ndays.diff[[s]] <- res$diffs
        tit.days.diff[[s]]        <- res$titles
      }

      if (all(c("eqm", "mbcn") %in% names(list.maps.ndays.Y[[s]]))) {
        res <- compute_ndays_diff(list.maps.ndays.Y[[s]], "eqm", rownames(df))
        list.maps.ndays.diff.eqm[[s]] <- res$diffs
        tit.days.diff.eqm[[s]]        <- res$titles
      }
    }
  }

  # --- Guardado del RDS para este modelo -------------------------------------
  results <- list(
    maps.dc                  = maps.dc,
    list.maps.ndays.Y        = list.maps.ndays.Y,
    tit.days                 = tit.days,
    list.maps.ndays.diff     = list.maps.ndays.diff,
    tit.days.diff            = tit.days.diff,
    list.maps.ndays.diff.eqm = list.maps.ndays.diff.eqm,
    tit.days.diff.eqm        = tit.days.diff.eqm,
    meta = list(
      model       = model,
      per         = per,
      per_hist    = per_hist,
      season      = season_to_compute,
      months      = months_to_compute,
      species.lab = species.lab,
      stress.cats = rownames(df),
      methods     = c("rcm", met),
      gwl         = gwl_to_compute
    )
  )

  out_file <- file.path(dir_gwl(model), sprintf("THI_deltas_results_%s_%s.rds", gwl_to_compute, season_to_compute))
  saveRDS(results, out_file)
  message(sprintf("  -> Guardado: %s", out_file))

  # Limpiar objetos grandes antes del siguiente modelo
  rm(grid_rcp, grid_rcm, methods.pred, methods.hist,
     maps.dc, list.maps.ndays.Y, tit.days,
     list.maps.ndays.diff, tit.days.diff,
     list.maps.ndays.diff.eqm, tit.days.diff.eqm,
     results)
  gc()
}

message("\nTodos los modelos procesados.")
