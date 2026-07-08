################################################################################
# 01_bias_correction.R
#
# Aplica la corrección de sesgo (validación histórica con cross-validation y
# corrección futura/GWL) a las simulaciones RCM para un modelo y GWL concretos.
#
# Entradas esperadas:
#   Data/grid_obs.rda
#   Data/<model>/grid_rcm.rda
#   Data/<model>/grid_rcp.rda
#
# Salidas:
#   Data/<model>/hist/<metodo>_hist.rds   (independiente del GWL)
#   Data/<model>/<gwl>/<metodo>_pred.rds
################################################################################


################################################################################
# 1. PARÁMETROS
################################################################################

# Modelo a procesar: "NorESM", "HadGEM" o "MPI"
model <- "NorESM"

# GWL a procesar: "GWL1.5", "GWL2" o "GWL3"
gwl <- "GWL3"

# Periodo histórico de referencia (observado y simulado)
years.hist <- 1986:2005

# Periodo futuro correspondiente a cada GWL y modelo.
gwl_periods <- list(
  HadGEM = list("GWL1.5" = 2014:2033, "GWL2" = 2026:2045, "GWL3" = 2045:2064),
  MPI    = list("GWL1.5" = 2008:2027, "GWL2" = 2028:2047, "GWL3" = 2052:2071),
  NorESM = list("GWL1.5" = 2023:2042, "GWL2" = 2039:2058, "GWL3" = 2063:2082)
)

stopifnot(model %in% names(gwl_periods))
stopifnot(gwl %in% c("GWL1.5", "GWL2", "GWL3"))

years.pred <- gwl_periods[[model]][[gwl]]

# Métodos de corrección de sesgo a aplicar
methods <- c("eqm", "qdm", "mbcn")

# Folds para la validación cruzada histórica (NULL = sin cross-validation)
folds <- 4

# Iteraciones para MBCn (y otros métodos MBC)
iter <- 10

# Qué partes ejecutar
run.hist <- TRUE
run.pred <- TRUE

# Directorio principal del proyecto
dir.project <- "/lustre/gmeteo/WORK/uribei/TFM"


################################################################################
# 2. DIRECTORIOS
################################################################################

dir       <- file.path(dir.project, "Data")
diro      <- file.path(dir, model)
diro.hist <- file.path(diro, "hist")
diro.bc   <- file.path(diro, gwl)

dir.create(diro.hist, recursive = TRUE, showWarnings = FALSE)
dir.create(diro.bc, recursive = TRUE, showWarnings = FALSE)


################################################################################
# 3. FUNCIONES
################################################################################

script.bc <- file.path(dir.project, "Scripts", "funciones_bc.R")
if (!file.exists(script.bc)) {
  script.bc <- file.path(dir.project, "R", "funciones_bc.R")
}
source(script.bc)


################################################################################
# 4. DATOS DE ENTRADA
################################################################################

load(file.path(dir, "grid_obs.rda"))    # grid.obs
load(file.path(diro, "grid_rcm.rda"))   # grid.rcm
load(file.path(diro, "grid_rcp.rda"))   # grid.rcp


################################################################################
# 5. VALIDACIÓN HISTÓRICA (con cross-validation)
################################################################################

if (run.hist) {

  met.hist <- methods_bc(
    grid_obs   = grid.obs,
    grid_rcm   = grid.rcm,
    grid_rcp   = grid.rcm,
    methods    = methods,
    years_hist = years.hist,
    years_pred = years.hist,
    folds      = folds,
    iter       = iter
  )

  save_bc_outputs(met.hist, diro = diro.hist, suffix = "hist")

  rm(met.hist)
  gc()
}


################################################################################
# 6. CORRECCIÓN FUTURA / GWL
################################################################################

if (run.pred) {

  met.pred <- methods_bc(
    grid_obs   = grid.obs,
    grid_rcm   = grid.rcm,
    grid_rcp   = grid.rcp,
    methods    = methods,
    years_hist = years.hist,
    years_pred = years.pred,
    folds      = NULL,
    iter       = iter
  )

  save_bc_outputs(met.pred, diro = diro.bc, suffix = "pred")

  rm(met.pred)
  gc()
}


################################################################################
# 7. COMPROBACIÓN DE SALIDAS
################################################################################

list.files(diro.hist, pattern = "\\.rds$", full.names = TRUE)
list.files(diro.bc, pattern = "\\.rds$", full.names = TRUE)
