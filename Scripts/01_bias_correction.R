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
# 3. LIBRERIAS Y FUNCIONES
################################################################################

library(transformeR)
library(downscaleR)

# Transformacion logit para variables acotadas entre 0 y 100, como hurs.
to_logit <- function(grid, eps = 1e-5) {
  x <- grid$Data / 100
  x <- pmin(pmax(x, eps), 1 - eps)
  grid$Data <- qlogis(x)
  grid
}

# Transformacion inversa del logit.
from_logit <- function(grid) {
  grid$Data <- plogis(grid$Data) * 100
  grid
}

# Guarda cada metodo de correccion de sesgo en un archivo .rds.
save_bc_outputs <- function(met, diro, suffix) {
  dir.create(diro, recursive = TRUE, showWarnings = FALSE)

  for (m in names(met)) {
    saveRDS(
      met[[m]],
      file = file.path(diro, paste0(m, "_", suffix, ".rds"))
    )
  }
}

# Aplica los metodos de correccion de sesgo (univariados y multivariados)
# a un modelo, para un periodo historico y un periodo futuro/GWL dados.
methods_bc <- function(grid_obs, grid_rcm, grid_rcp,
                       methods = c("eqm", "mbcn", "mbcr", "mbcp", "qdm"),
                       years_hist, years_pred,
                       LatLim = NULL, LonLim = NULL,
                       folds = NULL, iter = 20, diro = NULL) {

  vars0 <- getVarNames(grid_obs)
  vars <- tolower(vars0)
  methods <- tolower(methods)

  res <- list()

  # Datos historicos simulados.
  x_list <- lapply(vars, function(v) {
    subsetGrid(grid_rcm, var = v, years = years_hist)
  })
  names(x_list) <- vars

  # Observaciones historicas.
  y_list <- lapply(vars, function(v) {
    subsetGrid(grid_obs, var = v, years = years_hist)
  })
  names(y_list) <- vars

  # Datos a corregir, normalmente futuro/proyeccion.
  newdata <- lapply(vars, function(v) {
    subsetGrid(grid_rcp, var = v, years = years_pred)
  })
  names(newdata) <- vars

  # Metodos univariados.
  for (m in methods) {
    if (m %in% c("eqm", "qdm")) {
      corrected <- lapply(vars, function(v) {
        biasCorrection(
          y = y_list[[v]],
          x = x_list[[v]],
          newdata = newdata[[v]],
          wet.threshold = if (v == "pr") 0.1 else NULL,
          precipitation = (v == "pr"),
          method = m,
          cross.val = if (!is.null(folds)) "kfold" else NULL,
          folds = folds
        )
      })

      names(corrected) <- vars
      res[[m]] <- makeMultiGrid(corrected)
    }
  }

  # Metodos multivariados.
  mbc <- intersect(methods, c("mbcp", "mbcr", "mbcn"))

  for (m in mbc) {
    rot.seq <- replicate(iter, rot.random(length(vars)), simplify = FALSE)
    ratio.seq <- vars %in% c("pr", "windsfc")

    mbc.args <- if (m == "mbcn") {
      list(
        iter = iter,
        rot.seq = rot.seq,
        ratio.seq = ratio.seq,
        trace = 0.1
      )
    } else {
      list(
        iter = iter,
        ratio.seq = ratio.seq,
        trace = 0.1
      )
    }

    y_use <- y_list
    x_use <- x_list
    newdata_use <- newdata

    # Para MBCn, la humedad relativa se transforma a la escala real mediante logit.
    if (m == "mbcn" && "hurs" %in% vars) {
      y_use[["hurs"]] <- to_logit(y_use[["hurs"]])
      x_use[["hurs"]] <- to_logit(x_use[["hurs"]])
      newdata_use[["hurs"]] <- to_logit(newdata_use[["hurs"]])
    }

    bc <- biasCorrection(
      y = y_use,
      x = x_use,
      newdata = newdata_use,
      wet.threshold = if ("pr" %in% vars) 0.1 else NULL,
      precipitation = ("pr" %in% vars),
      mbc.args = mbc.args,
      method = m,
      cross.val = if (!is.null(folds)) "kfold" else NULL,
      folds = folds
    )

    names(bc) <- vars

    # Se vuelve a pasar hurs a porcentaje.
    if (m == "mbcn" && "hurs" %in% vars) {
      bc[["hurs"]] <- from_logit(bc[["hurs"]])
    }

    res[[m]] <- makeMultiGrid(bc)
  }

  res
}


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