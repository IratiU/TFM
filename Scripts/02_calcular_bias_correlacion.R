################################################################################
# 02_calcular_bias_correlacion.R
#
# Calcula el sesgo (absoluto y relativo) y la correlación entre variables para
# cada método de corrección de sesgo, comparando frente a las observaciones en
# el periodo histórico. No genera gráficos: los resultados se guardan en RDS
# para que 03_plot_bias_correlacion.R los pinte.
#
# Requiere que 01_bias_correction.R ya haya generado:
#   Data/<model>/hist/<metodo>_hist.rds
################################################################################


################################################################################
# 1. PARÁMETROS
################################################################################

model <- "NorESM"

# Periodo de validación (debe coincidir con years.hist de 01_bias_correction.R)
years.hist <- 1986:2005

# Métodos guardados en 01_bias_correction.R que se quieren evaluar
methods <- c("eqm", "qdm", "mbcn")

# Incluir la simulación RCM sin corregir como referencia ("raw")
incluir_raw <- TRUE

# Medidas de sesgo (revisa los measure.code/index.code disponibles en
# climate4R.value con ?valueMeasure; "biasrel" es el nombre habitual para el
# sesgo relativo, confírmalo en tu instalación)
measure.code.abs <- "bias"
measure.code.rel <- "biasrel"
index.code.bias  <- "mean"

# Medida de correlación entre variables
measure.code.cor <- "corr"
index.code.cor   <- NULL
correlation_type <- "Pearson"

# Directorio principal del proyecto
dir.project <- "/lustre/gmeteo/WORK/uribei/TFM"


################################################################################
# 2. DIRECTORIOS
################################################################################

dir       <- file.path(dir.project, "Data")
diro      <- file.path(dir, model)
diro.hist <- file.path(diro, "hist")


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
if (incluir_raw) {
  load(file.path(diro, "grid_rcm.rda")) # grid.rcm
}

# Corrección de sesgo histórica (validación cruzada) generada en 01_bias_correction.R
met.hist <- setNames(
  lapply(methods, function(m) {
    readRDS(file.path(diro.hist, paste0(m, "_hist.rds")))
  }),
  methods
)

if (incluir_raw) {
  met.hist[["raw"]] <- grid.rcm
}


################################################################################
# 5. CÁLCULO DE SESGO (ABSOLUTO Y RELATIVO)
################################################################################

bias_abs <- calcular_bias(
  methods      = met.hist,
  grid_obs     = grid.obs,
  years        = years.hist,
  measure.code = measure.code.abs,
  index.code   = index.code.bias
)

bias_rel <- calcular_bias(
  methods      = met.hist,
  grid_obs     = grid.obs,
  years        = years.hist,
  measure.code = measure.code.rel,
  index.code   = index.code.bias
)

saveRDS(bias_abs, file.path(diro.hist, "bias_abs.rds"))
saveRDS(bias_rel, file.path(diro.hist, "bias_rel.rds"))


################################################################################
# 6. CÁLCULO DE CORRELACIÓN ENTRE VARIABLES
################################################################################

# calcular_correlacion() dibuja un boxplot por cada par de variables como
# comprobación rápida; en modo batch se redirige a un dispositivo nulo para no
# depender de una pantalla gráfica ni generar un Rplots.pdf de sobra.
pdf(NULL)
correlaciones <- calcular_correlacion(
  methods          = met.hist,
  grid_obs         = grid.obs,
  years_hist       = years.hist,
  measure.code     = measure.code.cor,
  index.code       = index.code.cor,
  correlation_type = correlation_type
)
dev.off()

saveRDS(
  correlaciones,
  file.path(diro.hist, paste0("correlaciones_", tolower(correlation_type), ".rds"))
)


################################################################################
# 7. COMPROBACIÓN DE SALIDAS
################################################################################

list.files(diro.hist, pattern = "^(bias_|correlaciones_).*\\.rds$", full.names = TRUE)
