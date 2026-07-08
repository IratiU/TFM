################################################################################
# 03_plot_bias_correlacion.R
#
# Genera los mapas de sesgo (absoluto/relativo) y las comparativas de
# correlación calculados en 02_calcular_bias_correlacion.R.
################################################################################


################################################################################
# 1. PARÁMETROS
################################################################################

model <- "NorESM"

# Debe coincidir con el correlation_type usado en 02_calcular_bias_correlacion.R
correlation_type <- "Pearson"

# Shapefile opcional para superponer contornos en los mapas (NULL = sin shapefile)
shp <- NULL

dir.project <- "/lustre/gmeteo/WORK/uribei/TFM"


################################################################################
# 2. DIRECTORIOS
################################################################################

dir       <- file.path(dir.project, "Data")
diro      <- file.path(dir, model)
diro.hist <- file.path(diro, "hist")
diro.fig  <- file.path(diro.hist, "figuras")

dir.create(diro.fig, recursive = TRUE, showWarnings = FALSE)


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

bias_abs <- readRDS(file.path(diro.hist, "bias_abs.rds"))
bias_rel <- readRDS(file.path(diro.hist, "bias_rel.rds"))
correlaciones <- readRDS(
  file.path(diro.hist, paste0("correlaciones_", tolower(correlation_type), ".rds"))
)


################################################################################
# 5. MAPAS DE SESGO
################################################################################

plot_bias(
  bias_abs = bias_abs,
  bias_rel = bias_rel,
  main     = paste0(model, "_hist"),
  diro     = diro.fig,
  shp      = shp
)


################################################################################
# 6. COMPARATIVA DE CORRELACIONES
################################################################################

comparar_correlaciones(
  cors             = correlaciones,
  correlation_type = correlation_type,
  diro             = diro.fig
)


################################################################################
# 7. COMPROBACIÓN DE SALIDAS
################################################################################

list.files(diro.fig, pattern = "\\.png$", full.names = TRUE)
