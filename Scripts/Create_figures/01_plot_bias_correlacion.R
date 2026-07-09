################################################################################
# 01_plot_bias_correlacion.R
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
# 3. LIBRERIAS Y FUNCIONES
################################################################################

library(transformeR)
library(visualizeR)
library(gridExtra)

# Pone la primera letra de una cadena en mayuscula.
first_upper <- function(x) {
  paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
}

# Calcula una escala simetrica para los mapas de sesgo.
calc_at <- function(bias, n = 25, trim = NULL) {
  vars <- names(bias)
  at_list <- lapply(vars, function(v) {
    vals <- unlist(lapply(bias[[v]], function(g) as.vector(g$Data)))
    vals <- vals[is.finite(vals)]

    if (length(vals) == 0) {
      warning(paste("No hay valores finitos para", v))
      return(NULL)
    }

    if (!is.null(trim)) {
      q <- stats::quantile(vals, probs = c(trim, 1 - trim), na.rm = TRUE)
      vals <- vals[vals >= q[1] & vals <= q[2]]
    }

    m <- max(abs(vals), na.rm = TRUE)
    if (m == 0) m <- 1e-6

    seq(-m, m, length.out = n)
  })

  names(at_list) <- vars
  at_list
}

# Genera los mapas de sesgo absoluto y relativo para cada variable.
plot_bias <- function(bias_abs, bias_rel, main = "", diro = NULL, shp = NULL) {

  vars <- names(bias_abs)
  at_abs <- calc_at(bias_abs)
  at_rel <- calc_at(bias_rel)

  shp_layout <- if (is.null(shp)) {
    NULL
  } else {
    list(
      list("sp.polygons", shp, col = "black", lwd = 1, fill = "transparent")
    )
  }

  plots <- lapply(vars, function(v) {
    plot_abs <- spatialPlot(
      do.call(makeMultiGrid, bias_abs[[v]]),
      as.table = TRUE,
      rev.colors = TRUE,
      backdrop.theme = if (is.null(shp)) "countries" else "none",
      layout = c(length(bias_abs[[v]]), 1),
      at = at_abs[[v]],
      names.attr = names(bias_abs[[v]]),
      main = paste0(first_upper(main), " Abs ", first_upper(v)),
      sp.layout = shp_layout
    )

    plot_rel <- spatialPlot(
      do.call(makeMultiGrid, bias_rel[[v]]),
      as.table = TRUE,
      rev.colors = TRUE,
      backdrop.theme = if (is.null(shp)) "countries" else "none",
      layout = c(length(bias_rel[[v]]), 1),
      at = at_rel[[v]],
      names.attr = names(bias_rel[[v]]),
      main = paste0(first_upper(main), " Rel ", first_upper(v)),
      sp.layout = shp_layout
    )

    p <- grid.arrange(grobs = list(plot_abs, plot_rel))

    if (!is.null(diro)) {
      dir.create(diro, recursive = TRUE, showWarnings = FALSE)

      png(
        filename = file.path(diro, paste0("bias_", main, "_", v, ".png")),
        width = 1800,
        height = 1200,
        res = 100
      )
      grid::grid.draw(p)
      dev.off()
    }

    p
  })

  names(plots) <- vars
  invisible(plots)
}

# Genera las comparativas (scatter) de correlacion observada vs. cada metodo.
comparar_correlaciones <- function(cors, correlation_type = "", diro = NULL) {

  pares <- names(cors)
  methods <- setdiff(names(cors[[1]]), "obs")
  n_methods <- length(methods)
  ncol_plot <- 2
  nrow_plot <- ceiling(n_methods / ncol_plot)

  for (v in pares) {
    if (!is.null(diro)) {
      dir.create(diro, recursive = TRUE, showWarnings = FALSE)

      png(
        filename = file.path(diro, paste0("correlaciones_", correlation_type, "_", v, ".png")),
        width = 1800,
        height = 1200,
        res = 150
      )
    }

    par(mfrow = c(nrow_plot, ncol_plot), oma = c(0, 0, 3, 0))

    for (m in methods) {
      x <- as.vector(cors[[v]]$obs)
      y <- as.vector(cors[[v]][[m]])

      plot(
        x, y,
        main = toupper(m),
        xlab = "Correlacion observada",
        ylab = paste0("Correlacion ", toupper(m)),
        pch = 16,
        col = "black",
        cex = 0.8,
        las = 1,
        bty = "l"
      )

      abline(a = 0, b = 1, col = "red", lwd = 1, lty = 2)
    }

    mtext(
      paste0("Comparativa de correlaciones ", correlation_type, " - ", v),
      outer = TRUE,
      cex = 1.3,
      font = 2,
      line = 1
    )

    if (!is.null(diro)) {
      dev.off()
    }
  }
}


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