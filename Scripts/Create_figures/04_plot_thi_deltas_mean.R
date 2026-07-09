# =============================================================================
# 04_plot_thi_deltas_mean.R
# Lee el RDS generado por 06_compute_THI_deltas_ensemble_mean.R (ensemble mean
# + puntos de acuerdo entre GCMs) y genera un PNG por especie con:
#   (1) mapas principales (uno por metodo x categoria)
#   (2) diff MBCn - QDM
#   (3) diff MBCn - EQM
# =============================================================================

# ---- Librerias ------------------------------------------------------------------
library(visualizeR)
library(transformeR)
library(RColorBrewer)
library(gridExtra)
library(grid)

# ---- Parametros -------------------------------------------------------------------
model_to_plot  <- "Simulation Mean"
gwl_to_plot    <- "GWL3"     # "GWL1.5", "GWL2", "GWL3"
season_to_plot <- "JJA"      # "ANN", "DJF", "MAM", "JJA", "SON"

methods     <- c("rcm", "eqm", "qdm", "mbcn")
stress.cats <- c("mild", "moderate", "severe")

at.list <- list(mild = 60, moderate = 50, severe = 50)

pal.div  <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(30))
pal.diff <- rev(colorRampPalette(brewer.pal(11, "RdYlGn"))(30))

dir_project <- "/lustre/gmeteo/WORK/uribei/TFM"

dir_in  <- file.path(dir_project, "Data/ENSEMBLE_MEAN", gwl_to_plot)
dir_out <- file.path(dir_project, "Data/ENSEMBLE_MEAN", gwl_to_plot, season_to_plot)
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

f.in <- file.path(dir_in, sprintf("THI_ENSEMBLE_MEAN_%s_%s.rds", gwl_to_plot, season_to_plot))

f.tas  <- file.path(dir_project, "Data/tas_obs.rda")
f.hurs <- file.path(dir_project, "Data/hurs_obs.rda")

# ---- Comprobaciones ----------------------------------------------------------------
stopifnot(
  file.exists(f.in),
  file.exists(f.tas),
  file.exists(f.hurs)
)

# ---- Funciones auxiliares ------------------------------------------------------------

get_lim <- function(maps.list, step = 5, default = 5) {

  vals <- unlist(lapply(maps.list, function(x) as.vector(x$Data)))
  vals <- vals[is.finite(vals)]

  if (length(vals) == 0) return(default)

  lim <- ceiling(max(abs(vals), na.rm = TRUE) / step) * step

  if (!is.finite(lim) || lim == 0) lim <- default

  lim
}

# ---- Mascara espacial ----------------------------------------------------------------
load(f.tas)
load(f.hurs)

tas.era5$Variable$varName <- "tas"

grid_obs <- makeMultiGrid("tas" = tas.era5, "hurs" = hurs.era5)
grid_obs <- subsetGrid(grid_obs, years = 1986:2005)

land.mask <- subsetGrid(grid_obs, var = "tas")
land.mask <- climatology(land.mask)

land.mask$Data[!is.na(land.mask$Data)] <- 1
land.mask$Data[is.na(land.mask$Data)]  <- NA

rm(tas.era5, hurs.era5, grid_obs)
gc()

# ---- Carga de resultados --------------------------------------------------------------
message("Cargando resultados...")

res <- readRDS(f.in)

list.maps.ndays.Y        <- res$list.maps.ndays.Y
list.maps.ndays.diff     <- res$list.maps.ndays.diff
list.maps.ndays.diff.eqm <- res$list.maps.ndays.diff.eqm
agreement_pts            <- res$agreement_pts

# IMPORTANTE:
# No usar res$meta$methods ni res$meta$stress.cats para el ensemble mean.
# Los fijamos directamente (arriba, en parametros) para evitar listas vacias.
species.lab <- names(list.maps.ndays.Y)

season_title <- ifelse(season_to_plot == "ANN", "annual", season_to_plot)

# ---- Generacion de PNGs por especie -------------------------------------------------
for (s in species.lab) {

  message(sprintf("Generando PNG: %s", toupper(s)))

  # ---------------------------------------------------------------------------
  # 1) MAPAS PRINCIPALES: una fila por categoria
  # ---------------------------------------------------------------------------

  p.main.list <- list()

  for (cat in stress.cats) {

    maps.cat <- list()
    tits.cat <- character()

    for (m in methods) {

      map.tmp <- list.maps.ndays.Y[[s]][[m]][[cat]]
      map.tmp <- gridArithmetics(map.tmp, land.mask, operator = "*")

      key <- paste(m, cat, sep = "_")
      maps.cat[[key]] <- map.tmp

      tits.cat <- c(
        tits.cat,
        sprintf("%s | %s", toupper(m), tools::toTitleCase(cat))
      )
    }

    lim.cat <- at.list[[cat]]

    p.main.list[[cat]] <- spatialPlot(
      makeMultiGrid(maps.cat),
      backdrop.theme = "countries",
      names.attr = tits.cat,
      as.table = TRUE,
      layout = c(length(methods), 1),
      par.strip.text = list(cex = 0.75),
      col.regions = pal.div,
      set.min = -lim.cat,
      set.max = lim.cat,
      at = seq(-lim.cat, lim.cat, length.out = 21),
      xlab = "",
      ylab = "",
      colorkey = list(
        space = "right",
        labels = list(cex = 0.65)
      ),
      par.settings = list(
        panel.background = list(col = "white"),
        background = list(col = "white"),
        layout.heights = list(
          top.padding = 0,
          bottom.padding = 0,
          main = 1.2,
          main.key.padding = 0.1,
          key.axis.padding = 0,
          axis.xlab.padding = 0,
          xlab = 0
        ),
        layout.widths = list(
          left.padding = 0,
          right.padding = 0,
          ylab = 0,
          axis.ylab.padding = 0
        )
      )
    )
  }

  # ---------------------------------------------------------------------------
  # 2) DIFERENCIA MBCn - QDM
  # ---------------------------------------------------------------------------

  maps.diff <- list()

  for (cat in stress.cats) {
    map.tmp <- list.maps.ndays.diff[[s]][[cat]]
    map.tmp <- gridArithmetics(map.tmp, land.mask, operator = "*")
    maps.diff[[cat]] <- map.tmp
  }

  lim.diff <- get_lim(maps.diff, step = 5, default = 5)

  p.diff <- spatialPlot(
    makeMultiGrid(maps.diff),
    backdrop.theme = "countries",
    names.attr = c(
      "MBCn-QDM | Mild",
      "MBCn-QDM | Moderate",
      "MBCn-QDM | Severe"
    ),
    as.table = TRUE,
    layout = c(3, 1),
    par.strip.text = list(cex = 0.75),
    col.regions = pal.diff,
    set.min = -20,
    set.max = 20,
    at = seq(-20, 20, length.out = 21),
    main = list("Difference MBCn - QDM", cex = 0.8),
    xlab = "",
    ylab = "",
    colorkey = list(
      space = "right",
      labels = list(cex = 0.65)
    ),
    par.settings = list(
      panel.background = list(col = "white"),
      background = list(col = "white"),
      layout.heights = list(
        top.padding = 0,
        bottom.padding = 0,
        main = 1.2,
        main.key.padding = 0.5,
        key.axis.padding = 0,
        axis.xlab.padding = 0,
        xlab = 0
      ),
      layout.widths = list(
        left.padding = 0,
        right.padding = 0,
        ylab = 0,
        axis.ylab.padding = 0
      )
    )
  )

  # ---------------------------------------------------------------------------
  # 3) DIFERENCIA MBCn - EQM
  # ---------------------------------------------------------------------------

  maps.diff.eqm <- list()

  for (cat in stress.cats) {
    map.tmp <- list.maps.ndays.diff.eqm[[s]][[cat]]
    map.tmp <- gridArithmetics(map.tmp, land.mask, operator = "*")
    maps.diff.eqm[[cat]] <- map.tmp
  }

  lim.diff.eqm <- get_lim(maps.diff.eqm, step = 5, default = 5)

  p.diff.eqm <- spatialPlot(
    makeMultiGrid(maps.diff.eqm),
    backdrop.theme = "countries",
    names.attr = c(
      "MBCn-EQM | Mild",
      "MBCn-EQM | Moderate",
      "MBCn-EQM | Severe"
    ),
    as.table = TRUE,
    layout = c(3, 1),
    par.strip.text = list(cex = 0.75),
    col.regions = pal.diff,
    set.max = 20,
    set.min = -20,
    at = seq(-20, 20, length.out = 21),
    main = list("Difference MBCn - EQM", cex = 0.8),
    xlab = "",
    ylab = "",
    colorkey = list(
      space = "right",
      labels = list(cex = 0.65)
    ),
    par.settings = list(
      panel.background = list(col = "white"),
      background = list(col = "white"),
      layout.heights = list(
        top.padding = 0,
        bottom.padding = 0,
        main = 1.2,
        main.key.padding = 0.5,
        key.axis.padding = 0,
        axis.xlab.padding = 0,
        xlab = 0
      ),
      layout.widths = list(
        left.padding = 0,
        right.padding = 0,
        ylab = 0,
        axis.ylab.padding = 0
      )
    )
  )

  # ---------------------------------------------------------------------------
  # 4) COMPOSICION Y GUARDADO
  # ---------------------------------------------------------------------------

  out_file <- file.path(
    dir_out,
    sprintf(
      "THI_deltas_%s_%s_%s_%s.png",
      gsub(" ", "_", model_to_plot),
      gwl_to_plot,
      season_title,
      s
    )
  )

  grobs   <- c(p.main.list[stress.cats], list(p.diff, p.diff.eqm))
  heights <- rep(1, length(grobs))

  png(
    filename = out_file,
    width    = 16,
    height   = 3.1 * sum(heights),
    units    = "in",
    res      = 300,
    type     = "cairo"
  )

  grid.arrange(
    grobs = grobs,
    ncol = 1,
    heights = heights,
    top = textGrob(
      sprintf(
        "%s - %s | Change in number of heat stress days for %s",
        toupper(s),
        model_to_plot,
        gwl_to_plot
      ),
      gp = gpar(fontsize = 15, fontface = "bold")
    )
  )

  dev.off()

  message(sprintf("  -> Guardado: %s", out_file))
}