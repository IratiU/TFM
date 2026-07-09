# =============================================================================
# 03_plot_THI_deltas.R
# Lee los resultados calculados en 04_compute_THI_deltas.R y genera un PNG
# por especie con: (1) mapas principales, (2) diff MBCn-QDM, (3) diff MBCn-EQM.
# =============================================================================

library(visualizeR)
library(transformeR)
library(RColorBrewer)
library(gridExtra)
library(grid)

# -----------------------------------------------------------------------------
# CONFIGURACION
# -----------------------------------------------------------------------------

model_to_plot  <- "NorESM"   # "HadGEM", "MPI", "NorESM"
gwl_to_plot    <- "GWL3"     # "GWL1.5", "GWL2", "GWL3"
season_to_plot <- "JJA"      # "ANN", "DJF", "MAM", "JJA", "SON"

dir.project <- "/lustre/gmeteo/WORK/uribei/TFM"

dir_in  <- file.path(dir.project, "Data", model_to_plot, gwl_to_plot)
dir_out <- file.path(dir_in, season_to_plot)

dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)


# -----------------------------------------------------------------------------
# MASCARA DE TIERRA (a partir de las observaciones ERA5)
# -----------------------------------------------------------------------------

load(file.path(dir.project, "Data", "tas_obs.rda"))
load(file.path(dir.project, "Data", "hurs_obs.rda"))

tas.era5$Variable$varName <- "tas"

grid_obs <- makeMultiGrid("tas" = tas.era5, "hurs" = hurs.era5)
grid_obs <- subsetGrid(grid_obs, years = 1986:2005)

land.mask <- subsetGrid(grid_obs, var = "tas")
land.mask <- climatology(land.mask)

land.mask$Data[!is.na(land.mask$Data)] <- 1
land.mask$Data[is.na(land.mask$Data)]  <- NA


# -----------------------------------------------------------------------------
# CARGA DE RESULTADOS
# -----------------------------------------------------------------------------

message("Cargando resultados...")
res <- readRDS(file.path(dir_in, sprintf("THI_deltas_results_%s_%s.rds", gwl_to_plot, season_to_plot)))

list.maps.ndays.Y        <- res$list.maps.ndays.Y
tit.days                 <- res$tit.days
list.maps.ndays.diff     <- res$list.maps.ndays.diff
tit.days.diff            <- res$tit.days.diff
list.maps.ndays.diff.eqm <- res$list.maps.ndays.diff.eqm
tit.days.diff.eqm        <- res$tit.days.diff.eqm

per         <- res$meta$per
per_hist    <- res$meta$per_hist
species.lab <- res$meta$species.lab
stress.cats <- res$meta$stress.cats
methods     <- res$meta$methods   # c("rcm","eqm","qdm","mbcn")
model       <- res$meta$model
gwl         <- res$meta$gwl
season      <- res$meta$season

season_title <- ifelse(season == "ANN", "annual", season)

# Paleta divergente comun
pal.div <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(30))


# -----------------------------------------------------------------------------
# FUNCION: mapa de diferencia entre metodos (MBCn - comparador)
# -----------------------------------------------------------------------------
make_diff_plot <- function(diff_list, cats, method_label, land_mask) {

  maps.diff <- lapply(cats, function(cat) {
    gridArithmetics(diff_list[[cat]], land_mask, operator = "*")
  })
  names(maps.diff) <- cats

  vals.diff <- unlist(lapply(maps.diff, function(x) as.vector(x$Data)))
  lim.diff  <- ceiling(max(abs(range(vals.diff, na.rm = TRUE))) / 5) * 5
  if (!is.finite(lim.diff) || lim.diff == 0) lim.diff <- 5

  spatialPlot(
    makeMultiGrid(maps.diff),
    backdrop.theme = "countries",
    names.attr = sprintf("MBCn-%s | %s", method_label, tools::toTitleCase(cats)),
    as.table = TRUE,
    layout = c(length(cats), 1),
    par.strip.text = list(cex = 0.75),
    col.regions = rev(colorRampPalette(brewer.pal(11, "RdYlGn"))(30)),
    set.min = -lim.diff,
    set.max = lim.diff,
    at = seq(-lim.diff, lim.diff, length.out = 21),
    main = list(sprintf("Difference MBCn - %s", method_label), cex = 0.8),
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
}


# -----------------------------------------------------------------------------
# GENERACION DE PNGs POR ESPECIE
# Formato compacto: bloque principal 3x4 + diferencias debajo
# -----------------------------------------------------------------------------

for (s in names(list.maps.ndays.Y)) {

  message(sprintf("Generando PNG: %s", toupper(s)))

  ## ---------------------------------------------------------------------------
  ## 1) MAPAS PRINCIPALES: una fila por categoria
  ## ---------------------------------------------------------------------------

  p.main.list <- list()

  for (cat in stress.cats) {

    maps.cat <- list()
    tits.cat <- character()

    for (m in methods) {

      key <- paste(m, cat, sep = "_")

      map.tmp <- list.maps.ndays.Y[[s]][[m]][[cat]]
      map.tmp <- gridArithmetics(map.tmp, land.mask, operator = "*")

      maps.cat[[key]] <- map.tmp
      tits.cat <- c(tits.cat, tit.days[[s]][[m]][[cat]])
    }

    vals.cat <- unlist(lapply(maps.cat, function(x) as.vector(x$Data)))
    lim.cat <- ceiling(max(abs(range(vals.cat, na.rm = TRUE))) / 5) * 5

    if (!is.finite(lim.cat) || lim.cat == 0) {
      lim.cat <- 5
    }

    p.main.list[[cat]] <- spatialPlot(
      makeMultiGrid(maps.cat),
      backdrop.theme = "countries",
      names.attr = tits.cat,
      as.table = TRUE,
      layout = c(length(methods), 1),   # 4 columnas x 1 fila
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

  ## ---------------------------------------------------------------------------
  ## 2) DIFERENCIA MBCn - QDM y MBCn - EQM: una fila cada una
  ## ---------------------------------------------------------------------------

  p.diff     <- make_diff_plot(list.maps.ndays.diff[[s]],     stress.cats, "QDM", land.mask)
  p.diff.eqm <- make_diff_plot(list.maps.ndays.diff.eqm[[s]], stress.cats, "EQM", land.mask)

  ## ---------------------------------------------------------------------------
  ## 3) COMPOSICION DINAMICA Y GUARDADO EN PNG
  ## ---------------------------------------------------------------------------

  out_file <- file.path(
    dir_out,
    sprintf(
      "THI_deltas_%s_%s_%s_%s.png",
      model, gwl, season_title, s
    )
  )

  grobs <- list()
  heights <- c()

  # 3 filas principales: mild, moderate, severe
  for (cat in stress.cats) {
    grobs[[length(grobs) + 1]] <- p.main.list[[cat]]
    heights <- c(heights, 1)
  }

  # Fila 4: MBCn - QDM
  grobs[[length(grobs) + 1]] <- p.diff
  heights <- c(heights, 1)

  # Fila 5: MBCn - EQM
  grobs[[length(grobs) + 1]] <- p.diff.eqm
  heights <- c(heights, 1)

  # Como son 5 filas de mapas, controla aqui la altura por fila
  row.height <- 3.1

  png(
    filename = out_file,
    width    = 16,
    height   = row.height * sum(heights),
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
        "%s - %s | %s | Change in number of heat stress days (%d-%d wrt %d-%d)",
        toupper(s),
        model,
        season_title,
        min(per), max(per),
        min(per_hist), max(per_hist)
      ),
      gp = gpar(fontsize = 15, fontface = "bold")
    )
  )

  dev.off()

  message(sprintf("  -> Guardado: %s", out_file))
}
