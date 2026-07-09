
# =============================================================================
# 02_plot_THI_seasons.R
#
# Genera los mapas espaciales (PNG) del nº de dias bajo estres termico (THI)
# por estacion del anio y categoria de severidad, para cada especie, a partir
# del RDS generado por compute_THI_estaciones.R
# =============================================================================
 
# ---- Librerias --------------------------------------------------------------
library(transformeR)
library(visualizeR)
library(sp)
library(RColorBrewer)
library(gridExtra)
library(grid)
 
# ---- Parametros ---------------------------------------------------------------
dir.project <- "/lustre/gmeteo/WORK/uribei/TFM"
dir.obs     <- file.path(dir.project, "Data/obs")
 
per_hist <- 1986:2005
 
f.in <- file.path(dir.obs, sprintf("THI_obs_ndays_by_season_%d-%d.rds",
                                    min(per_hist), max(per_hist)))
 
at.list <- list(
  mild     = seq(0, 70, length.out = 21),
  moderate = seq(0, 25, length.out = 21),
  severe   = seq(0, 10, length.out = 21)
)
 
# ---- Comprobaciones -----------------------------------------------------------
stopifnot(file.exists(f.in))
 
# ---- Carga de resultados --------------------------------------------------------
obs_THI     <- readRDS(f.in)
maps.by.cat <- obs_THI$maps.by.cat
species.lab <- obs_THI$species.lab
categories  <- obs_THI$categories
per_hist    <- obs_THI$per_hist
 
# ---- Figuras ----------------------------------------------------------------------
for (sp in species.lab) {
 
  plots <- list()
 
  for (j in seq_along(categories)) {
 
    category <- categories[j]
    maps.cat <- maps.by.cat[[sp]][[category]]
    at.cat   <- at.list[[category]]
 
    cols <- colorRampPalette(
      c("white", brewer.pal(9, "YlOrRd"))
    )(length(at.cat) - 1)
 
    season.names <- names(maps.cat)
 
    plots[[j]] <- spatialPlot(
      climatology(makeMultiGrid(maps.cat, skip.temporal.check = TRUE)),
      names.attr     = toupper(season.names),
      main           = toupper(category),
      as.table       = TRUE,
      layout         = c(1, length(season.names)),
      backdrop.theme = "countries",
      col.regions    = cols,
      at             = at.cat,
      set.min        = min(at.cat),
      set.max        = max(at.cat)
    )
  }
 
  f.png <- file.path(dir.obs, sprintf("THI_obs_seasons-%s.png", sp))
 
  png(
    file   = f.png,
    width  = 14,
    height = 9,
    units  = "in",
    res    = 300,
    type   = "cairo"
  )
 
  grid.arrange(
    plots[[1]], plots[[2]], plots[[3]],
    ncol = 3,
    top = textGrob(
      sprintf(
        "%s -Number of days under heat stress by season (%d-%d)", toupper(sp),
        min(per_hist), max(per_hist)
      ),
      gp = gpar(fontsize = 16, fontface = "bold")
    )
  )
 
  dev.off()
 
  message("Guardado: ", f.png)
}