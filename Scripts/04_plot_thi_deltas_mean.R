# =============================================================================
# 04_plot_THI_deltas.R
# Lee los resultados calculados en 01_compute_THI_deltas.R y genera un PDF
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
model_to_plot  <- "Simulation Mean"
gwl_to_plot    <- "GWL1.5"     # "GWL1.5", "GWL2", "GWL3"
season_to_plot <- "JJA"      # "ANN", "DJF", "MAM", "JJA", "SON"

dir_in <- file.path(
  "~/lustre/gmeteo/WORK/uribei/TFM/Data/ENSEMBLE_MEAN",
  gwl_to_plot
)

source("~/lustre/gmeteo/WORK/uribei/TFM/Scripts/THI/04_aggrement_points.R")


dir_out <- file.path(
  "~/lustre/gmeteo/WORK/uribei/TFM/Data/ENSEMBLE_MEAN",
  gwl_to_plot,
  season_to_plot
)

dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# MASCARA ESPACIAL
# -----------------------------------------------------------------------------
load("~/lustre/gmeteo/WORK/uribei/TFM/Data/tas_obs.rda")
load("~/lustre/gmeteo/WORK/uribei/TFM/Data/hurs_obs.rda")

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

res <- readRDS(
  file.path(
    dir_in,
    sprintf(
      "THI_ENSEMBLE_MEAN_%s_%s.rds",
      gwl_to_plot,
      season_to_plot
    )
  )
)

list.maps.ndays.Y        <- res$list.maps.ndays.Y
tit.days                 <- res$tit.days
list.maps.ndays.diff     <- res$list.maps.ndays.diff
tit.days.diff            <- res$tit.days.diff
list.maps.ndays.diff.eqm <- res$list.maps.ndays.diff.eqm
tit.days.diff.eqm        <- res$tit.days.diff.eqm

# -----------------------------------------------------------------------------
# IMPORTANTE:
# No usar res$meta$methods ni res$meta$stress.cats para el ensemble mean.
# Los fijamos directamente para evitar listas vacias.
# -----------------------------------------------------------------------------
species.lab <- names(list.maps.ndays.Y)
methods     <- c("rcm", "eqm", "qdm", "mbcn")
stress.cats <- c("mild", "moderate", "severe")

season_title <- ifelse(season_to_plot == "ANN", "annual", season_to_plot)

model  <- model_to_plot
gwl    <- gwl_to_plot
season <- season_to_plot

pal.div <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(30))
pal.diff <- rev(colorRampPalette(brewer.pal(11, "RdYlGn"))(30))

# -----------------------------------------------------------------------------
# FUNCIONES AUXILIARES
# -----------------------------------------------------------------------------

get_lim <- function(maps.list, step = 5, default = 5) {
  
  vals <- unlist(lapply(maps.list, function(x) as.vector(x$Data)))
  vals <- vals[is.finite(vals)]
  
  if (length(vals) == 0) {
    return(default)
  }
  
  lim <- ceiling(max(abs(vals), na.rm = TRUE) / step) * step
  
  if (!is.finite(lim) || lim == 0) {
    lim <- default
  }
  
  return(lim)
}

get_agreement_layout <- function(gwl, sp, cat, methods) {
  sp.layout.cat <- list()
  for (i in seq_along(methods)) {
    m <- methods[i]
    pts <- agreement_pts[[gwl]][[sp]][[m]][[cat]]
    pts[["which"]] <- i
    sp.layout.cat[[length(sp.layout.cat) + 1]] <- pts
  }
  
  return(sp.layout.cat)
}
at.list <- list(mild = 60, moderate = 50, severe = 50)
# -----------------------------------------------------------------------------
# GENERACION DE PNGs POR ESPECIE
# -----------------------------------------------------------------------------

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
      
      # Titulos seguros, sin depender de que tit.days este perfecto
      tits.cat <- c(
        tits.cat,
        sprintf("%s | %s", toupper(m), tools::toTitleCase(cat))
      )
    }
    
    lim.cat <- at.list[[cat]]
    
    sp.layout.cat <- get_agreement_layout(
      gwl = gwl_to_plot,
      sp = s,
      cat = cat,
      methods = methods
    )
    
    p.main.list[[cat]] <- spatialPlot(
      makeMultiGrid(maps.cat),
      backdrop.theme = "countries",
      names.attr = tits.cat,
      as.table = TRUE,
      layout = c(length(methods), 1),
      par.strip.text = list(cex = 1),
      col.regions = pal.div,
      set.min = -lim.cat,
      set.max = lim.cat,
      at = seq(-lim.cat, lim.cat, length.out = 21),
      xlab = "",
      ylab = "",
      sp.layout = sp.layout.cat,
      colorkey = list(
        space = "right",
        labels = list(cex = 1)
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
    par.strip.text = list(cex = 1),
    col.regions = pal.diff,
    set.min = -20,
    set.max = 20,
    at = seq(-20, 20, length.out = 21),
    main = list("Difference MBCn - QDM", cex = 1.2),
    xlab = "",
    ylab = "",
    colorkey = list(
      space = "right",
      labels = list(cex = 1)
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
    par.strip.text = list(cex = 1),
    col.regions = pal.diff, set.max = 20, set.min = -20,
    at = seq(-20, 20, length.out = 21),
    main = list("Difference MBCn - EQM", cex = 1.2),
    xlab = "",
    ylab = "",
    colorkey = list(
      space = "right",
      labels = list(cex = 1)
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
      gsub(" ", "_", model),
      gwl,
      season_title,
      s
    )
  )
  
  grobs <- list()
  heights <- c()
  
  for (cat in stress.cats) {
    grobs[[length(grobs) + 1]] <- p.main.list[[cat]]
    heights <- c(heights, 1)
  }
  
  grobs[[length(grobs) + 1]] <- p.diff
  heights <- c(heights, 1)
  
  grobs[[length(grobs) + 1]] <- p.diff.eqm
  heights <- c(heights, 1)
  
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
        model,
        gwl
      ),
      gp = gpar(fontsize = 15, fontface = "bold")
    )
  )
  
  dev.off()
  
  message(sprintf("  -> Guardado: %s", out_file))
}