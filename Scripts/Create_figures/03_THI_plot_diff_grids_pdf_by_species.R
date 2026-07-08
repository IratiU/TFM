library(gridExtra)
library(grid)
library(convertR)
library(transformeR)
library(loadeR)
library(visualizeR)
library(RColorBrewer)
library(lattice)
source("/lustre/gmeteo/WORK/uribei/TFM/Scripts/THI/functions_THI.R")
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
## -----------------------------------------------------------------------------
## Paths
## -----------------------------------------------------------------------------
diro <- "/lustre/gmeteo/WORK/uribei/TFM/Data"
grid.out.dir <- file.path(diro, "THI_saved_grids")
plot.out.dir <- file.path(diro, "THI_plots_pdf")
dir.create(plot.out.dir, recursive = TRUE, showWarnings = FALSE)

## -----------------------------------------------------------------------------
## Basic configuration
## -----------------------------------------------------------------------------
per <- 1986:2005
season.name <- "JJA"

season.def <- list(DJF = c(12, 1, 2), MAM = c(3, 4, 5), JJA = c(6, 7, 8), SON = c(9, 10, 11))

months <- season.def[[season.name]]


simulations <- c("NorESM", "MPI", "HadGEM")
methods <- c("eqm", "qdm", "mbcr", "mbcp", "mbcn")
method.names <- c("EQM-Obs", "QDM-Obs", "MBCr-Obs", "MBCp-Obs", "MBCn-Obs")

species <- c("cattle", "cattle", "ruminants", "ruminants", "poultry", "swine")
species.lab <- c("dairy", "beef", "sheep", "goats", "poultry", "swine")

mild <- list(c(72, 78), c(75, 81), c(25, 30), c(27, 32), c(27.8, 28.8), c(23.33, 26.11))
moderate <- list(c(78, 88), c(81, 91), c(30, 35), c(32, 37), c(28.8, 30), c(26.11, 28.88))
severe <- list(c(88, Inf), c(91, Inf), c(35, Inf), c(37, Inf), c(30, Inf), c(28.88, Inf))

df <- data.frame(rbind(mild, moderate, severe))
colnames(df) <- species.lab

## Observed maps shown in the first row of each PDF
obs.categories.by.species <- function(sp_lab) {
  if (sp_lab == "swine") return(c("mild", "moderate", "severe"))
  if (sp_lab == "dairy") return(c("mild", "moderate", "severe"))
  return(c("mild"))
}

## Difference blocks shown below the observed row
plot.categories.by.species <- obs.categories.by.species

## Color scale for model minus observation differences
get_at_diff <- function(category) {
  if (category == "mild") return(seq(-10, 10, length.out = 21))
  return(seq(-4, 4, length.out = 21))
}

get_cols_diff <- function(category) {
  at.diff <- get_at_diff(category)
  rev(colorRampPalette(brewer.pal(11, "RdBu"))(length(at.diff) - 1))
}

## Color scale for observed annual number of THI days
get_at_obs <- function(category) {
  if (category == "mild") return(seq(0, 60, length.out = 21))
  if (category == "moderate") return(seq(0, 20, length.out = 21))
  if (category == "severe") return(seq(0, 10, length.out = 21))
  stop("Unknown THI category: ", category)
}

get_cols_obs <- function(category) {
  at.obs <- get_at_obs(category)
  colorRampPalette(brewer.pal(9, "YlOrRd"))(length(at.obs) - 1)
}
## -----------------------------------------------------------------------------
## Load observations and define mask
## -----------------------------------------------------------------------------
load(file.path(diro, "tas_obs.rda"))
load(file.path(diro, "hurs_obs.rda"))
tas.era5$Variable$varName <- "tas"

grid.obs <- makeMultiGrid("tas" = tas.era5, "hurs" = hurs.era5)
grid.obs <- subsetGrid(grid.obs, years = per)

obs.mask <- subsetGrid(grid.obs, var = "tas")
obs.mask <- climatology(obs.mask)
mask.na <- is.na(drop(obs.mask$Data))

## -----------------------------------------------------------------------------
## Auxiliary functions
## -----------------------------------------------------------------------------
compute_seasonal_thi <- function(grid, sp, per, months) {
  
  tas <- subsetGrid(grid, season = months, years = per, var = "tas")
  attributes(tas$Variable)$units <- "degC"
  
  hurs <- subsetGrid(grid, season = months, years = per, var = "hurs")
  
  computeTHI(tas, hurs, species = sp)
}

make_land_mask <- function(template, mask.na) {
  land.mask <- template
  land.mask$Data[] <- 1
  d <- dim(land.mask$Data)
  if (length(d) == 2) land.mask$Data[mask.na] <- NA
  if (length(d) == 3) {
    for (tt in seq_len(d[1])) {
      tmp <- land.mask$Data[tt, , ]
      tmp[mask.na] <- NA
      land.mask$Data[tt, , ] <- tmp
    }
  }
  land.mask
}

apply_obs_mask <- function(grid, mask.na) {
  out <- grid
  d <- dim(out$Data)
  if (length(d) == 2) out$Data[mask.na] <- NA
  if (length(d) == 3) {
    for (tt in seq_len(d[1])) {
      tmp <- out$Data[tt, , ]
      tmp[mask.na] <- NA
      out$Data[tt, , ] <- tmp
    }
  }
  out
}

make_obs_maps <- function(thi.obs, sp_lab, categories) {
  maps.obs <- list()
  for (cat in categories) {
    thi.Y.obs <- analyzeTHI_Y(
      thi = thi.obs,
      threshold.min = df[[sp_lab]][[cat]][1],
      threshold.max = df[[sp_lab]][[cat]][2]
    )
    maps.obs[[cat]] <- apply_obs_mask(thi.Y.obs$ndays, mask.na)
  }
  maps.obs
}
plot_obs_panel <- function(obs.map, sp_lab, category) {
  
  at.obs <- get_at_obs(category)
  cols.obs <- get_cols_obs(category)
  
  spatialPlot(
    climatology(obs.map),
    names.attr = paste("OBS", toupper(category), sep = " - "),
    as.table = TRUE,
    layout = c(1, 1),
    backdrop.theme = "countries",
    col.regions = cols.obs,
    set.min = min(at.obs),
    set.max = max(at.obs),
    at = at.obs,
    main = list(paste0("OBS - ", toupper(category)), cex = 1),
    par.strip.text = list(cex = 1),
    par.settings = list(
      panel.background = list(col = "white"),
      background = list(col = "white")
    )
  )
}

plot_category_block <- function(maps.diff.by.sim, obs.map, sp_lab, category) {
  
  obs.panel <- plot_obs_panel(obs.map, sp_lab, category)
  diff.block <- plot_diff_block(maps.diff.by.sim, category)
  
  arrangeGrob(
    obs.panel,
    diff.block,
    ncol = 2,
    widths = c(1.2, 4.2)
  )
}
read_diff_file <- function(category, sp_lab, season.name) {
  saved.file <- file.path(grid.out.dir, paste0("THI_diff_grids_", toupper(category), "_", season.name, "_", sp_lab, ".rds"))
  if (!file.exists(saved.file)) {
    warning("File not found: ", saved.file)
    return(NULL)
  }
  message("Loading: ", saved.file)
  readRDS(saved.file)
}
## -----------------------------------------------------------------------------
## Plot difference block: one category, three simulations
## -----------------------------------------------------------------------------
plot_diff_block <- function(maps.diff.by.sim, category) {
  at.diff <- get_at_diff(category)
  cols.diff <- get_cols_diff(category)
  plots <- list()
  
  for (i in seq_along(simulations)) {
    
    sim <- simulations[i]
    maps.sim <- maps.diff.by.sim[[sim]][methods]
    
    plots[[i]] <- spatialPlot(
      climatology(makeMultiGrid(maps.sim)),
      names.attr = method.names,
      as.table = TRUE,
      layout = c(length(maps.sim), 1),
      backdrop.theme = "countries",
      col.regions = cols.diff,
      set.min = min(at.diff),
      set.max = max(at.diff),
      at = at.diff,
      main = list(paste0(sim, "-REMO  ", toupper(category)), cex = 0.8),
      par.strip.text = list(cex = 1),
      par.settings = list(
        panel.background = list(col = "white"),
        background = list(col = "white"),
        layout.heights = list(
          top.padding = 0,
          bottom.padding = 0,
          main = 1.4,
          main.key.padding = 0.5,
          key.axis.padding = 0,
          axis.xlab.padding = 0,
          xlab = 0
        )
      )
    )
  }
  
  arrangeGrob(
    grobs = plots,
    ncol = 1,
    nrow = length(plots),
    padding = unit(0, "pt"),
    heights = unit(rep(1, length(plots)), "null")
  )
  }
## -----------------------------------------------------------------------------
## One PDF per species
## -----------------------------------------------------------------------------
for (s in seq_along(species.lab)) {
  
  sp_lab <- species.lab[s]
  sp <- species[s]
  categories.obs <- obs.categories.by.species(sp_lab)
  categories.plot <- plot.categories.by.species(sp_lab)
  
  message("\nProcessing PDF for ", sp_lab)
  
  thi.obs <- compute_seasonal_thi(grid.obs, sp, per, months = months)
  maps.obs <- make_obs_maps(thi.obs, sp_lab, categories.obs)
  
  grobs <- list()
  heights <- c()
  
  for (cat in categories.plot) {
    
    saved.obj <- read_diff_file(cat, sp_lab, season.name)
    if (is.null(saved.obj)) next
    
    grobs[[length(grobs) + 1]] <- plot_category_block(
      maps.diff.by.sim = saved.obj$maps.diff.by.sim,
      obs.map = maps.obs[[cat]],
      sp_lab = sp_lab,
      category = cat
    )
    
    heights <- c(heights, 2.25)
  }
  
  if (length(grobs) == 0) {
    warning("No difference grids found for ", sp_lab, ". No PDF will be created.")
    next
  }
  
  png.file <- file.path(
    plot.out.dir,
    paste0("THI_hist-",season.name, "-", sp_lab, ".png")
  )
  
  png(
    filename = png.file,
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
      paste0(toupper(sp_lab),
        "- Heat stress days: observations and methods for JJA (1986-2005)"
      ),
      gp = gpar(fontsize = 18, fontface = "bold")
    )
  )
  
  dev.off()
  
  message("Saved PNG: ", png.file)
  
}
