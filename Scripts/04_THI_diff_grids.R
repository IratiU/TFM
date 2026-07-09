library(convertR)
library(gridExtra)
source("/lustre/gmeteo/WORK/uribei/TFM/Scripts/THI/functions_THI.R")
library(transformeR)
library(loadeR)
library(visualizeR)
library(sp)
library(RColorBrewer)
library(lattice)

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

base.dir <- "/lustre/gmeteo/WORK/uribei/TFM/Data"
diro <- "/lustre/gmeteo/WORK/uribei/TFM/Data"
grid.out.dir <- file.path(diro, "THI_saved_grids")
dir.create(grid.out.dir, recursive = TRUE, showWarnings = FALSE)

per <- 1986:2005
simulations <- c("NorESM", "MPI", "HadGEM")
methods <- c("eqm", "qdm", "mbcr", "mbcp", "mbcn")
species <- c("cattle", "cattle", "ruminants", "ruminants", "poultry", "swine")
species.lab <- c("dairy", "beef", "sheep", "goats", "poultry", "swine")

mild <- list(c(72, 78), c(75, 81), c(25, 30), c(27, 32), c(27.8, 28.8), c(23.33, 26.11))
moderate <- list(c(78, 88), c(81, 91), c(30, 35), c(32, 37), c(28.8, 30), c(26.11, 28.88))
severe <- list(c(88, Inf), c(91, Inf), c(35, Inf), c(37, Inf), c(30, Inf), c(28.88, Inf))

df <- data.frame(rbind(mild, moderate, severe))
colnames(df) <- species.lab
category <- "mild"

# season.def <- list(DJF = c(12, 1, 2), MAM = c(3, 4, 5), JJA = c(6, 7, 8), SON = c(9, 10, 11))
season.def <- list(JJA = c(6,7,8))

load("/lustre/gmeteo/WORK/uribei/TFM/Data/tas_obs.rda")
load("/lustre/gmeteo/WORK/uribei/TFM/Data/hurs_obs.rda")
tas.era5$Variable$varName <- "tas"
grid_obs <- makeMultiGrid("tas" = tas.era5, "hurs" = hurs.era5)
grid_obs <- subsetGrid(grid_obs, years = per)

obs.mask <- subsetGrid(grid_obs, var = "tas")
obs.mask <- climatology(obs.mask)
mask.na <- is.na(drop(obs.mask$Data))

grids <- list()
for (sim in simulations) {
  grids[[sim]] <- list()
  for (m in methods) {
    file_path <- file.path(base.dir, sim, "grids", paste0(m, "_hist.rds"))
    grids[[sim]][[m]] <- readRDS(file_path)
  }
}

compute_seasonal_thi <- function(grid, sp, per, months) {
  
  thi <- vector("list", length(months))
  
  for (i in seq_along(months)) {
    
    M <- months[i]
    
    data <- subsetGrid(grid, season = M, years = per, var = "tas")
    attributes(data$Variable)$units <- "degC"
    tasmax <- data
    rm(data)
    
    data <- subsetGrid(grid, season = M, years = per, var = "hurs")
    hurs <- data
    rm(data)
    
    thi[[i]] <- computeTHI(tasmax, hurs, species = sp)
  }
  
  bindGrid(thi, dimension = "time")
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
  return(land.mask)
}
for (season.name in names(season.def)) {
  
  months <- season.def[[season.name]]
  
  message("\n==============================")
  message("Processing season: ", season.name)
  message("==============================")
  
  for (s in seq_along(species.lab)) {
    
    sp_lab <- species.lab[s]
    sp <- species[s]
    
    message("\nProcessing ", sp_lab, " - ", toupper(category), " - ", season.name)
    
    thi.obs <- compute_seasonal_thi(
      grid = grid_obs,
      sp = sp,
      per = per,
      months = months
    )
    
    thi.Y.obs <- analyzeTHI_Y(
      thi = thi.obs,
      threshold.min = df[[sp_lab]][[category]][1],
      threshold.max = df[[sp_lab]][[category]][2]
    )
    
    land.mask <- make_land_mask(thi.Y.obs$ndays, mask.na)
    
    maps.diff.by.sim <- setNames(lapply(simulations, function(x) list()), simulations)
    
    for (sim in simulations) {
      
      message("  Simulation: ", sim)
      
      for (m in methods) {
        
        message("    Method: ", m)
        
        thi.D <- compute_seasonal_thi(
          grid = grids[[sim]][[m]],
          sp = sp,
          per = per,
          months = months
        )
        
        thi.Y <- analyzeTHI_Y(
          thi = thi.D,
          threshold.min = df[[sp_lab]][[category]][1],
          threshold.max = df[[sp_lab]][[category]][2]
        )
        
        ndays.diff <- gridArithmetics(thi.Y$ndays, thi.Y.obs$ndays, operator = "-")
        
        ndays.diff <- gridArithmetics(ndays.diff, land.mask, operator = "*")
        
        maps.diff.by.sim[[sim]][[m]] <- ndays.diff
        
        message("      NAs after mask: ", sum(is.na(ndays.diff$Data)))
      }
    }
    
    out.file <- file.path(
      grid.out.dir,
      paste0(
        "THI_diff_grids_",
        toupper(category), "_",
        season.name, "_",
        sp_lab,
        ".rds"
      )
    )
    
    saveRDS(
      object = list(
        maps.diff.by.sim = maps.diff.by.sim,
        category = category,
        season = season.name,
        months = months,
        sp_lab = sp_lab,
        sp = sp,
        simulations = simulations,
        methods = methods,
        period = per
      ),
      file = out.file
    )
    
    message("Saved: ", out.file)
  }
}
