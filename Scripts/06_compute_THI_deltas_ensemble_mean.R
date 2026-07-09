# =============================================================================
# Media de deltas THI entre MPI, HadGEM y NorESM
# =============================================================================

library(transformeR)

models <- c("MPI", "HadGEM", "NorESM")
gwls   <- c("GWL1.5","GWL2", "GWL3")
# season <- "JJA"

dir_base <- "/lustre/gmeteo/WORK/uribei/TFM/Data"

# -----------------------------------------------------------------------------
# Media de 3 grids climate4R
# --------------------------.----------------------------------------------------

mean_grid <- function(g1, g2, g3) {
  g <- gridArithmetics(g1, g2, operator = "+")
  g <- gridArithmetics(g,  g3, operator = "+")
  g <- gridArithmetics(g,  3,  operator = "/")
  return(g)
}

# -----------------------------------------------------------------------------
# Media de listas que contienen grids
# -----------------------------------------------------------------------------

mean_list <- function(x1, x2, x3) {
  
  if ("Data" %in% names(x1)) {
    return(mean_grid(x1, x2, x3))
  }
  
  out <- list()
  
  for (nm in names(x1)) {
    out[[nm]] <- mean_list(
      x1[[nm]],
      x2[[nm]],
      x3[[nm]]
    )
  }
  
  return(out)
}

# -----------------------------------------------------------------------------
# Bucle por GWL
# -----------------------------------------------------------------------------

for (gwl in gwls) {
  
  message("Calculando media para ", gwl)
  
  f_MPI <- file.path(
    dir_base, "MPI", gwl,
    paste0("THI_deltas_results_", gwl, "_", season, ".rds")
  )
  
  f_HadGEM <- file.path(
    dir_base, "HadGEM", gwl,
    paste0("THI_deltas_results_", gwl, "_", season, ".rds")
  )
  
  f_NorESM <- file.path(
    dir_base, "NorESM", gwl,
    paste0("THI_deltas_results_", gwl, "_", season, ".rds")
  )
  
  r_MPI    <- readRDS(f_MPI)
  r_HadGEM <- readRDS(f_HadGEM)
  r_NorESM <- readRDS(f_NorESM)
  
  ensemble <- list()
  
  ensemble$maps.dc <- mean_list(
    r_MPI$maps.dc,
    r_HadGEM$maps.dc,
    r_NorESM$maps.dc
  )
  
  ensemble$list.maps.ndays.Y <- mean_list(
    r_MPI$list.maps.ndays.Y,
    r_HadGEM$list.maps.ndays.Y,
    r_NorESM$list.maps.ndays.Y
  )
  
  ensemble$list.maps.ndays.diff <- mean_list(
    r_MPI$list.maps.ndays.diff,
    r_HadGEM$list.maps.ndays.diff,
    r_NorESM$list.maps.ndays.diff
  )
  
  ensemble$list.maps.ndays.diff.eqm <- mean_list(
    r_MPI$list.maps.ndays.diff.eqm,
    r_HadGEM$list.maps.ndays.diff.eqm,
    r_NorESM$list.maps.ndays.diff.eqm
  )
  
  ensemble$tit.days          <- r_MPI$tit.days
  ensemble$tit.days.diff     <- r_MPI$tit.days.diff
  ensemble$tit.days.diff.eqm <- r_MPI$tit.days.diff.eqm
  
  out_dir <- file.path(dir_base, "ENSEMBLE_MEAN", gwl)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  out_file <- file.path(
    out_dir,
    paste0("THI_ENSEMBLE_MEAN_", gwl, "_", season, ".rds")
  )
  
  saveRDS(ensemble, out_file)
  
  message("Guardado: ", out_file)
}
  
  