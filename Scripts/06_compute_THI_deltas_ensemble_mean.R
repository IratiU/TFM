# =============================================================================
# 06_compute_THI_deltas_ensemble_mean.R
#
# A partir de las deltas THI corregidas de los 3 GCMs (MPI, HadGEM, NorESM):
#   1) Calcula la media del ensemble (ensemble mean) de deltas.
#   2) Calcula, para cada especie/metodo/categoria, los puntos donde los 3
#      GCMs coinciden en el signo de la delta (acuerdo / stippling).
#
# Guarda un unico RDS por GWL con ambos resultados, que consume
# 04_plot_thi_deltas_mean.R
# =============================================================================

# ---- Librerias ----------------------------------------------------------------
library(transformeR)
library(visualizeR)

# ---- Parametros -----------------------------------------------------------------
models <- c("MPI", "HadGEM", "NorESM")
gwls   <- c("GWL1.5", "GWL2", "GWL3")

season_to_analyse <- "JJA"      # "ANN", "DJF", "MAM", "JJA", "SON"

species    <- c("dairy", "beef", "sheep", "goats", "poultry", "swine")
methods    <- c("rcm", "eqm", "qdm", "mbcn")
categories <- c("mild", "moderate", "severe")

dir_project <- "/lustre/gmeteo/WORK/uribei/TFM"
dir_base    <- file.path(dir_project, "Data")

delta_file_model <- function(model, gwl) {
  file.path(
    dir_base, model, gwl,
    paste0("THI_deltas_results_", gwl, "_", season_to_analyse, ".rds")
  )
}

out_file_ensemble <- function(gwl) {
  out_dir <- file.path(dir_base, "ENSEMBLE_MEAN", gwl)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  file.path(out_dir, sprintf("THI_ENSEMBLE_MEAN_%s_%s.rds", gwl, season_to_analyse))
}

# ---- Funciones: media de grids / listas de grids ---------------------------------

mean_grid <- function(g1, g2, g3) {
  g <- gridArithmetics(g1, g2, operator = "+")
  g <- gridArithmetics(g,  g3, operator = "+")
  gridArithmetics(g, 3, operator = "/")
}

mean_list <- function(x1, x2, x3) {

  if ("Data" %in% names(x1)) {
    return(mean_grid(x1, x2, x3))
  }

  out <- list()
  for (nm in names(x1)) {
    out[[nm]] <- mean_list(x1[[nm]], x2[[nm]], x3[[nm]])
  }
  out
}

# ---- Funciones: acuerdo de signo entre los 3 GCMs --------------------------------

get_delta_grid <- function(obj, sp, method, category) {

  grid <- obj$list.maps.ndays.Y[[sp]][[method]][[category]]

  if (is.null(grid) || is.null(grid$Data) || length(grid$Data) == 0) {
    stop("No existe (o esta vacio) el grid para: ", sp, " - ", method, " - ", category)
  }

  attr(grid$Data, "climatology:fun") <- "none"
  grid
}

make_sign_agreement_grid <- function(delta_MPI, delta_HadGEM, delta_NorESM, eps = 0) {

  stopifnot(
    identical(dim(delta_MPI$Data), dim(delta_HadGEM$Data)),
    identical(dim(delta_MPI$Data), dim(delta_NorESM$Data))
  )

  same_positive <- delta_MPI$Data >  eps &
    delta_HadGEM$Data >  eps &
    delta_NorESM$Data >  eps

  same_negative <- delta_MPI$Data < -eps &
    delta_HadGEM$Data < -eps &
    delta_NorESM$Data < -eps

  same_sign <- same_positive | same_negative
  same_sign[is.na(same_sign)] <- FALSE

  agreement_grid <- delta_MPI
  agreement_grid$Data[] <- NA_real_
  agreement_grid$Data[same_sign] <- 1

  attr(agreement_grid$Data, "climatology:fun") <- "none"
  agreement_grid
}

# ---- Bucle por GWL ----------------------------------------------------------------
for (gwl in gwls) {

  message("Calculando ensemble mean y puntos de acuerdo para ", gwl)

  # Cargar los 3 modelos una unica vez: se reutilizan tanto para el ensemble
  # mean como para los puntos de acuerdo de signo (antes se leian dos veces)
  r <- list()
  for (m in models) {
    f <- delta_file_model(m, gwl)
    stopifnot(file.exists(f))
    r[[m]] <- readRDS(f)
  }

  # -- Ensemble mean ----------------------------------------------------------------
  ensemble <- list()

  ensemble$maps.dc                  <- mean_list(r$MPI$maps.dc, r$HadGEM$maps.dc, r$NorESM$maps.dc)
  ensemble$list.maps.ndays.Y        <- mean_list(r$MPI$list.maps.ndays.Y, r$HadGEM$list.maps.ndays.Y, r$NorESM$list.maps.ndays.Y)
  ensemble$list.maps.ndays.diff     <- mean_list(r$MPI$list.maps.ndays.diff, r$HadGEM$list.maps.ndays.diff, r$NorESM$list.maps.ndays.diff)
  ensemble$list.maps.ndays.diff.eqm <- mean_list(r$MPI$list.maps.ndays.diff.eqm, r$HadGEM$list.maps.ndays.diff.eqm, r$NorESM$list.maps.ndays.diff.eqm)

  ensemble$tit.days          <- r$MPI$tit.days
  ensemble$tit.days.diff     <- r$MPI$tit.days.diff
  ensemble$tit.days.diff.eqm <- r$MPI$tit.days.diff.eqm

  # -- Puntos de acuerdo de signo -----------------------------------------------------
  agreement_grids <- list()
  agreement_pts    <- list()

  for (sp in species) {

    agreement_grids[[sp]] <- list()
    agreement_pts[[sp]]    <- list()

    for (method in methods) {

      agreement_grids[[sp]][[method]] <- list()
      agreement_pts[[sp]][[method]]    <- list()

      for (cat in categories) {

        delta_MPI    <- get_delta_grid(r$MPI,    sp, method, cat)
        delta_HadGEM <- get_delta_grid(r$HadGEM, sp, method, cat)
        delta_NorESM <- get_delta_grid(r$NorESM, sp, method, cat)

        agreement_grids[[sp]][[method]][[cat]] <- make_sign_agreement_grid(
          delta_MPI, delta_HadGEM, delta_NorESM
        )

        agreement_pts[[sp]][[method]][[cat]] <- map.stippling(
          clim      = agreement_grids[[sp]][[method]][[cat]],
          threshold = 0.5,
          condition = "GT",
          pch = 19,
          col = "black",
          cex = 0.12
        )
      }
    }
  }

  ensemble$agreement_grids <- agreement_grids
  ensemble$agreement_pts   <- agreement_pts

  rm(r)
  gc()

  # -- Guardado -----------------------------------------------------------------------
  out_file <- out_file_ensemble(gwl)
  saveRDS(ensemble, out_file)
  message("Guardado: ", out_file)
}