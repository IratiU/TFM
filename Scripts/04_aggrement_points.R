library(visualizeR)
library(transformeR)

# ------------------------------------------------------------------
# Funcion para crear un grid donde las 3 simulaciones tienen
# el mismo signo en las deltas
# ------------------------------------------------------------------

make_sign_agreement_grid <- function(delta_MPI, delta_HadGEM, delta_NorESM, eps = 0) {
  
  stopifnot(identical(dim(delta_MPI$Data), dim(delta_HadGEM$Data)))
  stopifnot(identical(dim(delta_MPI$Data), dim(delta_NorESM$Data)))
  
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
  
  return(agreement_grid)
}


# ------------------------------------------------------------------
# Funcion para extraer el grid de una especie, metodo y categoria
# ------------------------------------------------------------------

get_delta_grid <- function(obj, method, sp, cat) {
  
  grid <- obj$list.maps.ndays.Y[[sp]][[method]][[cat]]
  
  if (is.null(grid)) {
    stop("No existe el grid para: ", sp, " - ", method, " - ", cat)
  }
  
  if (is.null(grid$Data) || length(grid$Data) == 0) {
    stop("Grid vacio para: ", sp, " - ", method, " - ", cat)
  }
  
  attr(grid$Data, "climatology:fun") <- "none"
  
  return(grid)
}


# ------------------------------------------------------------------
# Leer los RDS
# ------------------------------------------------------------------

models <- c("MPI", "HadGEM", "NorESM")
GWL <- c("GWL1.5", "GWL2", "GWL3")

grids <- list()

for (gwl in GWL) {
  
  grids[[gwl]] <- list()
  
  for (m in models) {
    
    file <- file.path(
      "~/lustre/gmeteo/WORK/uribei/TFM/Data",
      m,
      gwl,
      sprintf("THI_deltas_results_%s_%s.rds", gwl, season_to_plot)
    )
    
    grids[[gwl]][[m]] <- readRDS(file)
  }
}


# ------------------------------------------------------------------
# Crear grids y puntos de acuerdo de signo
# ------------------------------------------------------------------

species <- c("dairy", "beef", "sheep", "goats", "poultry", "swine")
methods <- c("rcm", "eqm", "qdm", "mbcn")
categories <- c("mild", "moderate", "severe")

agreement_grids <- list()
agreement_pts <- list()

for (gwl in GWL) {
  
  agreement_grids[[gwl]] <- list()
  agreement_pts[[gwl]] <- list()
  
  for (sp in species) {
    
    agreement_grids[[gwl]][[sp]] <- list()
    agreement_pts[[gwl]][[sp]] <- list()
    
    for (method in methods) {
      
      agreement_grids[[gwl]][[sp]][[method]] <- list()
      agreement_pts[[gwl]][[sp]][[method]] <- list()
      
      for (cat in categories) {
        
        delta_MPI <- get_delta_grid(
          grids[[gwl]][["MPI"]],
          method,
          sp,
          cat
        )
        
        delta_HadGEM <- get_delta_grid(
          grids[[gwl]][["HadGEM"]],
          method,
          sp,
          cat
        )
        
        delta_NorESM <- get_delta_grid(
          grids[[gwl]][["NorESM"]],
          method,
          sp,
          cat
        )
        
        agreement_grids[[gwl]][[sp]][[method]][[cat]] <-
          make_sign_agreement_grid(
            delta_MPI,
            delta_HadGEM,
            delta_NorESM
          )
        
        agreement_pts[[gwl]][[sp]][[method]][[cat]] <-
          map.stippling(
            clim      = agreement_grids[[gwl]][[sp]][[method]][[cat]],
            threshold = 0.5,
            condition = "GT",
            pch = 19,
            col = "black",
            cex = 0.12
          )
      }
    }
  }
}

