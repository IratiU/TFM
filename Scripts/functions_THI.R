#' @title Computation of THI on vectors
#' @description Function to compute THI for different livestock species
#' @details 
#' @param Tdb Dry-bulb temperature (in ºC)
#' @param RH Relative humidity (in %)
#' @param species Livestock species of interest: {'cattle', 'ruminants', 'swine', 'poultry'}
#' @author R. Manzanas
#' @export
computeTHI.1D <- function(Tdb, RH, species) {
  if (species == "cattle") {
    # THI for dairy cattle based on National Research Council (1971)
    return(1.8*Tdb + 32 - ((0.55 - 0.0055*RH)*(1.8*Tdb - 26.8)))
  } else if (species == "ruminants") {
    # THI for small ruminants based on LPHSI (1990) modified by Marai et al. (2001)
    return(Tdb - ((0.31 - 0.31*RH/100) * (Tdb - 14.4)))
  } else if (species == "swine") {
    # THI for swine based on Marai et al. (2001)
    Twb = Tdb*atan(0.151977*sqrt(RH + 8.313659)) + atan(Tdb + RH) - atan(RH - 1.676331) + 0.00391838*sqrt(RH^3)*atan(0.023101*RH) - 4.686035  # wet-bulb temperature, based on Wagner and Prub 2002 (http://link.aip.org/link/JPCRBU/v31/i2/p387/s1)
    return(0.75*Tdb + 0.25*Twb)
  } else if (species == "poultry") {
    # THI for poultry based on Zulovich and DeShazer (1990)
    Twb = Tdb*atan(0.151977*sqrt(RH + 8.313659)) + atan(Tdb + RH) - atan(RH - 1.676331) + 0.00391838*sqrt(RH^3)*atan(0.023101*RH) - 4.686035  # wet-bulb temperature, based on Wagner and Prub 2002 (http://link.aip.org/link/JPCRBU/v31/i2/p387/s1)
    return(Tdb*0.6 + Twb*0.4)
  }
}

#' @title Computation of THI on C4R grids
#' @description Function to compute THI for different livestock species
#' @details 
#' @param Tdb Dry-bulb temperature (in ºC)
#' @param RH Relative humidity (in %)
#' @param species Livestock species of interest: {'cattle', 'ruminants', 'swine', 'poultry'}
#' @author R. Manzanas
#' @export
computeTHI <- function(Tdb, RH, species) {
  if (identical(getShape(Tdb), getShape(RH))) {
    ntime = getShape(Tdb)["time"]
    nlat = getShape(Tdb)["lat"]
    nlon = getShape(Tdb)["lon"]
    
    thi.data = array(NA, getShape(Tdb))
    for (ilat in 1:nlat) {
      for (ilon in 1:nlon) {
        thi.data[, ilat, ilon] = computeTHI.1D(Tdb$Data[, ilat, ilon], 
                                               RH$Data[, ilat, ilon],
                                               species = species)
      }
    }
    thi = Tdb; thi$Data = thi.data
    attributes(thi$Variable)$description = "Temperature-Humidity Index"
    attributes(thi$Variable)$longname = "THI"
    attributes(thi$Variable)$units = ""
    attributes(thi$Data)$dimensions = attributes(Tdb$Data)$dimensions
  }
  return(thi)
}

#' @title Analysis of THI on a monthly time-scale
#' @description Function to analyze the results provided by computeTHI on a monthly time-scale
analyzeTHI_M <- function(thi, threshold.min, threshold.max) {
  # thi: Daily values for THI
  # threshold.min: Minimum value of interest
  # threshold.max: Maximum value of interest  
  
  thi.D.bin = gridArithmetics(binaryGrid(thi, condition = "GT", threshold = threshold.min),
                              binaryGrid(thi, condition = "LE", threshold = threshold.max), 
                              operator = "*")
  
  thi.M.ndays = aggregateGrid(thi.D.bin, 
                              aggr.m = list(FUN = "sum", na.rm = FALSE))
  
  ## output
  return(thi.M.ndays)
}

#' @title Analysis of THI on a yearly time-scale
#' @description Function to analyze the results provided by computeTHI on a yearly time-scale
analyzeTHI_Y <- function(thi, threshold.min, threshold.max) {
  # thi: Daily values for THI
  # threshold.min: Minimum value of interest
  # threshold.max: Maximum value of interest  

  thi.D.bin = gridArithmetics(binaryGrid(thi, condition = "GT", threshold = threshold.min),
                              binaryGrid(thi, condition = "LE", threshold = threshold.max), 
                              operator = "*")
  thi.D.bin$Dates$start <- as.POSIXct(thi.D.bin$Dates$start, tz = "GMT")
  thi.D.bin$Dates$end   <- as.POSIXct(thi.D.bin$Dates$end,   tz = "GMT")
  thi.Y.ndays = aggregateGrid(thi.D.bin, 
                              aggr.y = list(FUN = "sum", na.rm = TRUE))
 thi.Y.ndays.trend = climatology(thi.Y.ndays, clim.fun = list(FUN = "computeTrend"))
 thi.Y.ndays.sig = map.stippling(clim = climatology(thi.Y.ndays, clim.fun = list(FUN = "computeSigTrend")), 
                                threshold = 0.05, condition = "LT", 
                                pch = 19, cex = .1, col = "black")  # points exhibiting significant trends
  ## output

  out = list()
  out$ndays = thi.Y.ndays
  out$ndays.trend = thi.Y.ndays.trend
  out$ndays.sig = thi.Y.ndays.sig
  return(out)
}
