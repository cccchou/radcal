# ResAge R code for reservoir age offset calculation. Script colyear.R
# Guillaume Soulet <gsoulet@whoi.edu>
# see accompanying manual and Soulet 2015 (Quaternary Geochronology 29:97-103, doi: 10.1016/j.quageo.2015.05.023)
# updated with reservoir 14C disequilibria F14R (f.14R), delta14R (delta.14R) and dcp (dcp) 

colyear <- function(name="example_colyear", cc=1, cc1="IntCal13.14C", cc2="SHCal13.14C", cc3="IntCAL09.14C", ext=".csv", sep=",", yrsteps=1, delta.14R=FALSE, f.14R=FALSE, dcp=FALSE)
  .colyear(name, cc, cc1, cc2, cc3, ext, sep, yrsteps, delta.14R, f.14R, dcp)

.colyear <- function(name, cc, cc1, cc2, cc3, ext, sep, yrsteps, delta.14R, f.14R, dcp)
{
  dets <- suppressWarnings(read.csv(paste("InputData/", name, "/", name, ext, sep=""), sep=sep))
  
  # check for consistency 
  if(ncol(dets) != 4)
    stop("Your input file should contain 4 columns. Please, check the template designed for the colyear function.", call.=FALSE)
  
  if(min(dets[2]) < 0)
    stop("At least one of your radiocarbon data is negative.\n Colyear is not designed to work with post-bomb samples. Remove such data from the input file.", call.=FALSE)
  
  if(max(dets[4]) > 1950 || min(dets[4]) < 1500)
    stop("Colyear function is only designed for sample with calendar age between 1950 and 1500 AD. Remove such data from the intput file", call.=FALSE)
  
  if(any(dets[4] <= 450))
    stop("Did you already convert sample collection year from AD/BC timescale to the BP timescale?\n Please let collection year data in the AD/BC timescale.", call.=FALSE)
    
  if(yrsteps != 1)
    stop("The collection year (colyear) function requires the yrsteps argument is set to 1.\n In calling function, set yrsteps=1.", call.=FALSE)
 
  # set the calibration curve
  if(cc==1) calcurve <- read.table(cc1) else
    if(cc==2) calcurve <- read.table(cc2) else
      if(cc==3) calcurve <- read.table(cc3) else
      stop("I do not understand which calibration curve you mean.\n Available calibration curves are IntCal13 (cc=1), SHCal13 (cc=2), and IntCal09 (cc=3).\n Please check the manual.", call.=FALSE)
  
  if(cc==1) ccname <- "IntCal13" else
    if(cc==2) ccname <- "SHCal13" else
      if(cc==3) ccname <- "IntCal09" else
        stop()
  
  if(cc==1) cc <- "_using_IntCal13" else
    if(cc==2) cc <- "_using_SHCal13" else
      if(cc==3) cc <- "_using_IntCal09" else
      stop()
  
  cat("Reservoir age offsets are calculated using the", ccname, "calibration curve.\n")
  
  # interpolate calibration curve between 0 and 500 BP, i.e. between 1950 and 1500 AD
  cc.x <- seq(0, 450, by=yrsteps)
  cc.y <- approx(calcurve[,1], calcurve[,2], cc.x)$y
  cc.sd <- approx(calcurve[,1], calcurve[,3], cc.x)$y
  calcurve.interp <- cbind(cc.x, cc.y, cc.sd)
  
  # now work in the BP timescale
  colyr_AD <- dets$Collection_year_AD
  Collection_year_BP <- 1950 - colyr_AD
  
  # prepare for reservoir age calculation 
  rad_res <- dets$Radiocarbon_age
  rad_res_sd <- dets$Radiocarbon_age_error
  
  # calculate sample reservoir ages
  Atmospheric_14C_age <- round(calcurve.interp[Collection_year_BP + 1,2]) # find the 14C age of the atmosphere during collection year... 
  Atmospheric_14C_error <- round(calcurve.interp[Collection_year_BP + 1,3]) # ... and associated uncertainty
  Reservoir_age <- round(rad_res - Atmospheric_14C_age) # calculate reservoir age
  Reservoir_age_error <- round(sqrt(rad_res_sd*rad_res_sd + Atmospheric_14C_error*Atmospheric_14C_error)) # propagate errors through quadratic sum
  
  fm.res <- exp(rad_res/-8033)
  fm.res.sd <- exp(-(rad_res-rad_res_sd)/8033) - fm.res
  fm.atm <- exp(Atmospheric_14C_age/-8033)
  fm.atm.sd <- exp(-(Atmospheric_14C_age-Atmospheric_14C_error)/8033) - fm.atm
  
  F14R  <- round(fm.res/fm.atm,4) # calculate F14R rounded to the fourth digit
  F14R_error <- round(F14R*sqrt(fm.res.sd*fm.res.sd/(fm.res*fm.res)+fm.atm.sd*fm.atm.sd/(fm.atm*fm.atm)),4) # propagates errors
  delta14R <- 1000*(F14R-1) # calculate delta14R
  delta14R_error <- 1000*F14R_error # propagate errors
  DCP <- 100*(1-F14R) # calculate dcp
  DCP_error <- 100*F14R_error
  
  output <- cbind(dets, Collection_year_BP, Atmospheric_14C_age, Atmospheric_14C_error, Reservoir_age, Reservoir_age_error) 
  
  if(f.14R) output <- cbind(output, F14R, F14R_error)
  if(delta.14R) output <- cbind(output, delta14R, delta14R_error)
  if(dcp) output <- cbind(output, DCP, DCP_error)
  
  #View(output)
  
  # export data as .csv file
  write.csv(output,paste("InputData/", name, "/", name,cc,"_out", ext, sep=""), row.names = FALSE)
  
  # Done message
  cat("Work has been done successfully, check output in your folder.\n\n")
  
}

# list the available data
Data <- list.files("InputData/")

# Welcome
cat("Welcome, ResAge is here to help you with reservoir age calculation.\n")
cat("The colyear function is designed to work with pairs of reservoir-derived 14C age and perfectly known calendar age.\n\n")