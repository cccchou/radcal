# ResAge R code for reservoir age offset calculation. Script rad2.R
# Guillaume Soulet <gsoulet@whoi.edu>
# see accompanying manual and Soulet 2015 (Quaternary Geochronology 29:97-103, doi: 10.1016/j.quageo.2015.05.023)
# updated with reservoir 14C disequilibria F14R (f.14R), delta14R (delta.14R) and dcp (dcp) 

rad2 <- function(name="example_rad2", cc=1, cc1="IntCal20.14C", cc2="SHCal13.14C", cc3="IntCal09.14C", calibrate=FALSE, ext=".csv", sep=",", yrsteps=1, prob=0.95, threshold=1e-6, hpdsteps=1, storedat=FALSE, export.cal.pdf=FALSE, delta.14R=FALSE, f.14R=FALSE, dcp=FALSE)
  .rad2(name, cc, cc1, cc2, cc3, calibrate, ext, sep, yrsteps, prob, threshold, hpdsteps, storedat, export.cal.pdf, delta.14R, f.14R, dcp)

.rad2 <- function(name, cc, cc1, cc2, cc3, calibrate, ext, sep, yrsteps, prob, threshold, hpdsteps, storedat, export.cal.pdf, delta.14R, f.14R, dcp)
{
  dets <- suppressWarnings(read.csv(paste("InputData/", name, "/", name, ext, sep=""), sep=sep))
  

  # check for consistency 
  if(ncol(dets) != 5)
    stop("Your input file should contain 5 columns.\n Please, check the template designed for the radrad function.", call.=FALSE)
  
  if(min(dets[2]) < 0 || min(dets[4]) < 0)
    stop("At least one of your radiocarbon data is negative. ResAge is not designed to work with post-bomb samples.\n Remove such data from the input file.", call.=FALSE)
  
  # prepare for reservoir age calculation 
  id <- dets$id
  rad_res <- dets$Reservoir_derived_14C_age
  rad_res_sd <- dets$Reservoir_derived_14C_error
  rad_atm <- dets$Atmosphere_derived_14C_age
  rad_atm_sd <- dets$Atmosphere_derived_14C_error
  fm.res <- exp(rad_res/-8033)
  fm.res.sd <- exp(-(rad_res-rad_res_sd)/8033) - fm.res
  fm.atm <- exp(rad_atm/-8033)
  fm.atm.sd <- exp(-(rad_atm-rad_atm_sd)/8033) - fm.atm
  
  
  # calculate sample reservoir ages
  
  cat("Calculating reservoir age offsets...\n")
  
  Reservoir_age_offset <- round(rad_res - rad_atm) # calculate reservoir age
  Reservoir_age_offset_error <- round(sqrt(rad_res_sd*rad_res_sd + rad_atm_sd*rad_atm_sd)) # propagate errors through quadratic sum
  
  F14R  <- round(fm.res/fm.atm,4) # calculate F14R rounded to the fourth digit
  F14R_error <- round(F14R*sqrt(fm.res.sd*fm.res.sd/(fm.res*fm.res)+fm.atm.sd*fm.atm.sd/(fm.atm*fm.atm)),4) # propagates errors
  delta14R <- 1000*(F14R-1) # calculate delta14R
  delta14R_error <- 1000*F14R_error # propagate errors
  DCP <- 100*(1-F14R) # calculate dcp
  DCP_error <- 100*F14R_error
  
  output <- cbind(dets, Reservoir_age_offset, Reservoir_age_offset_error) 
  
  if(f.14R) output <- cbind(output, F14R, F14R_error)
  if(delta.14R) output <- cbind(output, delta14R, delta14R_error)
  if(dcp) output <- cbind(output, DCP, DCP_error)
  
# find the calibrated distributions of 14C dates
.caldist <- function(f.cage, f.error, theta, f.mu, f.sigma, yrsteps, threshold)
  {
    # calibrate; find how far f.cage (measurement) is from f.mu (calibration curve)
    cal <- cbind(theta, dnorm(f.mu, f.cage, sqrt(f.error^2+f.sigma^2)))
  
    # interpolate and normalise calibrated distribution to 1
    cal <- cal[min(which(cal[,2] > 0)):max(which(cal[,2] > 0)),] # remove unnecessary data
    cal <- approx(cal[,1], cal[,2], seq(min(cal[,1]), max(cal[,1]), by=yrsteps))
    cal <- cbind(cal$x, cal$y/sum(cal$y))
    
    # only report those normalised calibrated probabilities beyond a threshold
    cal[cal[,2] > threshold,]
  }
  
  
# find the highest posterior density (hpd) of the calibrated distribution
.hpd <- function(dat, prob, hpdsteps, yrsteps)
  {
    # interpolate and rank the ages according to their calibrated distribution probabilities
    dat <- approx(dat[,1], dat[,2], seq(min(dat[,1]), max(dat[,1]), by=yrsteps))
    o <- order(dat$y, decreasing=TRUE)
    dat <- cbind(dat$x[o], dat$y[o]/sum(dat$y))
    
    # only retain those ages with cumulative normalised probabilities within required percentage
    dat <- dat[which(cumsum(dat[,2]) <= prob),]
    dat <- dat[order(dat[,1]),]
    
    # identify any individual ranges within the hpd range and calculate their probability
    dif <- which(diff(dat[,1]) > hpdsteps)
    
    if(length(dif)==0)
      hpds <- cbind(min(dat[,1]), max(dat[,1]), 100*prob) else
      {
        dif <- c(dat[1,1], sort(c(dat[dif,1], dat[dif+1,1])), dat[nrow(dat),1])
        dif <- matrix(dif, ncol=2, byrow=TRUE)
        probs <- c()
        for(i in 1:nrow(dif))
          probs[i] <- round(100*sum(dat[which(dat[,1]==dif[i,1]):which(dat[,1]==dif[i,2]),2]), 1)
        hpds <- cbind(dif, probs)
      }
    hpds
  }
  
  # if calibrate=TRUE, let's calibrate your Atmosphere_derived_14C_ages
  if(calibrate)
    { 
          
      # set the calibration curve
      if(cc==1) calcurve <- read.table(cc1) else
        if(cc==2) calcurve <- read.table(cc2) else
          if(cc==3) calcurve <- read.table(cc3) else
            stop("I do not understand which calibration curve you mean.\n Available calibration curves are IntCal13 (cc=1), SHCal13 (cc=2), and IntCal09 (cc=3).\n Please check the manual.", call.=FALSE)
      
      if(cc==1) ccname <- "IntCal20" else
        if(cc==2) ccname <- "SHCal13" else
          if(cc==3) ccname <- "IntCal09" else
            stop()
      
      if(cc==1) cc <- "_using_IntCal20" else
        if(cc==2) cc <- "_using_SHCal13" else
          if(cc==3) cc <- "_using_IntCal09" else
            stop()
      
      cat("Calibrating atmospheric 14C ages using the", ccname, "calibration curve...\n")
      
      # check whether date lies partly or entirely beyond the calibration curve
      border <- 0
      if(any(rad_atm - rad_atm_sd < min(calcurve[,2] + calcurve[,3])))
        if(any(rad_atm + rad_atm_sd > min(calcurve[,2] - calcurve[,3])))
          border <- 1 else border <- 2
      if(any(rad_atm + rad_atm_sd > max(calcurve[,2] - calcurve[,3])))
        if(any(rad_atm - rad_atm_sd < max(calcurve[,2] + calcurve[,3])))
          border <- 1 else border <- 2
      if(border==1)
        cat("\nAt least one 14C age falls partly beyond calibration curve and will be truncated")
      if(border==2)
        stop("\nAt least one 14C age cannot be calibrated because falling beyond calibration curve! Please remove those\n\n", call.=FALSE)
      
      # prepare calcurve for calculation in F14C
      theta <- calcurve[,1]
      f.mu <- exp(-calcurve[,2]/8033)
      f.sigma <- exp(-(calcurve[,2]-calcurve[,3])/8033) - f.mu
  
      # prepare data for calibration
      f.cage <- exp(-rad_atm/8033)
      f.error <- exp(-(rad_atm-rad_atm_sd)/8033) - f.cage
      
      # calibrate the date and report its highest posterior density (hpd) range
      
      dat <- data.frame(cbind(f.cage, f.error), row.names = NULL)
      
      for(i in 1:length(f.cage))
        { 
          calib <- .caldist(dat$f.cage[[i]], dat$f.error[[i]], theta, f.mu, f.sigma, yrsteps, threshold)
          dat$calib[[i]] <- calib # resulting calibrated distribution
          dat$hpd[[i]] <- .hpd(calib, prob=prob, hpdsteps=hpdsteps, yrsteps=yrsteps) # highest probability density corresponding to prob
          dat$mid1[[i]] <- round((dat$hpd[[i]][1] + dat$hpd[[i]][2*nrow(dat$hpd[[i]])])/2) # midpoints of calibrated ranges
          yrs <- calib[,1]
          dat$mid2[[i]] <- round(mean(c(max(yrs), min(yrs)))) # midpoints of entire calibrated distributions (with probabilities beyond threshold) 
          dat$wmn[[i]] <- round(weighted.mean(calib[,1], 1/calib[,2])) # weighted means of calibrated ranges
          dat$med[[i]] <- calib[max(which(cumsum(calib[,2]) <= .5)),1] # medians of calibrated distributions
          dat$mode[[i]] <- calib[which(calib[,2] == max(calib[,2])),1][1] #maximum densities of calibrated distributions
        }
      
      if(storedat) dat <<- dat
      
     # export the calibrated probability density if you required it
     if(export.cal.pdf)
     {
       cat("Exporting your calibrated ages...\n")
       
       # find the min and max cal ages among all the cal age density probability
       mini <- c()
       maxi <- c()
       
       for(i in 1:nrow(dat))
       {
         mini[[i]] <- min(as.data.frame(dat$calib[[i]])[,1])
         maxi[[i]] <- max(as.data.frame(dat$calib[[i]])[,1])
       }
       min <- min(unlist(mini)) 
       max <- max(unlist(maxi))
       
       # create a matrix of zeros and paste inside all the cal age densities 
       cal.age <- matrix(0, ncol=1+nrow(dat), nrow=max-min+1)
       cal.age[1:nrow(cal.age), 1] <- seq(min,max, by=yrsteps)
       
       header <- c()
       
       for(i in 1:nrow(dat))
       {
         d <- max(diff(as.data.frame(dat$calib[i])[,1]))
         if (d > 1)
           stop("Please, decrease threshold parameter: e.g. threshold = 1e-14 in the comand line\n\n", call.=FALSE)
         
         cal.age[(min(as.data.frame(dat$calib[[i]])[,1])-min+1):(max(as.data.frame(dat$calib[[i]])[,1])-min+1), i+1] <- as.data.frame(dat$calib[[i]])[,2]
         header[[i]] <- matrix(id[i])
       }
       header <- header
       colnames(cal.age) <- c("cal_year_BP", header)
       
       write.csv(cal.age,paste("InputData/", name, "/", name,cc,"_","Calibrated_ages_pdfs", ext, sep=""), row.names = FALSE)
      
     }
     
    # export calibrated ranges and estimators
        
    for(i in 1:nrow(dat))
    {
      hpds <- dat$hpd[[i]]
      if(nrow(hpds) == 1) test <- 1 else test <- 2
      if(test == 2) break
    }
    
    if(test==1) 
    {
      export.calib <- matrix(0, ncol=5, nrow=nrow(dat))
      for(i in 1:nrow(dat)) 
      {
        hpds <- dat$hpd[[i]]
        export.calib[i,1] <- hpds[,1]
        export.calib[i,2] <- hpds[,2]
        export.calib[i,3] <- dat$mode[[i]]
        export.calib[i,4] <- dat$mid1[[i]] 
        export.calib[i,5] <- dat$med[[i]]
        
      }
      
      colnames(export.calib) <- c(paste("CalibratedMinRange","_at_",100*prob,"%", sep=""),paste("CalibratedMaxRange","_at_",100*prob,"%", sep=""), "Mode", "MidRange", "Median") 
      cal.output <- cbind(as.data.frame(id), export.calib)
      write.csv(cal.output,paste("InputData/", name, "/", name,cc,"_","Calibrated_age_ranges_at", 100*prob, "%", ext, sep=""), row.names = FALSE)
      
    } else
    {  
     ref <- "Calibrated_"
    
     hpd.file <- file(paste("InputData/", name, "/", name,cc,"_",ref,"ranges.txt", sep=""), "w")
     cat(paste(ref,"ranges at ", 100*prob, "% confidence intervals\n", sep=""), file=hpd.file)
     
     for(i in 1:nrow(dat))
     {
       cat(paste("\n\nIdentification_number: ", dets$id[[i]], "\nCal_year_MinRange\tCal_year_MaxRange\tprobability\n"), file=hpd.file)
       hpds <- dat$hpd[[i]]
       pdf.info <- cbind(dat$mode[[i]], dat$mid1[[i]], dat$med[[i]])
       
       for(j in 1:nrow(hpds))
       {
         for(k in 1:3) cat(hpds[j,k], "\t", file=hpd.file)
         cat("\n", file=hpd.file)
       }
       
       cat(paste("Mode\tMidRange\tMedian\n"), file=hpd.file)
       for(l in 1:ncol(pdf.info)) cat(pdf.info[l], "\t", file=hpd.file)
       cat("\n", file=hpd.file)
       
     }
     
     close(hpd.file)
    }
  
    }
  
  # export output as .csv file
  write.csv(output,paste("InputData/", name, "/", name,"_out", ext, sep=""), row.names = FALSE)
  
  #View(output)

  # Done message
  cat("Work has been done successfully, check outputs in your folder.\n\n")
}


# list the available data
Data <- list.files("InputData/")

#Welcome
cat("Hi, the ResAge package is here to help you with reservoir age calculation!\n")
cat("The rad2 function is designed to work with pairs of reservoir-derived and atmospheric 14C ages.\n\n")