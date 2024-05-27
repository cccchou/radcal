# ResAge R code for reservoir age offset calculation. Script radcal.R
# Guillaume Soulet <gsoulet@whoi.edu>
# see accompanying manual and Soulet 2015 (Quaternary Geochronology 29:97-103, doi: 10.1016/j.quageo.2015.05.023)
# updated with reservoir 14C disequilibria F14R (f.14R), delta14R (delta.14R) and dcp (dcp) 

radcal <- function(name="example_radcal_Southon2012", cc=1, cc1="IntCal13.14C", cc2="SHCal13.14C", cc3="IntCAL09.14C", ext=".csv", sep=",", yrsteps=1, uncalsteps=5, radsteps=1, convolsteps=1, prob=0.95, threshold=1e-6, hpdsteps=1, storedat=FALSE, export.uncal.pdf=FALSE, mixture.pdf=FALSE, delta.14R=FALSE, dcp=FALSE, f.14R=FALSE)
  .radcal(name, cc, cc1, cc2, cc3, ext, sep, yrsteps, uncalsteps, radsteps, convolsteps, prob, threshold, hpdsteps, storedat, export.uncal.pdf, mixture.pdf, delta.14R, dcp, f.14R)

.radcal <- function(name, cc, cc1, cc2, cc3, ext, sep, yrsteps, uncalsteps, radsteps, convolsteps, prob, threshold, hpdsteps, storedat, export.uncal.pdf, mixture.pdf, delta.14R, dcp, f.14R)
{
  dets <- suppressWarnings(read.csv(paste("InputData/", name, "/", name, ext, sep=""), sep=sep))
  
  # check for consistency 
  if(ncol(dets) != 5)
    stop("Your input file should contain 5 columns.\n Please, check the template designed for the radcal function.", call.=FALSE)
  
  if(min(dets[2]) < 0 || min(dets[4]) < 0)
    stop("At least one of your radiocarbon data or one of your calendar age is negative.\n ResAge is not designed to work with post-bomb samples.\n Remove such data from the input file.", call.=FALSE)
  
  if(delta.14R & dcp || delta.14R & f.14R || dcp & f.14R || dcp & f.14R & delta.14R)
    stop("Please, choose only one metric: reservoir age offset (default), f.14R, delta.14R, or dcp")
  
  # prepare for reservoir age calculation 
  id <- dets$id
  rad_res <- dets$Reservoir_derived_14C_age
  rad_res_sd <- dets$Reservoir_derived_14C_error
  cal.age <- dets$Calendar_age
  cal.sigma <- dets$Calendar_age_error
  
########################################################################################################################    
# the function to uncalibrate calendar ages
.uncaldist <- function(cal.age, cal.sigma, theta, f.mu, f.sigma, yrsteps, threshold, uncalsteps, radsteps)
  {
    # Let's reduce calculation time
    
    # 1) work over 6 time sigma of the sample calendar age 
    cal.span <- theta[which((theta >= cal.age - 6*cal.sigma) & (theta <= cal.age + 6*cal.sigma))]
    
    # 2) select to suitable range of the radiocarbon time scale to calibrate and convert it into F14C
    rad.span <- round(seq(min(rad[cal.span+1]) - 50*max(rad.sigma[cal.span+1]), max(rad[cal.span+1]) + 50*max(rad.sigma[cal.span+1]), by=uncalsteps))
    f.span <- exp(-rad.span/8033)
    
    # find how far each calendar time of cal.span is from cal.age (calendar measurement), this is also the normal pdf of cal.age over cal.span
    LKH.cal.age <- cbind(cal.span, dnorm(cal.span, cal.age, cal.sigma))
    
    # Uncalibrate !
    uncal <- c()

    for(i in 1:length(f.span))
      {
        # How far is a single F14C from the calibration curve f.mu  
        LKH.cc <- dnorm(f.mu[cal.span+1], f.span[[i]], f.sigma[cal.span+1]) 
        LKH.ri <- yrsteps*sum(LKH.cal.age[,2]*LKH.cc)
        uncal[[i]] <- LKH.ri
      }
    
    uncal <- cbind(rad.span, uncal)
    
    # interpolate and normalise uncalibrated distribution to 1
    uncal <- uncal[min(which(uncal[,2] > 0)):max(which(uncal[,2] > 0)),] # remove unnecessary data 
    uncal <- approx(uncal[,1], uncal[,2], seq(min(unlist(uncal[,1])), max(unlist(uncal[,1])), by=radsteps))
    uncal <- cbind(uncal$x, uncal$y/sum(uncal$y)) # normalise
      
    # only report those normalised uncalibrated probabilities beyond a threshold
    uncal[uncal[,2] > threshold,]
  }
  
########################################################################################################
# Find reservoir age: subtract uncalibrated age to reservoir-derived 14C age
.convolution <- function(dat, rad_res, rad_res_sd, threshold, convolsteps, radsteps)
{
  # prepare for convolution
  # create the normal distribution of the reservoir-derived 14C age
  spl.span <- seq(rad_res - 10*rad_res_sd, rad_res + 10*rad_res_sd, by=radsteps)
  spl.pdf <- cbind(spl.span, dnorm(spl.span, rad_res, rad_res_sd))
  spl.pdf <- cbind(spl.pdf[,1], spl.pdf[,2]/sum(spl.pdf[,2]))
  spl.pdf <- spl.pdf[spl.pdf[,2] > threshold,]
  
  # move the end of matrix dat (uncalibrated age) to the start of reservoir 14C age spl.pdf
  new.origin <- min(spl.pdf[,1])-max(dat[,1])
  uncal.conv <- cbind(c(dat[,1]+new.origin-1), dat[,2])
  
  # make matrices same dimension so that they overlap
  m <- nrow(spl.pdf)  
  n <- nrow(uncal.conv)
  
  spl.pdf.new <- rbind(cbind(uncal.conv[,1], rep(0,n)), spl.pdf)
  uncal.conv.new <- rbind(uncal.conv, cbind(spl.pdf[,1], rep(0,m)))

  spl.density <- spl.pdf.new[,2]
  uncal.density <- uncal.conv.new[,2]
  
  # convolution !
  
  LKH.res <- c()

  for(i in seq(1, (n+m+n), by=convolsteps))
    {
      spl.density <- c(spl.density, rep(0, convolsteps))
      uncal.density <- c(rep(0, convolsteps), uncal.density)
      LKH.res[[i]] <- sum(spl.density*uncal.density)    
    }

  res.age <- cbind(c(new.origin - 1 + seq(1,(n+m+n), by=convolsteps)), LKH.res)
  res.age <- res.age[min(which(res.age[,2] > 0)):max(which(res.age[,2] > 0)),] # remove unnecessary data 
  
  
  # only report those normalised calibrated probabilities beyond a threshold
  res.age[res.age[,2] > threshold,]
}


############################################################################################################################
# find the highest posterior density (hpd) of the calibrated distribution
.hpd <- function(dat, prob, hpdsteps, yrsteps)
{
  # interpolate and rank the ages according to their calibrated distribution probabilities
  dat <- approx(dat[,1], dat[,2], seq(min(unlist(dat[,1])), max(unlist(dat[,1])), by=yrsteps))
  o <- order(dat$y, decreasing=TRUE)
  dat <- cbind(dat$x[o], dat$y[o]/sum(dat$y))
  
  # only retain those ages with cumulative normalised probabilities within required percentage
  dat <- dat[which(cumsum(dat[,2]) <= prob),]
  dat <- dat[order(dat[,1]),]
  
  # identify any individual ranges within the hpd range and calculate their probability
  dif <- which(diff(dat[,1]) > hpdsteps)
  
  if(length(dif)==0)
    hpds <- cbind(min(unlist(dat[,1])), max(unlist(dat[,1])), 100*prob) else
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


#######################################################################################################
# Now, let's do the job now


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

cat("Uncalibrating calendar ages using the", ccname, "calibration curve and calculating reservoir age offsets...\n")

# prepare calcurve for calculation in F14C and interpolate
theta <- calcurve[,1]
f.mu <- exp(-calcurve[,2]/8033)
f.sigma <- exp(-(calcurve[,2]-calcurve[,3])/8033) - f.mu

# interpolate calibration each yrsteps (required for uncalibration function)
theta <- seq(0, max(theta), by=yrsteps)
rad <- approx(calcurve[,1], calcurve[,2], theta)$y
rad.sigma <- approx(calcurve[,1], calcurve[,3], theta)$y
f.mu <- approx(calcurve[,1], f.mu, theta)$y
f.sigma <- approx(calcurve[,1], f.sigma, theta)$y

dat <- data.frame(cbind(cal.age, cal.sigma, rad_res, rad_res_sd), row.names = NULL)

for(i in 1:length(cal.age))
{ cat(".")
  uncalib <- .uncaldist(dat$cal.age[[i]], dat$cal.sigma[[i]], theta, f.mu, f.sigma, yrsteps, threshold, uncalsteps, radsteps)
  dat$uncalib[[i]] <- uncalib # resulting uncalibrated (atm radiocarbon distribution) distribution
  dat$hpd_uncal[[i]] <- .hpd(uncalib, prob=prob, hpdsteps=hpdsteps, yrsteps=radsteps) # uncalibrated highest probability density corresponding to prob
  dat$mid1_uncal[[i]] <- round((dat$hpd_uncal[[i]][1] + dat$hpd_uncal[[i]][2*nrow(dat$hpd_uncal[[i]])])/2) # midpoints of uncalibrated ranges
  radyrs <- uncalib[,1]
  dat$mid2_uncal[[i]] <- round(mean(c(max(unlist(radyrs)), min(unlist(radyrs))))) # midpoints of entire uncalibrated distributions (with probabilities beyond threshold) 
  dat$wmn_uncal[[i]] <- round(weighted.mean(uncalib[,1], 1/uncalib[,2])) #round(weighted.mean(uncalib[,1], 1/uncalib[,2])) # weighted means of uncalibrated ranges
  dat$med_uncal[[i]] <- uncalib[max(which(cumsum(uncalib[,2]) <= .5)),1] # medians of uncalibrated distributions
  dat$mode_uncal[[i]] <- uncalib[which(uncalib[,2] == max(uncalib[,2])),1][1] #maximum densities of uncalibrated distributions
  
    
  resage <- .convolution(uncalib, dat$rad_res[[i]], dat$rad_res_sd[[i]], threshold, convolsteps, radsteps) 
  dat$resage[[i]] <- resage # resulting reservoir age 
  dat$hpd_resage[[i]] <- .hpd(resage, prob=prob, hpdsteps=hpdsteps, yrsteps=radsteps) # reservoir age highest probability density corresponding to prob
  dat$mid1_resage[[i]] <- round((dat$hpd_resage[[i]][1] + dat$hpd_resage[[i]][2*nrow(dat$hpd_resage[[i]])])/2) # midpoints of reservoir ranges
  resyrs <- resage[,1]
  
  dat$mid2_resage[[i]] <- round(mean(c(max(unlist(resyrs)), min(unlist(resyrs))))) # midpoints of entire reservoir age distributions (with probabilities beyond threshold) 
 
  #print(1/resage[,2])
  dat$wmn_resage[[i]] <- round(weighted.mean(unlist(resage[,1]), 1/unlist(resage[,2]))) # weighted means of reservoir ranges
  dat$med_resage[[i]] <- resage[max(which(cumsum(resage[,2]) <= .5)),1] # medians of reservoir distributions
  dat$mode_resage[[i]] <- resage[which(resage[,2] == max(unlist(resage[,2]))),1][1] #maximum densities of reservoir distributions
}

if(storedat) dat <<- dat

cat("\n")

#######################################################################################################
# export the uncalibrated probability density if you required it
if(export.uncal.pdf)
{
  cat("Exporting uncalibrated ages...\n")
  
  # find the min and max uncal ages among all the uncal age density probability
  mini <- c()
  maxi <- c()
  
  for(i in 1:nrow(dat))
  {
    mini[[i]] <- min(as.data.frame(dat$uncal[[i]])[,1])
    maxi[[i]] <- max(as.data.frame(dat$uncal[[i]])[,1])
  }
  min <- min(mini) 
  max <- max(maxi)
  
  # create a matrix of zeros and paste inside all the uncal age densities 
  uncal.age <- matrix(0, ncol=1+nrow(dat), nrow=max-min+1)
  uncal.age[1:nrow(uncal.age), 1] <- seq(min,max, by=radsteps)
  
  header <- c()
  
  
  for(i in 1:nrow(dat))
  {
    d <- max(diff(as.data.frame(dat$uncal[i])[,1]))
    if (d > 1)
      stop("Please, decrease threshold parameter: e.g. threshold = 1e-14 in the comand line\n\n", call.=FALSE)
    
    uncal.age[(min(as.data.frame(dat$uncal[[i]])[,1])-min+1):(max(as.data.frame(dat$uncal[[i]])[,1])-min+1), i+1] <- as.data.frame(dat$uncal[[i]])[,2]
    header[[i]] <- matrix(id[i])
  }
  header <- header
  colnames(uncal.age) <- c("Uncalibrated_age_14C_BP", header)
  
  write.csv(uncal.age,paste("InputData/", name, "/", name,cc,"_","uncalibrated_ages_pdfs", ext, sep=""), row.names = FALSE)
  
  
  # export uncalibrated ranges of all dates
  
  for(i in 1:nrow(dat))
  {
    hpds <- dat$hpd_uncal[[i]]
    if(nrow(hpds) == 1) test <- 1 else test <- 2
    if(test == 2) break
  }
  
  if(test==1) 
  {
    export.uncalib <- matrix(0, ncol=5, nrow=nrow(dat))
    for(i in 1:nrow(dat)) 
    {
      hpds <- dat$hpd_uncal[[i]]
      export.uncalib[i,1] <- hpds[,1]
      export.uncalib[i,2] <- hpds[,2]
      export.uncalib[i,3] <- dat$mode_uncal[[i]]
      export.uncalib[i,4] <- dat$mid1_uncal[[i]] 
      export.uncalib[i,5] <- dat$med_uncal[[i]]
      
    }
    
    colnames(export.uncalib) <- c(paste("UncalibratedMinRange","_at_",100*prob,"%", sep=""),paste("UncalibratedMaxRange","_at_",100*prob,"%", sep=""), "Mode", "MidRange", "Median") 
    output <- cbind(as.data.frame(id), export.uncalib)
    write.csv(output,paste("InputData/", name, "/", name,cc,"_","Uncalibrated_age_ranges_at", 100*prob, "%", ext, sep=""), row.names = FALSE)
    
  } else
  {
  
  hpd.file <- file(paste("InputData/", name, "/", name,cc,"_uncal_ranges.txt", sep=""), "w")
  cat(paste("Uncalibrated 14C age ranges at ", 100*prob, "% confidence intervals\n", sep=""), file=hpd.file)
  
  for(i in 1:nrow(dat))
  {
    cat(paste("\n\nIdentification_number: ", dets$id[[i]], "\nMinRange\tMaxRange\tprobability\n"), file=hpd.file)
    hpds <- dat$hpd_uncal[[i]]
    uncal.pdf.info <- cbind(dat$mode_uncal[[i]], dat$mid1_uncal[[i]], dat$med_uncal[[i]])
    
    for(j in 1:nrow(hpds))
    {
      for(k in 1:3) cat(hpds[j,k], "\t", file=hpd.file)
      cat("\n", file=hpd.file)
    }
    
    cat(paste("Mode\tMidRange\tMedian\n"), file=hpd.file)
    for(l in 1:ncol(uncal.pdf.info)) cat(uncal.pdf.info[l], "\t", file=hpd.file)
    cat("\n", file=hpd.file)
    
  }
  
  close(hpd.file)
  }
}

###########################################################################################################
# export the reservoir age distributions
cat("Exporting reservoir age offsets...\n")

# find the min and max reservoir ages among all the reservoir age density probability
mini <- c()
maxi <- c()
  
for(i in 1:nrow(dat))
{
  mini[[i]] <- min(unlist(as.data.frame(dat$resage[[i]])[,1]))
  maxi[[i]] <- max(unlist(as.data.frame(dat$resage[[i]])[,1]))
}
min <- min(unlist(mini)) 
max <- max(unlist(maxi))

  
# create a matrix of zeros and pastte inside all the reservoir age dennsities 
res.age <- matrix(0, ncol=1+nrow(dat), nrow=max-min+1)
#print(dim(res.age))
res.age[1:nrow(res.age), 1] <- seq(min,max, by=radsteps)

header <- c()

for(i in 1:nrow(dat))
{
  d <- max(diff(as.data.frame(dat$uncal[i])[,1]))
  
  if (d > 1)
    stop("Please, decrease threshold parameter: e.g threshold = 1e-14 in the command line\n", call.=FALSE)
  
  #print(as.data.frame(dat$resage[[i]])[,2])
  a=min(unlist(as.data.frame(dat$resage[[i]])[,1]))-min+1
  b=max(unlist(as.data.frame(dat$resage[[i]])[,1]))-min+1
  res.age[a:b, i+1] <- unlist(as.data.frame(dat$resage[[i]])[,2])
  
  header[[i]] <- matrix(id[i])
  
  
}


# in case user wants the mixture pdf
if(mixture.pdf) 
{
  mix <- cbind(rowSums(res.age[,2:ncol(res.age)]))
  mix <- mix/sum(mix)
  mix <- cbind(res.age[,1], mix)

  hpd.mix <- .hpd(mix, prob=prob, hpdsteps=hpdsteps, yrsteps=radsteps)
  
  mid1_mix <- round((hpd.mix[1] + hpd.mix[2*nrow(hpd.mix)])/2) # midpoints of reservoir ranges
  resyrs <- mix[,1]
  mid2_mix <- round(mean(c(max(resyrs), min(resyrs)))) # midpoints of entire reservoir age distributions (with probabilities beyond threshold) 
  wmn_mix <- round(weighted.mean(mix[,1], 1/mix[,2])) # weighted means of reservoir ranges
  med_mix <- mix[max(which(cumsum(mix[,2]) <= .5)),1] # medians of reservoir distributions
  mode_mix <- mix[which(mix[,2] == max(mix[,2])),1][1] #maximum densities of reservoir distributions
  
  mix.pdf.info <- cbind(mode_mix, mid1_mix, med_mix)
  
  res.age <- cbind(res.age, mix[,2])
  header <- c(header, "mixture_pdf")
  
  ref <- "Mixture_reservoir_age_offset_"
  
  if(f.14R)
  {
    mix.pdf.info <- round(exp(mix.pdf.info/-8033),4)
    hpd.mix.save <- hpd.mix[,1]
    hpd.mix[,1] <- round(exp(hpd.mix[,2]/-8033),4)
    hpd.mix[,2] <- round(exp(hpd.mix.save/-8033),4)
    ref <- "Mixture_F14R_percent_"
  }
    
  if(delta.14R)
  {
    mix.pdf.info <- round(1000*(exp(mix.pdf.info/-8033)-1),1)
    hpd.mix.save <- hpd.mix[,1]
    hpd.mix[,1] <- round(1000*(exp(hpd.mix[,2]/-8033)-1),1)
    hpd.mix[,2] <- round(1000*(exp(hpd.mix.save/-8033)-1),1)
    ref <- "Mixture_delta14R_permil_"
  }
  
  if(dcp)
  {
    mix.pdf.info <- round(100*(1-exp(mix.pdf.info/-8033)),2)
    hpd.mix[,1:2] <- round(100*(1-exp(hpd.mix[,1:2]/-8033)),2)
    ref <- "Mixture_dcp_percent_"
  }

  # export file with information about the mixture pdf
  
  hpd.file <- file(paste("InputData/", name, "/", name,cc,"_", ref,"Ranges.txt", sep=""), "w")
  cat(paste(ref, "at ", 100*prob, "% confidence intervals\n", sep=""), file=hpd.file)
  cat(paste("\n\nMixture_pdf: ", "\nMinRange\tMaxRange\tprobability\n"), file=hpd.file)
    
  
  for(j in 1:nrow(hpd.mix))
  {
    for(k in 1:3) cat(hpd.mix[j,k], "\t", file=hpd.file)
    cat("\n", file=hpd.file)
  }
  
  cat(paste("Mode\tMidRange\tMedian\n"), file=hpd.file)
  for(l in 1:ncol(mix.pdf.info)) cat(mix.pdf.info[l], "\t", file=hpd.file)
  cat("\n", file=hpd.file)
  
  close(hpd.file)
  
}

  header <- header



# Depending on Reservoir age offset, F14R, delta14R or dcp metrics

  ref <- "Reservoir_14C_years"
  ref2 <- "Reservoir_age_offset_"

  if(delta.14R)
    {res.age[,1] <- 1000*(exp(res.age[,1]/-8033)-1)
    ref <- "delta14R_permil"
    ref2 <- "delta14R_"
    }

  if(f.14R)
    {res.age[,1] <- exp(res.age[,1]/-8033)
     ref <- "F14R"
     ref2 <- "F14R_"
    }
  
  if(dcp)
    {res.age[,1] <- 100*(1-exp(res.age[,1]/-8033))
    ref <- "dcp_percent"
    ref2 <- "dcp_"
    }


  colnames(res.age) <- c(ref, header)
  res.age <- res.age
      
  write.csv(res.age,paste("InputData/", name, "/", name,cc,"_",ref2,"pdfs", ext, sep=""), row.names = FALSE)


# Export Reservoir age ranges of all dates

for(i in 1:nrow(dat))
{
  hpds <- dat$hpd_resage[[i]]
  if(nrow(hpds) == 1) test <- 1 else test <- 2
  if(test == 2) break
}

if(test == 1)
{
   export.resage <- matrix(0, ncol=5, nrow=nrow(dat))
   for(i in 1:nrow(dat)) 
   {
     hpds <- dat$hpd_resage[[i]]
     export.resage[i,1] <- hpds[,1]
     export.resage[i,2] <- hpds[,2]
     export.resage[i,3] <- dat$mode_resage[[i]]
     export.resage[i,4] <- dat$mid1_resage[[i]] 
     export.resage[i,5] <- dat$med_resage[[i]]
      
     
     ref <- "Reservoir_Min_Range"
     ref2 <- "Reservoir_Max_Range"
     ref3 <- "Reservoir_age_offset_14C_yrs_"
     
     if(f.14R)
     {export.resage[i,1] <- round(exp(hpds[,2]/-8033),4)
      export.resage[i,2] <- round(exp(hpds[,1]/-8033),4)
      export.resage[i,3] <- round(exp(dat$mode_resage[[i]]/-8033),4)
      export.resage[i,4] <- round(exp(dat$mid1_resage[[i]]/-8033),4)
      export.resage[i,5] <- round(exp(dat$med_resage[[i]]/-8033),4)
      ref <- "F14R_Min_Range"
      ref2 <- "F14R_Max_Range"
      ref3 <- "F14R_"
     }
     
     if(delta.14R)
       {export.resage[i,1] <- round(1000*(exp(hpds[,2]/-8033)-1),1)
       export.resage[i,2] <- round(1000*(exp(hpds[,1]/-8033)-1),1)
       export.resage[i,3] <- round(1000*(exp(dat$mode_resage[[i]]/-8033)-1),1)
       export.resage[i,4] <- round(1000*(exp(dat$mid1_resage[[i]]/-8033)-1),1)
       export.resage[i,5] <- round(1000*(exp(dat$med_resage[[i]]/-8033)-1),1)
       ref <- "delta14R_Min_Range"
       ref2 <- "delta14R_Max_Range"
       ref3 <- "delta14R_permil_"
       }
     
     if(dcp)
     {export.resage[i,1] <- round(100*(1-exp(hpds[,1]/-8033)),2)
      export.resage[i,2] <- round(100*(1-exp(hpds[,2]/-8033)),2)
      export.resage[i,3] <- round(100*(1-exp(dat$mode_resage[[i]]/-8033)),2)
      export.resage[i,4] <- round(100*(1-exp(dat$mid1_resage[[i]]/-8033)),2)
      export.resage[i,5] <- round(100*(1-exp(dat$med_resage[[i]]/-8033)),2)
      ref <- "dcp_Min_Range"
      ref2 <- "dcp_Max_Range"
      ref3 <- "dcp_percent_"
      }
   }
   
   
   header <- c(paste(ref,"_at_",100*prob,"%", sep=""),paste(ref2,"_at_",100*prob,"%", sep=""), "Mode", "MidRange", "Median")
   colnames(export.resage) <- header
   output <- cbind(as.data.frame(id), as.data.frame(cal.age), as.data.frame(cal.sigma), export.resage)
   colnames(output) <- c("id","Calendar_age", "Calendar_age_error", header)
   write.csv(output,paste("InputData/", name, "/", name,cc,"_",ref3,"ranges_at", 100*prob, "%", ext, sep=""), row.names = FALSE)
   
} else
{
  ref <- "Reservoir_age_offset_14C_yrs_"
  
  if(f.14R)
    ref <- "F14R_"
  
  if(delta.14R)
    ref <- "delta14R_permil_"
  
  if(dcp)
    ref <- "dcp_percent_"
  
hpd.file <- file(paste("InputData/", name, "/", name,cc,"_",ref,"ranges.txt", sep=""), "w")
cat(paste(ref,"ranges at ", 100*prob, "% confidence intervals\n", sep=""), file=hpd.file)

for(i in 1:nrow(dat))
{
  cat(paste("\n\nIdentification_number: ", dets$id[[i]], "\nMinRange\tMaxRange\tprobability\n"), file=hpd.file)
  hpds <- dat$hpd_resage[[i]]
  pdf.info <- cbind(dat$mode_resage[[i]], dat$mid1_resage[[i]], dat$med_resage[[i]])
  
  if(f.14R)
  {
    hpds_save <- hpds[,1]
    hpds[,1] <- round(exp(hpds[,2]/-8033),4)
    hpds[,2] <- round(exp(hpds_save/-8033),4)
    pdf.info <- round(exp(pdf.info/-8033),4)
  }
  
  if(delta.14R)
  {
    hpds_save <- hpds[,1]
    hpds[,1] <- round(1000*(exp(hpds[,2]/-8033)-1),1)
    hpds[,2] <- round(1000*(exp(hpds_save/-8033)-1),1)
    pdf.info <- round(1000*(exp(pdf.info/-8033)-1),1)
  }
  
  if(dcp)
  {
    hpds[,1] <- round(100*(1-exp(hpds[,1]/-8033)),2)
    hpds[,2] <- round(100*(1-exp(hpds[,2]/-8033)),2)
    pdf.info <- round(100*(1-exp(pdf.info/-8033)),2)
  }
  
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

#####################################################################################################
# Done message
cat("Work has been done successfully, check outputs in your folder.\n\n")

}


# list the available data
Data <- list.files("InputData/")

#Welcome
cat("Hi, ResAge is here to help you with reservoir age offset calculation!\n")
cat("The radcal function is designed to work with pairs of reservoir-derived 14C age and corresponding calendar age.\n\n")