## Frescalo, based on Hill (2012) paper and original fortran code (https://www.brc.ac.uk/biblio/frescalo-computer-program-analyse-your-biological-records)
## O.L. Pescott, 30/06/2020
## v0.0 -- First pass
# Load example data used in https://github.com/BiologicalRecordsCentre/sparta
#rm(list=ls())
library(useful)
load(file = "data/unicorns.rda")
wgts <- read.delim(file = "data/GB_LC_Wts.txt", header = F, sep = "")
#periods <- data.frame(t = c(1,2), start = as.Date(c("1980-01-01", "1990-01-01"), "%Y-%m-%d"), 
#                      end = as.Date(c("1989-12-31", "1999-12-31"), "%Y-%m-%d"))
periods <- data.frame(start = c(1980, 1990), end = c(1989,1999))
head(periods)

unicorns$year <- format(unicorns$TO_STARTDATE, "%Y")
unicorns$year <- as.numeric(unicorns$year)

# Add time period classification to data
# data.table::foverlaps probably a lot quicker...
for (d in 1:nrow(unicorns)) {
  unicorns$period[d] <- if(unicorns$year[d] >= 1980 && unicorns$year[d] <= 1989) {
    1
  } else if (unicorns$year[d] >= 1990 && unicorns$year[d] <= 1999) {
    2
  } else {
    NA
  }
}

## Only keep relevant time period data
dat <- unicorns[which(!is.na(unicorns$period)),]

## Do the data contain locations that have no weights defined? Exclude them.
dat <- dat[dat$hectad %in% wgts$V1,] # lost 1 Location
dat$datIndex <- 1:nrow(dat)

# useful lists etc.
uniSpp <- unique(dat$CONCEPT) # unique species list
sppDF <- data.frame(spIndex = 1:length(uniSpp), species = uniSpp, fLevel = as.numeric(uniSpp))
datM <- merge(dat, sppDF, by.x = "CONCEPT", by.y = "species")
datM <- datM[order(datM$datIndex),]
uniSites <- unique(datM$hectad) # unique data sites
siteDF <- data.frame(siteIndex = 1:length(uniSites), sites = uniSites)
siteNWgts <- wgts[wgts$V1 %in% uniSites,] # just keep relevant neighbourhoods

##########################################
## Calculate local weighted frequencies ##
##########################################
lwf <- matrix(nrow = length(uniSpp), ncol = length(uniSites))
for (j in 1:length(uniSpp)){ # for each species j
 for (i in 1:length(uniSites)) { # in each hectad i
# Lookup the weights in a neighbourhood where the species is present & sum:
# sum(siteNWgts[siteNWgts$V1==uniSites[i] 
#             & siteNWgts$V2 %in% 
#               dat[dat$CONCEPT==uniSpp[j],]$hectad,]$V3)
# Divide by summed weights in all hectads in that neighbourhood:
# sum(siteNWgts[siteNWgts$V1==uniSites[i],]$V3)
   lwf[j,i] <- sum(siteNWgts[siteNWgts$V1==uniSites[i] 
                         & siteNWgts$V2 %in% 
                           datM[datM$CONCEPT==uniSpp[j],]$hectad,]$V3)/
               sum(siteNWgts[siteNWgts$V1==uniSites[i],]$V3+1.0E-10) # adding this small amount to avoid zeros is in original fortran code
 }
}
# Ad hoc checks
#i = 1
#j = 1
#uniSpp[j]; uniSites[1]
#test <- siteNWgts[siteNWgts$V1==uniSites[i],] # all neighbours of i
#sppTest <- dat[dat$CONCEPT =="Species 18",] # all hectads with j
#sppTest$hectad[sppTest$hectad %in% test$V2] # list of neighbours with species

## Apply fmin and fmax in prep for later -log transformation (fortran code lines 479-480)
lwf <- apply(lwf, 1:2, function(x) ifelse(x < 1.0E-10, 1.0E-10, x))
lwf <- apply(lwf, 1:2, function(x) ifelse(x > 0.99999, 0.99999, x))
# Apply negative log transformation here (fortran code line 481)
lwfL <- apply(lwf, 1:2, function(x) -log(1-x))
#save(lwf, file = "outputs/lwf.rda") # raw weighted freqs
#save(lwfL, file = "outputs/lwfL.rda") # neg log transformed weighted freqs
#load(file = "outputs/lwf.rda")

##############################################
## Translation of fortran subroutine fresca ##
##############################################
# This function will be applied on a per site/neighbourhood basis (i.e. across i's)
krepmx <- 100 # max number of iterations
phibig <- 0.74 # target phi
tol <- 0.0003 # tolerance for approximation success
blmdef <- 0.2703 # default limit for benchmark species
## While the absolute deviation of phi from target phi is > tolerance and iteration number <= 100
# continue with iterative approximation process
# lwfL[j,i] # species x site rescaled (and transformed) frequency matrix
alpha <- rep(1, length(uniSites))
phi <- rep(0, length(uniSites))
tot <- tot2 <- spnum <- an2 <- abtot <- rep(NA, length(uniSites))
jDat <- jDat2 <- lwfRscd <- lwfRanked <- bench <- benchName <- matrix(nrow = length(uniSpp), ncol = length(uniSites))
siteBench <- list()
for (i in 1:length(uniSites)) {
  while (isTRUE(abs(phi[i] - phibig) > tol)) {
    for (k in 1:krepmx) { # k iterations
      for (j in 1:length(uniSpp)) { # over j species in neighbourhood i
          jDat[j,i] <- lwfL[j,i]
          jDat2[j,i] <- 1-exp(-jDat[j,i]*alpha[i]) # jDat2 = rescaled species frequencies
        }
      # Values for site (& neighbourhood) i
      tot[i] <- sum(jDat2[,i])
      tot2[i] <- sum(jDat2[,i]^2)
      phi[i] <- tot2[i]/tot[i]
      an2[i] <- tot[i]^2/tot2[i]
      spnum[i] <- tot[i]
      if (k < 20){
        alpha[i] <- alpha[i]*exp(1.86*(log(1-phi[i])-log(1-phibig)))
      } else {
        alpha[i] <- alpha[i]*(phibig/phi[i])
      }
    }
  }
  # Now calculate rescaled (j) species ranks in (i) site neighbourhoods
  # in order to label local benchmark species for time factor calcs
  for (j in 1:length(uniSpp)) {
    lwfRanked[j,i] <- -jDat2[j,i] + j*1.0E-12 # a way of sorting out ties amongst zeros
    }
  # Rank across all species within a neighbourhood and then normalise by expected (predicted) species number
  lwfRscd[,i] <- rank(lwfRanked[,i])/spnum[i]
  # Then label all species as benchmarks or not within each neighbourhood
  for (j in 1:length(uniSpp)) {
  # indicate which species are benchmarks in each neighbourhood
  # Either species where scaled rank < blmdef or where rank = 1 even though scaled rank > blmdef (for small samples)
    #bench[j,i] <- ifelse(lwfRscd[j,i] < blmdef,1,0) # 1 = benchmark
    bench[j,i] <- ifelse(lwfRscd[j,i] < blmdef | rank(lwfRanked[,i])[j] == 1,1,0) # 1 = benchmark
    benchName[j,i] <- ifelse(bench[j,i] == 1, sppDF[sppDF$spIndex==j,]$species, 0) # species as factor levels
  }
  # Total sum of benchmark spp at site
  abtot[i] <- sum(bench[,i])
  # List of benchmark species' factor levels at each site
  siteBench[[i]] <- benchName[,i][benchName[,i]!=0] # list of benchmark species (as factor levels) at each site
}
#save(alpha, file = "outputs/alpha_v0.0.rda") # save initial versions for quick comparison in vignetteEg.R
#save(spnum, file = "outputs/spnum_v0.0.rda") # save initial versions for quick comparison in vignetteEg.R
#save(bench, file = "outputs/bench_v0.0.rda") # save initial versions for quick comparison in vignetteEg.R
#load(file = "outputs/alpha_v0.0.rda")
#load(file = "outputs/spnum_v0.0.rda")
#load(file = "outputs/bench_v0.0.rda")

##############################################
## Translation of fortran subroutine tfcalc ##
##############################################
# Need some species x site x time period arrays for this bit [j,i,t]
# Following adpated from original Frescalo README
# sptot - Number of occurrences, given reduced weight of locations having very low sampling effort
# jtot - Actual number of occurrences of species j at time t
# esttot - Estimated number of occurrences of species j at time t, given fitted model
# iocc - 1 if species j is found at location i at time t, 0 otherwise
# smpint - the benchmark-based sampling intensity at location i and time t
sampint <- wgt <- matrix(data = 1.0E-7, nrow = nrow(periods), ncol = length(uniSites)) # [t,i]

## Create the time period-specific sampling intensity measures based on benchmarks recorded in site i in time period t
for (t in 1:nrow(periods)) {
  for (i in 1:length(uniSites)) {
    tiBenchRecs <- datM[which(datM$period == t & datM$hectad == uniSites[i] 
                          & datM$fLevel %in% siteBench[[i]]),]
    # Note that this approach does not currently incorporate the option to downweight particular named benchmark species
    sampint[t,i] <- nrow(unique(tiBenchRecs[,c("fLevel","hectad")]))/abtot[i]
    # Now create benchmark/time period weights as described in Frescalo README (Note 3, L80)
    wgt[t,i] <- ifelse(sampint[t,i] < 0.0995, sampint[t,i]*10 + 0.005, 1)
  }
}

# Species j pres/abs at time t, site i
iocc <- array(dim = c(length(uniSpp), length(uniSites), nrow(periods))) # [j,i,t]
for (j in 1:length(uniSpp)) {
  for (i in 1:length(uniSites)) {
    for (t in 1:nrow(periods)) {
      iocc[j,i,t] <- ifelse(nrow(datM[which(datM$period == t & datM$hectad == uniSites[i] 
                                      & datM$spIndex == j),]) >= 1, 1, 0)
    }
  }
}

# Calculate time factors
## I think we actually need to remove the zeros otherwise the TFactors are not comparable with the MOH results
plog <- pfac <- estval <- array(data = 0, dim = c(length(uniSpp), length(uniSites), nrow(periods))) # [j,i,t]
#plog <- pfac <- estval <- diff <- sptot <- matrix(nrow = nrow(periods), ncol = length(uniSpp)) # [t,j]
esttot <- sptot <- matrix(data = 0, nrow = nrow(periods), ncol = length(uniSpp)) # [t,j]
tf <- matrix(data = 1, nrow = nrow(periods), ncol = length(uniSpp)) # [t,j]
tf1 <- matrix(data = NA, nrow = nrow(periods), ncol = length(uniSpp)) # [t,j]
krepTf <- 100
for (t in 1:nrow(periods)) {
  for (j in 1:length(uniSpp)) {
    # Actual per species total across sites in each time period
    sptot[t,j] <- sum(wgt[t,] * iocc[j,,t]) # (with downweighting of sites with very little sampling)
    for (i in 1:length(uniSites)) {
      # adjust rescaled all-time-period weighted freqs by benchmark-based site sampling intensity
      pfac[j,i,t] <- jDat2[j,i] * sampint[t,i] 
      if(pfac[j,i,t] > 0.98) {pfac[j,i,t] <- 0.98} # avoid pfac == 1 situation
      plog[j,i,t] <- -log(1-pfac[j,i,t])
      # calculate initial estimated values across sites
      estval[j,,t] <- 1 - exp(-plog[j,,t] * tf[t,j])
      # sum to get total across sites (with down-weighting of under-sampled neighbourhoods)
      esttot[t,j] <- sum(wgt[t,] * estval[j,,t])
      # compare with empirical value sptot
      tf1[t,j] <- sptot[t,j]/(esttot[t,j]+0.0000001)
      # then iterate to find value of time factor that adjusts sptot to esttot
      # iterate as long as sptot and esttot are divergent within a t/j combination
      for (k in 1:krepTf) {
        if(abs(sptot[t,j]-esttot[t,j]) > 0.0005) {
          estval[j,,t] <- 1 - exp(-plog[j,,t] * tf1[t,j])
          esttot[t,j] <- sum(wgt[t,] * estval[j,,t])
          tf1[t,j] <- tf1[t,j] * (sptot[t,j]/(esttot[t,j]+0.0000001))
        }
      }
    }
  }
}
## END
# Time factor standard deviations not implemented

########################################
# Sparta/original Frescalo comparisons #
########################################
library(sparta)
library(reshape2)
myFresPath <- 'C:\\Frescalo_3a_windows.exe'
myFolder = 'C:\\Users\\olipes\\Desktop\\frescaOLP\\outputs\\fresVignette'
myTimePeriods2 <- data.frame(start = c(1980, 1990), end = c(1989,1999))
unicorn_TF <- frescalo(Data = unicorns, 
                       frespath = myFresPath,
                       Fres_weights = 'LCGB', # British Land Cover Map weights
                       time_periods = myTimePeriods2,
                       site_col = 'hectad',
                       sp_col = 'CONCEPT',
                       start_col = 'TO_STARTDATE',
                       end_col = 'Date',
                       sinkdir = myFolder)
head(unicorn_TF$trend)

#####################
## Compare results ##
#####################
#par(mfrow=c(1,2))
par(mfrow=c(2,4))
## Site-based stuff ##
# Alpha
head(unicorn_TF$stat)
topleft(alpha)
alphaDF <- data.frame(alpha = alpha, Location = siteDF$sites)
alphaCompare <- merge(unicorn_TF$stat, alphaDF, by = c("Location"))
head(alphaCompare)
cor(alphaCompare$alpha, alphaCompare$Alpha) # 0.993
plot(alphaCompare$alpha, alphaCompare$Alpha, main = "Alpha") # 0.993
abline(a = 0, b = 1)

# Estimated species richness
head(unicorn_TF$stat)
topleft(spnum)
spnumDF <- data.frame(spnum = spnum, Location = siteDF$sites)
spnumCompare <- merge(unicorn_TF$stat, spnumDF, by = c("Location"))
head(spnumCompare)
cor(spnumCompare$spnum, spnumCompare$Spnum_out) # 0.998
plot(spnumCompare$spnum, spnumCompare$Spnum_out, main = "Predicted \n site richness")
abline(a = 0, b = 1)

## Species x site based stuff ##
# Local frequency
head(unicorn_TF$freq)
topleft(lwf)
dimnames(lwf)[[1]] <- sppDF$species
dimnames(lwf)[[2]] <- siteDF$sites
lwfDF <- melt(lwf); names(lwfDF)[1:2] <- c("Species", "Location")
lwfCompare <- merge(unicorn_TF$freq, lwfDF, by = c("Species", "Location"))
head(lwfCompare)
cor(lwfCompare$value, lwfCompare$Freq) # 0.998
plot(lwfCompare$value, lwfCompare$Freq, main = "Raw species' weighted \n freqs") # 0.998
abline(a = 0, b = 1)

# Local frequency (log transformation)
head(unicorn_TF$freq)
topleft(lwfL)
dimnames(lwfL)[[1]] <- sppDF$species
dimnames(lwfL)[[2]] <- siteDF$sites
lwfLDF <- melt(lwfL); names(lwfLDF)[1:2] <- c("Species", "Location")
lwfLCompare <- merge(unicorn_TF$freq, lwfLDF, by = c("Species", "Location"))
head(lwfLCompare)
cor(lwfLCompare$value, lwfLCompare$Freq) # 0.992
abline(a = 0, b = 1)

# Rescaled rank
dimnames(lwfRscd)[[1]] <- sppDF$species
dimnames(lwfRscd)[[2]] <- siteDF$sites
lwfRscdDF <- melt(lwfRscd); names(lwfRscdDF)[1:2] <- c("Species", "Location")
lwfRscdCompare <- merge(unicorn_TF$freq, lwfRscdDF, by = c("Species", "Location"))
head(lwfRscdCompare)
cor(lwfRscdCompare$value, lwfCompare$Rank1) # 0.994
plot(lwfRscdCompare$value, lwfCompare$Rank1, main = "Rescaled species \nranks (R_ij)") # 0.994
abline(a = 0, b = 1)

# Rescaled freq
dimnames(jDat2)[[1]] <- sppDF$species
dimnames(jDat2)[[2]] <- siteDF$sites
jDat2DF <- melt(jDat2); names(jDat2DF)[1:2] <- c("Species", "Location")
jDat2DFCompare <- merge(unicorn_TF$freq, jDat2DF, by = c("Species", "Location"))
head(jDat2DFCompare)
cor(jDat2DFCompare$value, jDat2DFCompare$Freq1) # 0.996
plot(jDat2DFCompare$value, jDat2DFCompare$Freq1, main = "Rescaled species' \n weighted frequencies (f_ij)") # 0.996
abline(a = 0, b = 1)

## Species x time stuff ##
# Time factor [t,j]
dimnames(tf1)[[1]] <-  c(1984.5,1994.5)
dimnames(tf1)[[2]] <-  sppDF$species
tf1DF <- melt(tf1); names(tf1DF)[1:2] <- c("Time", "Species")
tf1DFCompare <- merge(unicorn_TF$trend, tf1DF, by = c("Time", "Species"))
head(tf1DFCompare)
tf1DFCompare <- tf1DFCompare[order(tf1DFCompare$TFactor),]
cor(tf1DFCompare$value, tf1DFCompare$TFactor) # 0.974
plot(tf1DFCompare$value, tf1DFCompare$TFactor, main = "Time factors") # 
abline(a = 0, b = 1)

# Estimated species total (across sites) after rescaling
dimnames(esttot)[[1]] <-  c(1984.5,1994.5)
dimnames(esttot)[[2]] <-  sppDF$species
esttotDF <- melt(esttot); names(esttotDF)[1:2] <- c("Time", "Species")
esttotDFCompare <- merge(unicorn_TF$trend, esttotDF, by = c("Time", "Species"))
head(esttotDFCompare)
esttotDFCompare <- esttotDFCompare[order(esttotDFCompare$TFactor),]
cor(esttotDFCompare$value, esttotDFCompare$Xest) # 0.997
plot(esttotDFCompare$value, esttotDFCompare$Xest, main = "Est. per species total \nacross neighbourhoods \n(Sigma_i P_ijt)") # 
abline(a = 0, b = 1)

## END
