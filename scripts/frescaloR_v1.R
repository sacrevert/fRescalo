## Frescalo, based on Hill (2012) paper and original fortran code (https://www.brc.ac.uk/biblio/frescalo-computer-program-analyse-your-biological-records)
## O.L. Pescott, 30/06/2020
## v0.0 -- First pass
## v1.0 -- add time factor SDs
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
# For yearlsey comparison (proj 413)
#s <- dat[,c(3,5,7)]
#names(s) <- c("location", "species", "time")
#write.csv(s, file = "data/clusterTestDat.csv")

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
      #sum(siteNWgts[siteNWgts$V1==uniSites[i],]$V3+1.0E-10) # adding this small amount to avoid zeros is in original fortran code ## wrong
      #######################################################
      sum(siteNWgts[siteNWgts$V1==uniSites[i],]$V3)+1.0E-10 # this small amount is actually only added once in Hill -- June 2024 correction
      #######################################################
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
#save(lwf, file = "outputs/lwf_Jun24.rda") # raw weighted freqs -- with correction L65 - Jun 24
#save(lwfL, file = "outputs/lwfL_Jun24.rda") # neg log transformed weighted freqs -- with correction L65 - Jun 24
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
jDat <- jDat2 <- lwfRscd <- lwfRank1 <- lwfRanked <- bench <- benchName <- matrix(nrow = length(uniSpp), ncol = length(uniSites))
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
  lwfRank1 [,i] <- rank(lwfRanked[,i])
  lwfRscd[,i] <- rank(lwfRanked[,i])/spnum[i]
  # Then label all species as benchmarks or not within each neighbourhood
  for (j in 1:length(uniSpp)) {
    # indicate which species are benchmarks in each neighbourhood
    # Either species where scaled rank < blmdef or where rank = 1 even though scaled rank > blmdef (for small samples)
    #bench[j,i] <- ifelse(lwfRscd[j,i] < blmdef,1,0) # 1 = benchmark
    #bench[j,i] <- ifelse(lwfRscd[j,i] < blmdef | rank(lwfRanked[,i])[j] == 1,1,0) # 1 = benchmark
    ####################################################################################
    bench[j,i] <- ifelse(lwfRscd[j,i] < blmdef | rank(lwfRanked[,i])[j] == 1, 1, 1.0E-7) # Minor change to bring in line with fortran, Jun 24 (no effect on correlations with fortran in this e.g.)
    ####################################################################################
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
#save(alpha, file = "outputs/alpha_vJun24.rda") # with L65 June 24 correction
#save(spnum, file = "outputs/spnum_vJun24.rda") # with L65 June 24 correction
#save(bench, file = "outputs/bench_vJun24.rda") # with L65 June 24 correction
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
plog <- pfac <- estval <- array(data = 0, dim = c(length(uniSpp), length(uniSites), nrow(periods))) # [j,i,t]
#plog <- pfac <- estval <- diff <- sptot <- matrix(nrow = nrow(periods), ncol = length(uniSpp)) # [t,j]
estvar <- esttot <- sptot1 <- sptot <- matrix(data = 0, nrow = nrow(periods), ncol = length(uniSpp)) # [t,j]
tf <- matrix(data = 1, nrow = nrow(periods), ncol = length(uniSpp)) # [t,j]
StDev <- tfSD <- tf1 <- matrix(data = NA, nrow = nrow(periods), ncol = length(uniSpp)) # [t,j]
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
    estvar[t,j] <- estvar[t,j] + sum(wgt[t,] * wgt[t,] * estval[j,,t] * (1-estval[j,,t]))
    sptot1[t,j] <- sptot[t,j] + sqrt(estvar[t,j])
# Calculate time factor SDs
    for (i in 1:length(uniSites)) {
      # No point redoing this step, use existing objects
      #pfac[j,i,t] <- jDat2[j,i] * sampint[t,i] 
      #if(pfac[j,i,t] > 0.98) {pfac[j,i,t] <- 0.98} # avoid pfac == 1 situation
      #plog[j,i,t] <- -log(1-pfac[j,i,t])
      # calculate initial estimated values across sites
      estval[j,,t] <- 1 - exp(-plog[j,,t] * tf1[t,j]) # start with previously estimated time factors
      # sum to get total across sites (with down-weighting of under-sampled neighbourhoods)
      esttot[t,j] <- sum(wgt[t,] * estval[j,,t])
      # compare with empirical value sptot1 (1 SD from sptot)
      tfSD[t,j] <- sptot1[t,j]/(esttot[t,j]+0.0000001)
      # then iterate to find value of time factor that adjusts sptot1 to esttot
      # iterate as long as sptot1 and esttot are divergent within a t/j combination
      for (k in 1:krepTf) {
        if(abs(sptot1[t,j]-esttot[t,j]) > 0.0005) {
          estval[j,,t] <- 1 - exp(-plog[j,,t] * tfSD[t,j])
          esttot[t,j] <- sum(wgt[t,] * estval[j,,t])
          tfSD[t,j] <- tfSD[t,j] * (sptot1[t,j]/(esttot[t,j]+0.0000001))
        }
      }
      StDev[t,j] <- tfSD[t,j] - tf1[t,j]
    }
  }
}
## END

########################################
# Sparta/original Frescalo comparisons #
########################################
# library(sparta)
library(reshape2)
# myFresPath <- 'C:\\analyses\\Frescalo_3a_windows.exe' # old laptop path
# myFresPath <- 'C:\\Analyses\\Frescalo_3a_windows.exe' # new laptop June 2024
# myFolder = 'C:\\Users\\olipes\\Desktop\\frescaOLP\\outputs\\fresVignette'
# myTimePeriods2 <- data.frame(start = c(1980, 1990), end = c(1989,1999))
unicorn_TF <- frescalo(Data = unicorns,
                       frespath = myFresPath,
                       Fres_weights = 'LCGB', # British Land Cover Map weights
                       time_periods = myTimePeriods2,
                       site_col = 'hectad',
                       sp_col = 'CONCEPT',
                       start_col = 'TO_STARTDATE',
                       end_col = 'Date',
                       #sinkdir = myFolder)
                       sinkdir = "C:\\Analyses\\test")
# head(unicorn_TF$trend)
#save(unicorn_TF, file = "outputs/unicorn_TF.rda")
load(file = "outputs/unicorn_TF.rda")

#####################
## Compare results ##
#####################
library(useful)
library(reshape2)
#par(mfrow=c(1,2))
par(mfrow=c(2,4))
## Site-based stuff ##
# Alpha
#par(mfrow=c(1,1))
head(unicorn_TF$stat)
topleft(alpha)
alphaDF <- data.frame(alpha = alpha, Location = siteDF$sites)
alphaCompare <- merge(unicorn_TF$stat, alphaDF, by = c("Location"))
head(alphaCompare)
cor(alphaCompare$alpha, alphaCompare$Alpha) # 0.993 -- original; Jun 24 correction to L65: 0.9998
plot(alphaCompare$alpha, alphaCompare$Alpha, main = "Alpha") # 0.993
#plot(alphaCompare$alpha, alphaCompare$Alpha, main = "Alpha", xlim = c(15,30), ylim = c(15,30)) # 0.993
abline(a = 0, b = 1)

# Estimated species richness
#par(mfrow=c(1,1))
head(unicorn_TF$stat)
topleft(spnum)
spnumDF <- data.frame(spnum = spnum, Location = siteDF$sites)
spnumCompare <- merge(unicorn_TF$stat, spnumDF, by = c("Location"))
head(spnumCompare)
cor(spnumCompare$spnum, spnumCompare$Spnum_out) # 0.998 -- original; Jun 24 correction to L65: 0.9999
plot(spnumCompare$spnum, spnumCompare$Spnum_out, main = "Predicted \n site richness")
plot(spnumCompare$spnum, spnumCompare$Spnum_out, main = "Predicted \n site richness", 
     xlim = c(15,18), ylim = c(15,18))
abline(a = 0, b = 1)

## Species x site based stuff ##
# Local frequency
#par(mfrow=c(1,1))
head(unicorn_TF$freq)
topleft(lwf)
dimnames(lwf)[[1]] <- sppDF$species
dimnames(lwf)[[2]] <- siteDF$sites
lwfDF <- melt(lwf); names(lwfDF)[1:2] <- c("Species", "Location")
lwfCompare <- merge(unicorn_TF$freq, lwfDF, by = c("Species", "Location"))
head(lwfCompare)
cor(lwfCompare$value, lwfCompare$Freq) # 0.998 -- original; Jun 24 correction to L65: 0.9997
plot(lwfCompare$value, lwfCompare$Freq, main = "Raw species' weighted \n freqs") # 0.998
abline(a = 0, b = 1)

# Ranks
#par(mfrow=c(1,1))
dimnames(lwfRank1)[[1]] <- sppDF$species
dimnames(lwfRank1)[[2]] <- siteDF$sites
lwfRank1DF <- melt(lwfRank1); names(lwfRank1DF)[1:2] <- c("Species", "Location")
lwfRank1Compare <- merge(unicorn_TF$freq, lwfRank1DF, by = c("Species", "Location"))
head(lwfRank1Compare)
cor(lwfRank1Compare$value, lwfCompare$Rank) #
plot(lwfRank1Compare$value, lwfCompare$Rank, main = "Species ranks") # inflated ranks presumablt reflected in rescaled rank below
abline(a = 0, b = 1)
# Discrepancies
lwfRank1Compare$diff <- abs(lwfRank1Compare$Rank - lwfRank1Compare$value)
lwfRank1Compare[order(lwfRank1Compare$diff, decreasing = T),]
sum(lwfRank1Compare$diff==0)/nrow(lwfRank1Compare)*100 # 92.8% agreement

# Rescaled rank
#par(mfrow=c(1,1))
dimnames(lwfRscd)[[1]] <- sppDF$species
dimnames(lwfRscd)[[2]] <- siteDF$sites
lwfRscdDF <- melt(lwfRscd); names(lwfRscdDF)[1:2] <- c("Species", "Location")
lwfRscdCompare <- merge(unicorn_TF$freq, lwfRscdDF, by = c("Species", "Location"))
head(lwfRscdCompare)
cor(lwfRscdCompare$value, lwfCompare$Rank1) # 0.994 -- original; Jun 24 correction to L65: 0.9987
plot(lwfRscdCompare$value, lwfCompare$Rank1, main = "Rescaled species \nranks (R_ij)") # 0.994
abline(a = 0, b = 1)
# Discrepancies
lwfRscdCompare$diff <- abs(lwfRscdCompare$Rank1 - lwfRscdCompare$value)
lwfRscdCompare[order(lwfRscdCompare$diff, decreasing = T),]

# Rescaled freq
#par(mfrow=c(1,1))
dimnames(jDat2)[[1]] <- sppDF$species
dimnames(jDat2)[[2]] <- siteDF$sites
jDat2DF <- melt(jDat2); names(jDat2DF)[1:2] <- c("Species", "Location")
jDat2DFCompare <- merge(unicorn_TF$freq, jDat2DF, by = c("Species", "Location"))
head(jDat2DFCompare)
cor(jDat2DFCompare$value, jDat2DFCompare$Freq1) # 0.996 -- original; Jun 24 correction to L65: 0.9994
plot(jDat2DFCompare$value, jDat2DFCompare$Freq1, main = "Rescaled species' \n weighted frequencies (f_ij)") # 0.996
abline(a = 0, b = 1)

## Species x time stuff ##
# Time factor [t,j]
#par(mfrow=c(1,1))
dimnames(tf1)[[1]] <-  c(1984.5,1994.5)
dimnames(tf1)[[2]] <-  sppDF$species
tf1DF <- melt(tf1); names(tf1DF)[1:2] <- c("Time", "Species")
tf1DFCompare <- merge(unicorn_TF$trend, tf1DF, by = c("Time", "Species"))
head(tf1DFCompare)
tf1DFCompare <- tf1DFCompare[order(tf1DFCompare$TFactor),]
cor(tf1DFCompare$value, tf1DFCompare$TFactor) # 0.974 -- original; June 24, L65 correction: 0.9984
plot(tf1DFCompare$value, tf1DFCompare$TFactor, main = "Time factors") # 
abline(a = 0, b = 1)
# Discrepancies
tf1DFCompare$diff <- tf1DFCompare$TFactor - tf1DFCompare$value
head(tf1DFCompare[order(abs(tf1DFCompare$diff), decreasing = T),])

# Estimated species total (across sites) after rescaling
#par(mfrow=c(1,1))
dimnames(esttot)[[1]] <-  c(1984.5,1994.5)
dimnames(esttot)[[2]] <-  sppDF$species
esttotDF <- melt(esttot); names(esttotDF)[1:2] <- c("Time", "Species")
esttotDFCompare <- merge(unicorn_TF$trend, esttotDF, by = c("Time", "Species"))
head(esttotDFCompare)
esttotDFCompare <- esttotDFCompare[order(esttotDFCompare$TFactor),]
cor(esttotDFCompare$value, esttotDFCompare$Xest) # 0.997 -- orig; Jun 24, L65 correction: 0.9995
plot(esttotDFCompare$value, esttotDFCompare$Xest, main = "Est. per species total \nacross neighbourhoods \n(Sigma_i P_ijt)") # 
abline(a = 0, b = 1)
plot(esttotDFCompare$value, esttotDFCompare$Xest, main = "Est. per species total \nacross neighbourhoods \n(Sigma_i P_ijt)",
     xlim = c(0,50), ylim = c(0,50)) # 
abline(a = 0, b = 1)

### Added for v1
# Time factor SDs [t,j]
#par(mfrow=c(1,1))
dimnames(StDev)[[1]] <-  c(1984.5,1994.5)
dimnames(StDev)[[2]] <-  sppDF$species
StDevDF <- melt(StDev); names(StDevDF)[1:2] <- c("Time", "Species")
StDevDFCompare <- merge(unicorn_TF$trend, StDevDF, by = c("Time", "Species"))
head(StDevDFCompare)
StDevDFCompare <- StDevDFCompare[order(StDevDFCompare$TFactor),]
cor(StDevDFCompare$value, StDevDFCompare$StDev) # 0.999 -- orig; Jun 24, L65 correction: 0.9990576
plot(StDevDFCompare$value, StDevDFCompare$StDev, main = "Time factor SDs") # 
abline(a = 0, b = 1)
# Look at SD for spp with biggest TF discrepancy
StDevDFCompare[StDevDFCompare$Species == "Species 18" & StDevDFCompare$Time == 1984.5,]
# General discrepancies
StDevDFCompare$diff <- StDevDFCompare$StDev - StDevDFCompare$value
head(StDevDFCompare[order(abs(StDevDFCompare$diff), decreasing = T),])

## END
