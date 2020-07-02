## Frescalo, based on Hill (2012) paper and original fortran code (https://www.brc.ac.uk/biblio/frescalo-computer-program-analyse-your-biological-records)
## O.L. Pescott, 30/06/2020
## v0.0 -- First pass
# Load example data used in https://github.com/BiologicalRecordsCentre/sparta
#rm(list=ls())
library(useful)
load(file = "data/unicorns.rda")
wgts <- read.delim(file = "data/GB_LC_Wts.txt", header = F, sep = "")
periods <- data.frame(start = c(1980, 1990), end = c(1989, 1999))
head(periods)

## Do the data contain locations that have no weights defined? Exclude them.
dat <- unicorns[unicorns$hectad %in% wgts$V1,] # lost 4 Locations
missingLocs <- unicorns[!(unicorns$hectad %in% wgts$V1),]

# useful lists etc.
uniSpp <- unique(dat$CONCEPT) # unique species list
sppDF <- data.frame(index = 1:length(uniSpp), species = uniSpp, fLevel = as.numeric(uniSpp))
uniSites <- unique(dat$hectad) # unique data sites
siteDF <- data.frame(index = 1:length(uniSites), sites = uniSites)
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
                           dat[dat$CONCEPT==uniSpp[j],]$hectad,]$V3)/
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
          jDat2[j,i] <- 1-exp(-jDat[j,i]*alpha[i])
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
    benchName[j,i] <- ifelse(bench[j,i] == 1, sppDF[sppDF$index==j,]$species, 0) # species as factor levels
  }
  # Total sum of benchmark spp at sites
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
iocc <- array(dim = c(length(uniSpp), length(uniSites), nrow(periods))) # [j,i,t]
sampint <- matrix(data = 1.0E-7, nrow = nrow(periods), ncol = length(uniSites)) # [t,i]
#sampef(i,iit)=sampef(i,iit)+ibench(i,j)*bwght(j)/abtot(i)
for (t in 1:nrow(periods)) {
  for (i in 1:length(uniSites)) {
    sampint[t,i] <- 
  }
}

