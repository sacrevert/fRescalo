## Frescalo, based on Hill (2012) paper and original fortran code (https://www.brc.ac.uk/biblio/frescalo-computer-program-analyse-your-biological-records)
## O.L. Pescott, 30/06/2020
## v0.0 -- First pass
# Load example data used in https://github.com/BiologicalRecordsCentre/sparta
load(file = "data/unicorns.rda")
wgts <- read.delim(file = "data/GB_LC_Wts.txt", header = F, sep = "")

## Do the data contain locations that have no weights defined? Exclude them.
dat <- unicorns[unicorns$hectad %in% wgts$V1,] # lost 4 Locations
missingLocs <- unicorns[!(unicorns$hectad %in% wgts$V1),]

# useful lists
uniSpp <- unique(dat$CONCEPT) # unique species list
uniSites <- unique(dat$hectad) # unique data sites
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

##############################################
## Translation of fortran subroutine fresca ##
##############################################
# This function will be applied on a per site/neighbourhood basis (i.e. across i's)
krepmx <- 100 # max number of iterations
phibig <- 0.74 # target phi
tol <- 0.0003 # tolerance for approximation success
## While the absolute deviation of phi from target phi is > tolerance and iteration number <= 100
# continue with iterative approximation process
# lwfL[j,i] # species x site rescaled (and transformed) frequency matrix
alpha <- rep(1, length(uniSites))
phi <- rep(0, length(uniSites))
tot <- tot2 <- spnum <- rep(NA, length(uniSites))
for (i in 1:length(uniSites)) {
  while (isTRUE(abs(phi[i] - phibig) > tol)) {
    jDat <- jDat2 <- rep(NA, length(uniSpp))
    for (k in 1:krepmx) { # k iterations
      for (j in 1:length(uniSpp)) { # over j species in neighbourhood i
          jDat[j] <- lwfL[j,i]
          jDat2[j] <- 1-exp(-jDat[j]*alpha[i])
          tot[i] <- sum(jDat2)
          tot2[i] <- sum(jDat2^2)
        }
      # Values for neighbourhood i
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
}
#save(alpha, file = "outputs/alpha_v0.0.rda") # save initial versions for quick comparison in vignetteEg.R
#save(spnum, file = "outputs/spnum_v0.0.rda") # save initial versions for quick comparison in vignetteEg.R



##############################################
## Translation of fortran subroutine tfcalc ##
##############################################

