library(sparta)
myFresPath <- 'C:\\Frescalo_3a_windows.exe'
myTimePeriods <- data.frame(start = c(1800), end = c(2018))
head(myTimePeriods)
myFolder = 'C:\\Users\\olipes\\Desktop\\frescaOLP\\outputs\\fresVignette'
unicorn_results <- frescalo(Data = unicorns, 
                            frespath = myFresPath,
                            Fres_weights = 'LCGB',
                            time_periods = myTimePeriods,
                            site_col = 'hectad',
                            sp_col = 'CONCEPT',
                            start_col = 'TO_STARTDATE',
                            end_col = 'Date',
                            sinkdir = myFolder)

head(unicorn_results$stat)

## Compare with alpha from my R version 0.0
wgts <- read.delim(file = "data/GB_LC_Wts.txt", header = F, sep = "")
## Do the data contain locations that have no weights defined? Exclude them.
dat <- unicorns[unicorns$hectad %in% wgts$V1,] # lost 4 Locations
load(file = "outputs/alpha_v0.0.rda")
load(file = "outputs/spnum_v0.0.rda")
alphTemp <- alpha
alphTemp[alphTemp > 999.99] <- 999.99
alphDF <- data.frame(alphTemp, 1:length(alphTemp))
alphDF$spnum <- spnum
alphDF$hectad <- unique(dat$hectad)
alphDF <- merge(alphDF, unicorn_results$stat, by = "hectad", by.y = "Location")
plot(alphDF$alphTemp, alphDF$Alpha)
cor(alphDF$alphTemp, alphDF$Alpha) # 0.99923324
plot(alphDF$spnum, alphDF$Spnum_out)
cor(alphDF$spnum, alphDF$Spnum_out) # 0.9991448
