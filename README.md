# fRescalo
FreScaLO (Frequency Scaling using Local Occupancy) was a method invented by Dr Mark O. Hill (of Hill number, DCA, and other fame) to analyse biological recording data (i.e. species occurrence data) collected with time- and space-varying effort.

R rewrite of Hill (2012; DOI:10.1111/j.2041-210X.2011.00146.x) fortran code available at https://www.brc.ac.uk/biblio/frescalo-computer-program-analyse-your-biological-records

This code was mainly written to aid understanding of the paper. Note that Dr Jon Yearsley has also provided what looks like a much faster R version of Frescalo here: https://github.com/DrJonYearsley/frescaloR

Also see the useful insights of R.J. Bijlsma (2013), including VBA re-writes of Hill's fortran code, here: https://www.blwg.nl/mossen/onderzoek/rapporten/BLWGRapport15.pdf

## Update March 2024
Added the calculation of time factor standard deviations in script `frescaloR_v1.R`

## Update June 2024
Minor change L66 script `frescaloR_v1.R` in line with Hill fortran code marginally increased correlation with that program's outputs (although it was already always more than 0.99 for any given output)
