curDir<- getwd()
library(Hmisc)
getRs('convertL2M.r', put='source')   # gets from github
# setwd("/media/alex/Seagate Expansion Drive/Treballs/2011-Materials en elaboracio per assignatura ANALISI DE DADES D'ALT RENDIMENT/1.1. Introduccio Als microarrays-c1")
sourceFileName <- "Selecting_GRM_with_Lheuristic.Rnw"
destFileName <- "Selecting_GRM_with_Lheuristic.Rmd"
convertL2M(sourceFileName, destFileName)
setwd(curDir)
