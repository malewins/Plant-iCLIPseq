# apply package to AtGRP7 iCLIP with 3 replicates
#library(tidyverse)
library(rtracklayer)

# AtGRP7 iCLIP binding sites ####
peakFile <- "05_called_peaks/AtGRP7-GFP.bed"
peaks <- 
  import(con = peakFile,
         format = "BED",
         extraCols=c("additionalScores" = "character"))

# filter lowest scores
quants = quantile(peaks$score, probs = seq(0,1, by = 0.05))
peaksFilter = peaks[peaks$score >= quants[2]]

# import iCLIP signals from bigwigs
bigwig_folder <- "04_mapped/01_deduped_BAMs/beds/01_bigwigs/" 
clipFilesP <- 
  list.files(
    bigwig_folder, 
    pattern = "[1-3].\\+.bw$", 
    full.names = T)

clipFilesM <- 
  list.files(
    bigwig_folder, 
    pattern = "[1-3].-.bw$", 
    full.names = T)

# binding site definition 
library(BindingSiteFinder)
meta = data.frame(
  id = c(1:3),
  condition = factor(rep("AtGRP7-GFP", 3)), 
  clPlus = clipFilesP, clMinus = clipFilesM)
bds <- 
  BSFDataSetFromBigWig(
    ranges = peaksFilter, 
    meta = meta, 
    silent = F)
bds_sites <- 
  makeBindingSites(
    object = bds, 
    bsSize = 5, 
    minWidth = 1,
    minCrosslinks = 2, 
    minClSites = 2)

# apply reproducibility filter
bdsFinal <-
  reproducibilityFilter(
    bds_sites, 
    cutoff = 0.3, 
    n.reps = 2, 
    min.crosslinks = 2)

# prepare binding sites for export and save to BED 
bdsFinal = annotateWithScore(bdsFinal, getRanges(bds))
finalRanges = getRanges(bdsFinal)
names(finalRanges) = paste0("BS", seq_along(finalRanges))
bdsFinal = setRanges(bdsFinal, finalRanges)
bindingSites = getRanges(bdsFinal)
export(object = bindingSites, con = "06_binding_sites/AtGRP7-GFP_repbsites.bed", format = "BED")
