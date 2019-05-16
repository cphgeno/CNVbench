#!/usr/bin/env Rscript
setwd('/home/projects/cu_10047/data/analyzed_samples/project_hierarchy/research_projects/CNV_benchmark/concordance/')

suppressPackageStartupMessages(library("rtracklayer", warn.conflicts = TRUE))
library("yaml")

RESULTSDIR <- "truth_set"
TRUTHDIR   <- "combined_calls"
TOOLS     <- yaml::yaml.load_file("../tools.yaml")


## Check tool results ---------------------------------
callset <- data.frame()

include_samples <- c(paste0('GB-WGS-', sprintf("%02d", 1:38)), paste0('GB-WES-', sprintf("%02d", 1:8)))

# Get results for each sample across all tools
for (sample in include_samples) {
  message("NOW: ", sample)
  
  # Read combined set
  truth_file <- paste(TRUTHDIR, paste0(sample, '_merged.bed'), sep = "/")
  truth    <- import(truth_file, genome = "GRCh37")

    
  # Find out if tool was run for sample
  filename <- paste0(RESULTSDIR, '/', sample, '.bed')
  bed      <- import(filename, genome = "GRCh37")
  
  # Find overlaps
  hits     <- findOverlaps(bed, truth)
  
  # Check that CNV type is the same for call set and truth (truth$name == "CNV" is for NA12878)
  hits     <- hits[bed$name[queryHits(hits)] == truth$name[subjectHits(hits)]]
  
  # Get called CNVs
  calls <- rep(0, length(truth))
  calls[subjectHits(hits)] <- 1

  
  callset <- rbind.data.frame(callset, data.frame(calls, row.names = paste(sample, as.character(truth), truth$name, sep = "_")))
}

write.table(callset, file = "cytoscan_called.txt", sep = "\t", quote = FALSE, row.names = TRUE)
