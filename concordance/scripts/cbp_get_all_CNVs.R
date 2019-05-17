#!/usr/bin/env Rscript
setwd('/home/projects/cu_10047/data/analyzed_samples/project_hierarchy/research_projects/CNV_benchmark/concordance/')

suppressPackageStartupMessages(library("rtracklayer", warn.conflicts = TRUE))
library("yaml")

RESULTSDIR <- "results_filtered"
TRUTHDIR   <- "combined_calls"
TOOLS     <- yaml::yaml.load_file("../tools.yaml")


## Check tool results ---------------------------------
callset <- data.frame()

tool_dirs <- list.dirs(RESULTSDIR, recursive=FALSE, full.names = FALSE)
tool_dirs <- subset(tool_dirs, tool_dirs %in% c(TOOLS$WES, TOOLS$WGS))

include_samples <- c(paste0('GB-WGS-', sprintf("%02d", 1:38)), paste0('GB-WES-', sprintf("%02d", 1:8)))

# Get results for each sample across all tools
for (sample in include_samples) {
  message("NOW: ", sample)
  
  # Read combined set
  truth_file <- paste(TRUTHDIR, paste0(sample, '_merged.bed'), sep = "/")
  truth    <- import(truth_file, genome = "GRCh37")
  
  # Make a df ready
  sample_table <- data.frame(row.names = paste(sample, as.character(truth), truth$name, sep = "_"))
  
  for (tool in tool_dirs) {
    
    # Find out if tool was run for sample
    filename <- paste0(RESULTSDIR, '/', tool, '/', sample, '.bed')
    if (file.exists(filename)) {
      bed      <- import(filename, genome = "GRCh37", extraCols = c(pval = "numeric"))
      
      # Find overlaps
      hits     <- findOverlaps(bed[bed$score!='-Inf'], truth)
      
      # Check that CNV type is the same for call set and truth (truth$name == "CNV" is for NA12878)
      hits     <- hits[bed$name[queryHits(hits)] == truth$name[subjectHits(hits)]]
      
      # Get called CNVs
      calls <- rep(0, length(truth))
      calls[subjectHits(hits)] <- 1
      
      sample_table <- cbind.data.frame(sample_table, calls)
      
    } else {
      # If tool was not run, add NAs
      sample_table <- cbind.data.frame(sample_table, rep(NA, length(truth)))
    }
  }
  colnames(sample_table) <- tool_dirs
  
  # Check for errors
  if (any(rowSums(sample_table)==0, na.rm = T)) {
    warning('Error in calling - please check file.')
  }
  
  callset <- rbind.data.frame(callset, sample_table)
}

write.table(callset, file = "all_called.txt", sep = "\t", quote = FALSE, row.names = TRUE)
