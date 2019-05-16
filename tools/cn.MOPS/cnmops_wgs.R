#!/usr/bin/env Rscript

library(cn.mops)
library(magrittr)

options(scipen = 999)

# Get arguments.
args <- commandArgs(trailingOnly = TRUE)
result_dir <- args[1]
bin_size <- as.integer(args[2])
threads <- as.integer(args[3])
bam_files <- tail(args, -3)

# Drop MT (fails otherwise) and get read counts in bin_size windows. Windows should contain 50-100 reads each.
seq_names <- c(as.character(seq(22)), "X", "Y")
bam_data_ranges <- getReadCountsFromBAM(bam_files, refSeqNames = seq_names, WL = bin_size, parallel = threads)

# Call CNVS and calculate integer copy numbers.
results <- cn.mops(bam_data_ranges, parallel = threads) %>% calcIntegerCopyNumbers()
cnvs <- cnvs(results)

# Function to format cn.mops result GRanges to bed-like data.frame.
granges_to_bed <- function(gr) {
    bed <- data.frame(chrom = as.character(seqnames(gr)), chromStart = start(ranges(gr)), chromEnd = end(ranges(gr)), copyNumber = sub("^CN", "", gr$CN), median = gr$median)
    return(bed)
}

# Function to write bed-like data.frame to bed file.
write_bed <- function(bed, file_name) {
    write.table(bed, file = file_name, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

# Function to write a particular sample in cn.mops result GRanges to bed file.
write_sample_cnvs <- function(sample_name, cnv_results) {
    sample_cnvs <- cnv_results[cnv_results$sampleName == sample_name]
    bed <- granges_to_bed(sample_cnvs)
    file_name <- sub(".bam$", ".bed", sample_name)
    path <- paste0(result_dir, file_name)
    write_bed(bed, path)
}

# Output CNV bed files for all input samples.
sapply(basename(bam_files), write_sample_cnvs, cnv_results = cnvs)
# save.image(file = paste0(result_dir, "cnmops_wgs.RData"))
