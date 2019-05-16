#!/usr/bin/env Rscript

library(CODEX2)

load("preprocess_codex2.RData")

# Get args
args <- commandArgs(trailingOnly = TRUE)
result_dir <- args[1]
chr <- as.integer(args[2])

# Run CODEX2.
chr_index <- which(seqnames(ref_qc) == chr)
norm_obj_null <- normalize_null(Y_qc = Y_qc[chr_index,],
                               gc_qc = gc_qc[chr_index],
                               K = 1:5, 
                               N = N)
Yhat_null <- norm_obj_null$Yhat
AIC_null <- norm_obj_null$AIC 
BIC_null <- norm_obj_null$BIC
RSS_null <- norm_obj_null$RSS

finalcall_CBS <- segmentCBS(Y_qc[chr_index,],
                            Yhat_null, 
                            optK = which.max(BIC_null),
                            K = 1:10,
                            sampname_qc = sample_names,#paste('sample', 1:ncol(Y_qc), sep=""),
                            ref_qc = ranges(ref_qc)[chr_index],
                            chr = chr, 
                            lmax = 400, 
                            mode = "integer")

# Filtering
filter1 <- finalcall_CBS$length_kb <= 200
filter2 <- finalcall_CBS$length_kb / (finalcall_CBS$ed_exon - finalcall_CBS$st_exon + 1) < 50
finalcall_CBS_filter <- finalcall_CBS[filter1 & filter2, ]

filter3 <- finalcall_CBS_filter$lratio>40
filter4 <- (finalcall_CBS_filter$ed_exon - finalcall_CBS_filter$st_exon) > 1
finalcall_CBS_filter <- finalcall_CBS_filter[filter3|filter4,]

# Print bed
df_to_bed <- function(df) {
  bed <- data.frame(chrom = df$chr, chromStart =  df$st_bp - 1, chromEnd = df$ed_bp, copyNumber = df$copy_no, likelihood_ratio = df$lratio)
  rownames(bed) <- c()
  return(bed)
}

write_bed <- function(bed, file_name) {
  write.table(bed, file = file_name, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

write_sample_cnvs <- function(sample_name, cnv_results, filter, chr, dir) {
  sample_cnvs <- cnv_results[cnv_results$sample_name == sample_name,]
  bed <- df_to_bed(sample_cnvs)
  if (filter == TRUE) {
    file_name <- paste0(sample_name, "_chr", chr, "_filtered", ".bed")
  } else {
    file_name <- paste0(sample_name, "_chr", chr, ".bed")
  }
  path <- file.path(dir, file_name)
  write_bed(bed, path)
}

sapply(input_sample_names, write_sample_cnvs, cnv_results = finalcall_CBS, filter = FALSE, chr = chr, dir = result_dir)
sapply(input_sample_names, write_sample_cnvs, cnv_results = finalcall_CBS_filter, filter = TRUE, chr = chr, dir = result_dir)
