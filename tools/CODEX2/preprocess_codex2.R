#!/usr/bin/env Rscript

library(CODEX2)

# Get args
args <- commandArgs(trailingOnly = TRUE)
input_bam_dir <- args[1]
pon_bam_dir   <- args[2]
bed_file <- args[3]
outdir   <- args[4]

# Get bam files from input and PoN directories and extract sample names. Use 50 PoN samples.
input_bam_files <- list.files(input_bam_dir, pattern = "GB-WES-.*.bam$")
input_bam_paths <- gsub("//", "/", file.path(input_bam_dir, input_bam_files))
input_sample_names <- sub("[.]bam", "", input_bam_files)
pon_bam_files <- head(list.files(pon_bam_dir, pattern = "RHS-.*-Blood-WES_.*.bam$"), 50)
pon_bam_paths <- gsub("//", "/", file.path(pon_bam_dir, pon_bam_files))
pon_sample_names <- sub("[.]bam", "", pon_bam_files)
bam_files <- c(input_bam_files, pon_bam_files)
bam_paths <- c(input_bam_paths, pon_bam_paths)
sample_names <- c(input_sample_names, pon_sample_names)

# Pre-processing
project_name <- "CODEX2"
bam_bed_obj <- getbambed(bamdir = bam_paths,
                         bedFile = bed_file,
                         sampname = sample_names,
                         projectname = project_name)
ref <- bam_bed_obj$ref

# GC content and mappability
genome <- BSgenome.Hsapiens.UCSC.hg19
gc <- getgc(ref, genome = genome)
mapp <- getmapp(ref, genome = genome)
values(ref) <- cbind(values(ref), DataFrame(gc, mapp))

# Raw read depth
coverage_obj <- getcoverage(bam_bed_obj, mapqthres = 20)
Y <- coverage_obj$Y

# Quality control
qc_obj <- qc(Y, sample_names, ref,
             cov_thresh = c(20, 4000),
             length_thresh = c(20, 2000),
             mapp_thresh = 0.9,
             gc_thresh = c(20, 80))
Y_qc <- qc_obj$Y_qc
sample_name_qc <- qc_obj$sampname_qc
ref_qc <- qc_obj$ref_qc
qcmat <- qc_obj$qcmat
gc_qc <- ref_qc$gc

# Running CODEX2
Y_nonzero <- Y_qc[apply(Y_qc, 1, function(x){!any(x == 0)}),]
pseudo_sample <- apply(Y_nonzero, 1, function(x){exp(1 / length(x) * sum(log(x)))})
N <- apply(apply(Y_nonzero, 2, function(x){x/pseudo_sample}), 2, median)

# Save.
save.image(file = paste(outdir, "preprocess_codex2.RData", sep = "/"))
