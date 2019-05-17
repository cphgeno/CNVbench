#!/usr/bin/env Rscript
get_script_dir <- function() {
	initial.options <- commandArgs(trailingOnly = FALSE)
	file.arg.name <- "--file="
	script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
	return(dirname(script.name))
}
setwd(paste(get_script_dir(), "..", sep ="/"))

suppressPackageStartupMessages(library("data.table"))
library("yaml")
library("lubridate", warn.conflicts = FALSE)

TRUTHDIR   <- "truth_set"
RESULTSDIR <- "results_filtered"
BENCHMARKDIR <- "../benchmarks"
TOOLS     <- yaml::yaml.load_file("../tools.yaml")

## Functions -----------------------------
flatten_mapply_results <- function(mapply_out.list) {
	mapply_out.df   <- lapply(mapply_out.list, as.data.frame)
	mapply_out.df$make.row.names = FALSE
	mapply_out      <- do.call("rbind", mapply_out.df)
	return(mapply_out)
}

extract_results <- function(sample, tool, resultsdir = RESULTSDIR, truthdir = TRUTHDIR) {
	res <- list(
		sample  = sample,
		tool    = tool,
		library = NA,
		TP      = NA,
		FP      = NA,
		FN      = NA,
		TP_FN   = NA,
		N_truth = NA,
		N_DEL   = NA,
		N_DUP   = NA
	)
	if (startsWith(sample, "GB-WGS")) {res$library = "WGS"} else {res$library = "WES"}
	if (tool == "Cytoscan") {
		bed_file <- sprintf("%s/%s.bed", truthdir, sample)
	} else {
		bed_file <- sprintf("%s/%s/%s.bed", resultsdir, tool, sample)
	}
	if (file.exists(bed_file)) {
		bed <- fread(bed_file, colClasses=list(character = 1))
		if (tool == "Cytoscan") {
			res$N_DEL <- sum(bed$V4 == "DEL")
			res$N_DUP <- sum(bed$V4 == "DUP")
		} else {
			res$N_DEL <- sum(bed$V4 == "DEL" & Inf > bed$V5 & bed$V5 >= 0)
			res$N_DUP <- sum(bed$V4 == "DUP" & Inf > bed$V5 & bed$V5 >= 0)
			truth_file <- sprintf("%s/%s.bed", truthdir, sample)
			if (file.exists(truth_file)) {
				# res$TP <- sum(0 < bed$V5 & bed$V5 <= 1) + sum(bed$V5[1 < bed$V5 & bed$V5 < Inf])
				# res$FN <- sum(bed$V5 == -1)
				# res$FP <- sum(bed$V5 ==  0)
				res$TP <- sum(bed$V5 ==  Inf)
				res$FN <- sum(bed$V5 == -Inf)
				res$FP <- with(res, N_DUP + N_DEL - TP)
				res$TP_FN <- with(res, TP + FN)
				res$N_truth <- R.utils::countLines(truth_file)[1]
			}
		}
	}
	return(res)
}

extract_hist_data <- function(sample, tool, breaks, resultsdir = RESULTSDIR, truthdir = TRUTHDIR) {
	res <- list(
		sample  = sample,
		tool    = tool,
		library = NA,
		bin     = NA,
		counts  = NA
	)
	if (startsWith(sample, "GB-WGS")) {res$library = "WGS"} else {res$library = "WES"}
	if (tool == "Cytoscan") {
		bed_file <- sprintf("%s/%s.bed", truthdir, sample)
	} else {
		bed_file <- sprintf("%s/%s/%s.bed", resultsdir, tool, sample)
	}
	if (file.exists(bed_file)) {
		bed <- fread(bed_file, colClasses=list(character = 1))
		if (tool != "Cytoscan") {
			# bed <- subset(bed, bed$V5 >= 0)
			bed <- subset(bed, Inf > bed$V5 & bed$V5 >= 0)
		}
		bed$width <- bed$V3 - bed$V2
		h <- hist(bed$width, breaks, plot = FALSE)
		res$bin    <- sprintf("%.0f-%.0f", head(breaks, -1), tail(breaks, -1))
		res$counts <- h$counts
	}
	return(res)
}

extract_mlpa_data <- function(sample, tool, resultsdir = RESULTSDIR, truthdir = TRUTHDIR) {
	res <- list(
		sample  = sample,
		tool    = tool,
		library = NA,
		overlap = NA
	)
	if (startsWith(sample, "GB-WGS")) {res$library = "WGS"} else {res$library = "WES"}
	bed_file <- sprintf("%s/%s/%s.bed", resultsdir, tool, sample)
	if (file.exists(bed_file)) {
		bed <- fread(bed_file, colClasses=list(character = 1))
		# if (any(bed$V5 == -1) || all(bed$V5) == 0) {
		# 	res$overlap = 0
		# } else {
		# 	res$overlap <- max(bed$V5)
		# }
		bed <- subset(bed, bed$V5 <= 1)
		res$overlap <- max(bed$V5, 0)
	}
	return(res)
}

extract_benchmarks <- function(sample, tool, benchdir = BENCHMARKDIR) {
	res <- list(
		sample  = sample,
		tool    = tool,
		library = NA,
		cputime = NA,
		peak_memory = NA
	)
	if (startsWith(sample, "GB-WGS")) {res$library = "WGS"} else {res$library = "WES"}
	if (tool == "CLC") {
		input <- sprintf("%s/%s/cpu_time.tsv", benchdir, tool)
		if (file.exists(input)) {
			df <- fread(input)
			df <- subset(df, df$sample == res$sample)
			res$cputime <- as.numeric(hm(df$`elapsed time`))
		}
	} else {
		if (tool %in% c("CNVkit", "CODEX2") || (tool == "cn.MOPS" && res$library == "WES")) {
			input <- Sys.glob(sprintf("%s/%s*/*%s*.txt", benchdir, tool, res$library))
		} else {
			input <- Sys.glob(sprintf("%s/%s*/*%s*.txt", benchdir, tool, sample))
		}
		if (length(input) > 0) {
			df <- do.call("rbind", lapply(input, fread))
			res$cputime     <- sum(df$s)
			res$peak_memory <- max(df$max_rss)
		}
	}
	return(res)
}

## Initialize ---------------------------
# tool_dirs   <- list.dirs(RESULTSDIR, recursive=FALSE, full.names = FALSE)
tool_dirs   <- union(TOOLS$WGS, TOOLS$WES)
cohorts <- fread("../cohorts.txt")
samples <- list(
	all = subset(cohorts, cohorts$cohort != "mlpa")$sample
)
samples$WES <- subset(samples$all, grepl("WES", samples$all))
samples$WGS <- subset(samples$all, grepl("WGS", samples$all))

## Summarize results -----------------------------
message("Creating summary.txt")
l <- ifelse(tool_dirs %in% TOOLS$WGS & tool_dirs %in% TOOLS$WES, list(c(samples$WGS, samples$WES)), ifelse(tool_dirs %in% TOOLS$WES, list(samples$WES), list(samples$WGS)))
map_args <- list(
	sample = do.call("c", l),
	tool   = unlist(mapply(rep, tool_dirs, sapply(l, length)), use.names = FALSE)
)

summary <- flatten_mapply_results(mapply(extract_results, map_args$sample, map_args$tool, SIMPLIFY = FALSE))
summary$precision <- with(summary, TP / (TP + FP))
summary$recall    <- with(summary, TP / (TP + FN))
fwrite(summary, file = "summary.txt", sep = "\t", na = ".", quote = FALSE)


## Make histogram data ---------------------------
message("Creating hist_data.txt")
hist_breaks <- c(0L, 500L, 1000L, 10000L, 100000L, Inf)
hist_data <- flatten_mapply_results(mapply(extract_hist_data, map_args$sample, map_args$tool, MoreArgs = list(breaks = hist_breaks), SIMPLIFY = FALSE))
fwrite(hist_data, file = "hist_data.txt", sep = "\t", na = ".", quote = FALSE)


## Summarize MLPA results ----------------------------
message("Creating mlpa_overlap.txt")
MLPATRUTHDIR <- "truth_set/mlpa"
mlpa_samples <- list.files(MLPATRUTHDIR)

samples$MLPA = subset(cohorts, cohorts$cohort == "mlpa")$sample
samples$MLPA_WES <- subset(samples$MLPA, grepl("WES", samples$MLPA))
samples$MLPA_WGS <- subset(samples$MLPA, grepl("WGS", samples$MLPA))

tool_dirs <- subset(tool_dirs, tool_dirs != "Cytoscan")
l <- ifelse(tool_dirs %in% TOOLS$WGS & tool_dirs %in% TOOLS$WES, list(c(samples$MLPA_WGS, samples$MLPA_WES)), ifelse(tool_dirs %in% TOOLS$WES, list(samples$MLPA_WES), list(samples$MLPA_WGS)))
mlpa_args <- list(
	sample = do.call("c", l),
	tool   = unlist(mapply(rep, tool_dirs, sapply(l, length)), use.names = FALSE)
)

mlpa_overlap    <- flatten_mapply_results(mapply(extract_mlpa_data, mlpa_args$sample, mlpa_args$tool, SIMPLIFY = FALSE))
fwrite(mlpa_overlap, file = "mlpa_overlap.txt", sep = "\t", na = ".", quote = FALSE)


## Summarize computational benchmark results --------------------------
message("Creating benchmarks.txt")
samples$REF = subset(cohorts, cohorts$cohort == "reference")$sample
samples$REF_WES <- subset(samples$REF, grepl("WES", samples$REF))
samples$REF_WGS <- subset(samples$REF, grepl("WGS", samples$REF))

benchmark_args <- list(
	sample = c(rep(samples$REF_WES, length(TOOLS$WES)), rep(samples$REF_WGS, length(TOOLS$WGS))),
	tool   = c(TOOLS$WES, TOOLS$WGS)
)
benchmarks <- flatten_mapply_results(mapply(extract_benchmarks, benchmark_args$sample, benchmark_args$tool, SIMPLIFY = FALSE))
fwrite(benchmarks, file = "benchmarks.txt", sep = "\t", na = ".", quote = FALSE)

