#!/usr/bin/env snakemake

import sys
import glob

SAMPLES, = glob_wildcards("reads/{sample}_R1.fastq.gz")
SAMPLES_WES = list(filter(lambda x: x.startswith("GB-WES-"), SAMPLES))
SAMPLES_WGS = list(filter(lambda x: x.startswith("GB-WGS-"), SAMPLES))

configfile: "config.yaml"

rule all:
	input:
		# align       = expand("results/alignments/{sample}.{ext}", sample = SAMPLES, ext = ["cram", "crai", "bam", "bai"]),
		align         = expand("stats/mosdepth/{sample}.regions.bed.gz", sample = SAMPLES),
		benchmark_wes = expand("results/{tool}/{sample}.bed", tool = config["TOOLS"]["WES"], sample = ["GB-WES-NA12878"]),
		benchmark_wgs = expand("results/{tool}/{sample}.bed", tool = config["TOOLS"]["WGS"], sample = ["GB-WGS-NA12878"]),
		cnvs_wes      = expand("results/{tool}/{sample}.bed", tool = config["TOOLS"]["WES"], sample = SAMPLES_WES),
		cnvs_wgs      = expand("results/{tool}/{sample}.bed", tool = config["TOOLS"]["WGS"], sample = SAMPLES_WGS),
		concordance   = expand("concordance/{file}.txt", file = ["summary", "mlpa_data", "hist_data"]),
	message: "Pipeline complete"

rule align:
	input: rules.all.input.align
	message: "Alignments generated"

rule benchmark:
	input:
		wes = rules.all.input.benchmark_wes,
		wgs = rules.all.input.benchmark_wgs,
	message: "Benchmark complete"

rule cnvs:
	input:
		WES = rules.all.input.cnvs_wes,
		WGS = rules.all.input.cnvs_wgs
	message: "CNV calling complete"

rule concordance:
	input:
		rules.all.input.concordance
	message: "Concordance calculation complete"

include: "tools/align.smk"
include: "tools/concordance.smk"
TOOLS = set(tool for library in config["TOOLS"] for tool in config["TOOLS"][library])
for tool in TOOLS:
	include: f"tools/{tool}/rules.smk"
