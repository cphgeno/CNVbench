#!/usr/bin/env snakemake

import sys


SAMPLES, = glob_wildcards("reads/{sample}_R1.fastq.gz")

configfile: "config.yaml"

# def eprint(*args, **kwargs):
# 	print(*args, **kwargs, file = sys.stderr)

rule all:
	input:
		WES = expand("results/{tool}/{sample}.bed", tool = config["TOOLS"]["WES"], sample = filter(lambda x: x.startswith("GB-WES-"), SAMPLES)),
		WGS = expand("results/{tool}/{sample}.bed", tool = config["TOOLS"]["WGS"], sample = filter(lambda x: x.startswith("GB-WGS-"), SAMPLES)),
	message: "Pipeline complete"

rule align:
	input: expand("results/alignments/{sample}.{ext}", sample = SAMPLES, ext = ["cram", "crai"])
	message: "Alignments generated"

# rule stats:
# 	input:  expand(directory("stats/{tool}/{sample}/"), tool = config["TOOLS"], sample = SAMPLES)
# 	output: "multiqc_report.html"
# 	log:    "logs/multiqc.log"
# 	shell:
# 		"(module unload anaconda3; module load anaconda2/4.4.0; "
# 		"multiqc --no-data-dir {input} --name {output}"
# 		") 2> {log}"

include: "tools/align.smk"
TOOLS = set(tool for library in config["TOOLS"] for tool in config["TOOLS"][library])
for tool in TOOLS:
	include: f"tools/{tool}/rules.smk"
