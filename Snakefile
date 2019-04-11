#!/usr/bin/env snakemake

import sys


SAMPLES, = glob_wildcards("reads/{sample}_R1.fastq.gz")

configfile: "config.yaml"

# def eprint(*args, **kwargs):
# 	print(*args, **kwargs, file = sys.stderr)

rule all:
	input: expand("results/{tool}/{sample}.bed", tool = config["TOOLS"], sample = SAMPLES)
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
for tool in config["TOOLS"]:
	include: f"tools/{tool}/rules.smk"
