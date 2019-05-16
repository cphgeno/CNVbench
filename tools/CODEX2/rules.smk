localrules: CODEX2_postprocess

rule CODEX2:
	input:
		beds = expand("results/CODEX2/{sample}.bed", sample = SAMPLES_WES),

rule CODEX2_preprocess:
	input:
		bam = expand("results/alignments/{sample}.bam", sample = SAMPLES_WES),
		bai = expand("results/alignments/{sample}.bai", sample = SAMPLES_WES),
	output:
		preprocess_data = temp("temp/CODEX2/preprocess_codex2.RData")
	params:
		input_bam_dir = "results/alignments",
		pon_bam_dir   = config["pon"]["WES"],
		bed_file = config["bed"]["WES_probes"],
		outdir   = "temp/CODEX2/"
	log: "logs/codex2/CODEX2_preprocess.log"
	benchmark: "benchmarks/CODEX2_preprocess/WES_samples.txt"
	threads: 28
	resources:
		mem_gb = 120,
		walltime_h = 30
	shell:
		"(module load tools intel/redist/2019_update2 intel/perflibs/64 R/3.5.0; "
		"Rscript tools/CODEX2/install_codex2.R; "
		"Rscript tools/CODEX2/preprocess_codex2.R {params.input_bam_dir} {params.pon_bam_dir} {params.bed_file} {params.outdir}"
		";) 2> {log}"


rule CODEX2_calling:
	input:
		preprocess_data = rules.CODEX2_preprocess.output.preprocess_data
	output:
		chr_beds = temp(expand("temp/CODEX2/scatter_gather/{sample}_chr{{chr}}.bed", sample = SAMPLES_WES)),
		filtered_chr_beds = temp(expand("temp/CODEX2/scatter_gather/{sample}_chr{{chr}}_filtered.bed", sample = SAMPLES_WES))
	params:
		outdir = "temp/CODEX2/scatter_gather/"
	log: "logs/codex2/run_codex2_{chr}.log"
	benchmark: "benchmarks/CODEX2_calling/WES_samples_{chr}.txt"
	threads: 1
	resources:
		mem_gb = 5,
		walltime_h = 10
	shell:
		"(module load tools intel/redist/2019_update2 intel/perflibs/64 R/3.5.0; "
		"Rscript tools/CODEX2/codex2.R {input} {params.outdir} {wildcards.chr}"
		";) 2> {log}"

rule CODEX2_postprocess:
	input:
		chr_beds = expand("temp/CODEX2/scatter_gather/{{sample}}_chr{chr}.bed", chr = list(range(1, 23))),
		filtered_chr_beds = expand("temp/CODEX2/scatter_gather/{{sample}}_chr{chr}_filtered.bed", chr = list(range(1, 23)))
	output:
		bed = "temp/CODEX2/{sample}.bed",
		filtered_bed = "results/CODEX2/{sample}.bed",
	wildcard_constraints:
		sample = r"[\w\-]+",
		chr    = r"\d+",
	log: "logs/codex2/cat_chr_beds_{sample}.log"
	threads: 1
	shell:
		"(cat {input.chr_beds} > {output.bed}; "
		"cat {input.filtered_chr_beds} > {output.filtered_bed} "
		";) 2> {log}"

