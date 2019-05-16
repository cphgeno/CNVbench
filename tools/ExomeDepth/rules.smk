localrules:
	ExomeDepth_convert

rule ExomeDepth:
	input:
		expand("results/{tool}/{sample}.bed", tool = "ExomeDepth", sample = SAMPLES_WES)

rule ExomeDepth_calling_WES:
	input:
		bam = "results/alignments/{sample}.bam",
		bai = "results/alignments/{sample}.bai",
		ref = config["human_reference"],
		bed = config["bed"]["WES_probes"],
		pon = config["pon"]["WES"],
	output:
		tsv = "results/ExomeDepth/{sample}.tsv"
	params:
		outdir = "results/ExomeDepth/",
	resources:
		mem_gb     = 20,
		walltime_h = 16
	threads: 1
	log: "logs/ExomeDepth_calling/{sample}.log"
	benchmark: "benchmarks/ExomeDepth/{sample}.txt"
	shell:
		"(module load shared tools ngs intel/compiler/64/2019_update3 gcc/8.2.0 R/3.5.0; "
		"Rscript --vanilla --slave tools/ExomeDepth/ExomeDepth_PoN.R --ref {input.ref} --bed {input.bed} --pon {input.pon} --out {params.outdir} {input.bam}; "
		") 2> {log}"

rule ExomeDepth_convert:
	input:
		tsv = rules.ExomeDepth_calling_WES.output.tsv
	output:
		bed = "results/ExomeDepth/{sample}.bed",
	resources:
		mem_gb     = 1,
		walltime_h = 1
	threads: 1
	log: "logs/ExomeDepth_convert/{sample}.log"
	shell:
		"(module load shared tools ngs anaconda3/4.4.0 perl/5.24.0; "
		"perl -MMath::Round -a -n -e 'chomp(); print($_.\"\t\".round($F[11]*2).\"\n\")' {input.tsv} | cut -f 5-7,9,15 | csvcut --tabs -c 3,1,2,5,4 | csvformat -T | tail -n +2 > {output.bed}; "
		") 2> {log}"
