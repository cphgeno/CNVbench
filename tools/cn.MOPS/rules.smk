def split(a, n):
	"""Split a in n chunks"""
	n = round(n)
	k, m = divmod(len(a), n)
	return (a[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(n))

## Get batches with approx. size of 10 samples
cnmops_wgs_batches = split(SAMPLES_WGS, len(SAMPLES_WGS) / 10)

rule cnmops:
	input:
		expand("results/cn.MOPS/{sample}.bed", sample = SAMPLES_WGS),
		expand("results/cn.MOPS/{sample}.bed", sample = SAMPLES_WES)

for batch_samples in cnmops_wgs_batches:
	rule:
		input:
			input_bams = expand("results/alignments/{sample}.bam", sample = batch_samples),
		output:
			cnv_beds = expand("results/cn.MOPS/{sample}.bed", sample = batch_samples)
		params:
			result_dir = "results/cn.MOPS/",
			bin_size = 500 # Assuming 150bp reads, 30x depth and 100 reads per bin, c.f. https://support.bioconductor.org/p/75969/
		log: "logs/cn.mops_wgs/samples_{}.log".format("-".join(batch_samples))
		benchmark: "benchmarks/cn.MOPS/{}.txt".format('_'.join(batch_samples))
		threads: 28
		resources:
			mem_gb = 120,
			walltime_h = 50
		shell:
			"(module load tools intel/redist/2019_update2 intel/perflibs/64 R/3.5.0; "
			"Rscript cnmops_wgs.R {params.result_dir} {params.bin_size} {threads} {input.input_bams}; "
			") 2> {log}"

rule cnmops_wes:
	input:
		input_bams = expand("results/alignments/{sample}.bam", sample = SAMPLES_WES),
	output:
		cnv_beds = expand("results/cn.MOPS/{sample}.bed", sample = SAMPLES_WES)
	params:
		result_dir = "results/cn.MOPS/",
		exon_bed = "/home/projects/cu_10047/data/references/hg19/regions/S04380110_Covered_BED.bed"
	log: "logs/cn.mops_wes/WES-samples.log"
	benchmark: "benchmarks/cn.MOPS/WES-samples.txt"
	threads: 28
	resources:
		mem_gb = 120,
		walltime_h = 50
	shell:
		"(module load tools intel/redist/2019_update2 intel/perflibs/64 R/3.5.0; "
		"Rscript cnmops_wes.R {params.result_dir} {params.exon_bed} {threads} {input.input_bams}; "
		") 2> {log}"

