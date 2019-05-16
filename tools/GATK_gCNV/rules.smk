localrules: GATK_gCNV_convert

gCNV_PoN = {
	"WES": "resources/hg19/PoN/gCNV/20190414_WES",
	"WGS": "resources/hg19/PoN/gCNV/20190413_WGS_whole"
}
def get_gCNV_PoN(wildcards, pon_dict = gCNV_PoN):
	library = "WGS" if "WGS" in wildcards.sample else "WES"
	return pon_dict[library]

gCNV_shards = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22", "X", "Y"]

rule GATK_gCNV:
	input:
		expand("results/{tool}/{sample}.bed", tool = "GATK_gCNV", sample = SAMPLES)

rule GATK_gCNV_CollectReadCounts:
	input:
		cram = "results/alignments/{sample}.cram",
		ref = config["human_reference"],
		bed = lambda wildcards: "tools/GermlineCNVCaller/Intervals/WGS.filtered.interval_list" if "WGS" in wildcards.sample else config["bed"]["WES_probes"],
		PoN = get_gCNV_PoN,
	output:
		hdf5 = temp("temp/GATK_gCNV/{sample}.hdf5"),
	params:
		prefix = "{sample}.contig_ploidy",
		outdir = "temp/GATK_gCNV/contig_ploidy_calls"
	log: "logs/GATK_gCNV_CollectReadCounts/{sample}.log"
	benchmark: "benchmarks/GATK_gCNV_CollectReadCounts/{sample}.txt"
	threads: 1
	version: "4.1.0.0"
	shell:
		"(module load tools anaconda3/4.4.0 gcc/8.2.0 java/1.8.0 gatk/{version}; \n"
		"gatk CollectReadCounts "
		"-I {input.cram} "
		"-R {input.ref} "
		"-L {input.bed} "
		"--interval-merging-rule OVERLAPPING_ONLY "
		"-O {output.hdf5}; \n"
		"set +u; source activate {config[gatk_conda_env]}; set -u; \n"
		"mkdir -p {params.outdir}; "
		"gatk DetermineGermlineContigPloidy "
		"-I {output.hdf5} "
		"--model {input.PoN}/contig_ploidy-model "
		"--output {params.outdir} "
		"--output-prefix {params.prefix}; "
		") 2> {log}"


rule GATK_gCNV_calling:
	input:
		hdf5 = rules.GATK_gCNV_CollectReadCounts.output.hdf5,
		shards= "resources/GATK_gCNV/beds/{shard_number}.bed",
		PoN = get_gCNV_PoN,
	output:
		tracking = directory("temp/GATK_gCNV/shards/{sample}_shard_{shard_number}-tracking"),
		#model = directory("temp/GATK_gCNV/shards/{sample}_shard_{shard_number}-model"),
		calls = directory("temp/GATK_gCNV/shards/{sample}_shard_{shard_number}-calls"),
	params:
		prefix = "{sample}_shard_{shard_number}",
		outdir = "temp/GATK_gCNV/shards"
	version: "4.1.0.0"
	log:  "logs/GATK_gCNV_calling/{sample}_{shard_number}.log"
	benchmark: "benchmarks/GATK_gCNV_calling/{sample}_{shard_number}.txt"
	resources:
		mem_gb     = 50,
		walltime_h = 48
	shell:
		"(module load tools anaconda3/4.4.0 java/1.8.0 gatk/{version}; \n"
		"set +u; source activate {config[gatk_conda_env]}; set -u;"
		"gatk GermlineCNVCaller "
		"-I {input.hdf5} "
		"--run-mode CASE "
		"--model {input.PoN}/shard_{wildcards.shard_number}-model "
		"--contig-ploidy-calls temp/GATK_gCNV/contig_ploidy_calls/{wildcards.sample}.contig_ploidy-calls "
		"--interval-merging-rule OVERLAPPING_ONLY "
		"--output {params.outdir} "
		"--output-prefix {params.prefix}; "
		") 2> {log}"

#give argument file instead
rule GATK_gCNV_PostprocessGermlineCNVCalls:
	input:
		tracking = expand("temp/GATK_gCNV/shards/{{sample}}_shard_{shard_number}-tracking", shard_number=gCNV_shards),
		#model = expand("temp/GATK_gCNV/shards/{{sample}}_shard_{shard_number}-model", shard_number=gCNV_shards),
		calls = expand("temp/GATK_gCNV/shards/{{sample}}_shard_{shard_number}-calls", shard_number=gCNV_shards),
		dic   = config["sequence_dict"],
	output:
		vcf = temp("temp/GATK_gCNV/{sample}_GATK_gCNV_unsort.vcf.gz"),
		vcf_sort = "temp/GATK_gCNV/{sample}_GATK_gCNV.vcf.gz",
		tbi = "temp/GATK_gCNV/{sample}_GATK_gCNV.vcf.gz.tbi",
	params:
		model_shards = lambda wildcards: expand("--model-shard-path {PoN}/shard_{shard}-model ", PoN = get_gCNV_PoN(wildcards), shard = gCNV_shards),
		calls_shards = lambda wildcards: expand("--calls-shard-path temp/GATK_gCNV/shards/{sample}_shard_{shard}-calls ", sample = wildcards.sample, shard = gCNV_shards),
	version: "4.1.0.0"
	log:  "logs/GATK_gCNV_PostprocessGermlineCNVCalls/{sample}.log"
	benchmark: "benchmarks/GATK_gCNV_PostprocessGermlineCNVCalls/{sample}.txt"
	shell:
		"module load tools anaconda3/4.4.0 java/1.8.0 gatk/{version}; \n"
		"source activate {config[gatk_conda_env]}; \n"
		"gatk PostprocessGermlineCNVCalls "
		"{params.model_shards} "
		"{params.calls_shards} "
		"--contig-ploidy-calls temp/GATK_gCNV/shards/{wildcards.sample}.contig_ploidy-calls "
		"--autosomal-ref-copy-number 2 --allosomal-contig X --allosomal-contig Y "
		"--sample-index 0 --output-genotyped-intervals /dev/stdout "
		"| perl -p -e 's    /<DEL>,<DUP>/<CNV>/' | bgzip > {output.vcf}; \n"
		"conda deactivate; \n"
		"gatk SortVcf -I {output.vcf} --SEQUENCE_DICTIONARY {input.dic} --COMPRESSION_LEVEL 9 --CREATE_INDEX true -O {output.vcf_sort} "


rule GATK_gCNV_convert:
	input:
		vcf = rules.GATK_gCNV_PostprocessGermlineCNVCalls.output.vcf_sort
	output:
		bed = "results/GATK_gCNV/{sample}.bed"
	run:
		BED = []
		import gzip
		with gzip.open(input.vcf, 'rb') as f:
			for line in f:
				line = line.decode('ascii')
				if line.startswith("#"):
					continue
				l = line.split()
				chrom = l[0]
				start = int(l[1])
				end = l[7].strip("END=")
				status = l[-1].split(":")[0]
				if status == "0":
					continue
				copy_number = l[-1].split(":")[1]
				cnq = l[-1].split(":")[-1]
				BED.append([chrom,str(start),end,copy_number,cnq,"\n"])
			with open(output.bed,"w") as out:
				for entry in BED:
					line = "\t".join(entry)
					out.write(line)

