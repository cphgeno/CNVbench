localrules: ControlFREEC_postprocess
rule ControlFREEC_estimate_gender:
	input:
		bam = "results/alignments/{sample}.bam",
		bai = "results/alignments/{sample}.bai",
	output:
		chrX   = temp("temp/ControlFREEC/{sample}/mosdepth/chrX.mosdepth.global.dist.txt"),
		chrY   = temp("temp/ControlFREEC/{sample}/mosdepth/chrY.mosdepth.global.dist.txt"),
		gender = temp("temp/ControlFREEC/{sample}_gender.txt"),
	params:
		prefixX = "temp/ControlFREEC/{sample}/mosdepth/chrX",
		prefixY = "temp/ControlFREEC/{sample}/mosdepth/chrY",
	threads: 1
	resources:
		walltime_h = 5,
		mem_gb     = 4,
	log:    "logs/ControlFREEC_estimate_gender/{sample}.log"
	benchmark: "benchmarks/ControlFREEC/{sample}.txt"
	shell:
		"(module load tools htslib/1.9 mosdepth/0.2.4;\n"
		">&2 echo 'Estimate chrX coverage';\n"
		"mosdepth -n --fast-mode --chrom X {params.prefixX} {input.bam};\n"
		">&2 echo 'Estimate chrY coverage';\n"
		"mosdepth -n --fast-mode --chrom Y {params.prefixY} {input.bam};\n"
		"XCOV=$(grep 'total\t10\t' {output.chrX} | cut -f 3);\n"
		"YCOV=$(grep 'total\t10\t' {output.chrY} | cut -f 3);\n"
		"if (( $(echo \"$XCOV > 0.1 && ! $YCOV > 0.1\" | bc) )); then "
		"gender=XX; else gender=XY; fi; "
		"echo $gender > {output.gender} "
		") 2> {log}"

rule ControlFREEC_calling:
	input:
		bam    = "results/alignments/{sample}.bam",
		bai    = "results/alignments/{sample}.bai",
		gender = "temp/ControlFREEC/{sample}_gender.txt",
		chrdir = config["chromosome_sequences_dir"],
		chrlen = config["chromosome_lengths"],
	output:
		cnv   = temp("temp/ControlFREEC/{sample}/{sample}.bam_CNVs"),
		ratio = temp("temp/ControlFREEC/{sample}/{sample}.bam_ratio.txt"),
	params:
		outdir = "temp/ControlFREEC/{sample}",
	threads: 28
	resources:
		walltime_h =  24,
		mem_gb     = 120,
	log:    "logs/ControlFREEC_calling/{sample}.log"
	shell:
		"(module load tools ngs gcc/8.2.0 anaconda3/4.4.0 "
		"intel/perflibs/64/2019_update3 intel/redist/2019_update2 R/3.5.0 "
		"samtools/1.9 bedtools/2.27.1 sambamba/0.6.7 control-freec/11.5; "
		"tools/ControlFREEC/wrapper.py "
		"--input {input.bam} "
		"--sex $(cat {input.gender}) "
		"--outputDir {params.outdir} "
		"--maxThreads {threads} "
		"--sambamba $(which sambamba) "
		"--SambambaThreads 2 "
		"--coefficientOfVariation 0.05 "
		"--minCNAlength 1 "
		"--forceGCcontentNormalization 0 "
		"--ploidy 2 "
		"--chrFiles {input.chrdir} "
		"--chrLenFile {input.chrlen} "
		">&2 "
		") 2> {log}"

rule ControlFREEC_significance:
	input:
		cnv   = "temp/ControlFREEC/{sample}/{sample}.bam_CNVs",
		ratio = "temp/ControlFREEC/{sample}/{sample}.bam_ratio.txt",
	output:
		txt   = "temp/ControlFREEC/{sample}/{sample}.bam_CNVs.p.value.txt"
	threads: 1
	resources:
		walltime_h = 10,
		mem_gb     =  4,
	log: "logs/ControlFREEC_significance/{sample}.log"
	shell:
		"(module load tools ngs gcc/8.2.0 intel/perflibs/64/2019_update3 intel/redist/2019_update2 R/3.5.0; "
		# "module load module load rdxplorer/3.2; "
		"Rscript tools/ControlFREEC/assess_significance.R {input.cnv} {input.ratio} "
		") 2> {log}"

rule ControlFREEC_filter:
	input:
		"temp/ControlFREEC/{sample}/{sample}.bam_CNVs.p.value.txt",
	output:
		p1 = "temp/ControlFREEC/{sample}/{sample}_1.txt",
		p5 = "temp/ControlFREEC/{sample}/{sample}_5.txt",
	threads: 1
	log: "logs/ControlFREEC_filter/{sample}.log"
	shell:
		"perl -ane 'print if $F[5] <= 0.05' < {input} > {output.p5}; "
		"perl -ane 'print if $F[5] <= 0.01' < {input} > {output.p1}; "

rule ControlFREEC_postprocess:
	input:
		txt = "temp/ControlFREEC/{sample}/{sample}.bam_CNVs.p.value.txt"
	output:
		bed = "results/ControlFREEC/{sample}.bed"
	log: "logs/ControlFREEC_postprocess/{sample}.log"
	shell:
		"("
		"tail -n +2 {input.txt} | cut -f 1-4,6 > {output.bed} "
		") 2> {log}"

rule ControlFREEC:
	input:
		expand("results/{tool}/{sample}.bed", tool = "ControlFREEC", sample = SAMPLES_WGS)

