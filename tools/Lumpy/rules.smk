localrules: Lumpy_convert

rule Lumpy:
	input:
		expand("results/{tool}/{sample}.bed", tool = "Lumpy", sample = SAMPLES)

rule Lumpy_get_discordants:
	input:
		bam = "results/alignments/{sample}.bam"
	output:
		discordant_bam = "temp/Lumpy/discordants/{sample}-discordants.bam"
	log: "logs/Lumpy/get_discordants/{sample}.log"
	benchmark: "benchmarks/Lumpy/{sample}_get_discordants.txt"
	threads: 7
	resources:
		mem_gb = 30,
		walltime_h = 10
	shell:
		"(module load tools samtools/1.9; "
		"samtools view -b -F 1294 {input.bam} | samtools sort -@ {threads} -o {output.discordant_bam} - "
		") 2> {log}"

rule Lumpy_get_splitters:
	input:
		bam = "results/alignments/{sample}.bam"
	output:
		splitter_bam = "temp/Lumpy/splitters/{sample}-splitters.bam"
	params:
		script = "/services/tools/lumpy/20180605/scripts/extractSplitReads_BwaMem"
	log: "logs/Lumpy/get_splitters/{sample}.log"
	benchmark: "benchmarks/Lumpy/{sample}_get_splitters.txt"
	threads: 7
	resources:
		mem_gb = 30,
		walltime_h = 10
	shell:
		"(module load tools anaconda2/4.4.0 samtools/1.9; "
		"samtools view -h {input.bam} | python2.7 {params.script} -i stdin | samtools sort -@ {threads} -o {output.splitter_bam} - "
		") 2> {log}"

rule Lumpy_express:
	input:
		bam = "results/alignments/{sample}.bam",
		discordant_bam = rules.Lumpy_get_discordants.output.discordant_bam,
		splitter_bam = rules.Lumpy_get_splitters.output.splitter_bam
	output:
		vcf = "temp/Lumpy/vcf/{sample}.vcf"
	params:
		config = "tools/Lumpy/Lumpy.config"
	log: "logs/Lumpy/lumpyexpress/{sample}.log"
	benchmark: "benchmarks/Lumpy/{sample}_express.txt"
	threads: 7
	resources:
		mem_gb = 30,
		walltime_h = 10
	shell:
		"(module load tools samblaster/0.1.24 sambamba/0.6.7 samtools/1.2 anaconda2/2.2.0 lumpy/20180605; "
		"lumpyexpress -K {params.config} "
		"-B {input.bam} "
		"-S {input.splitter_bam} "
		"-D {input.discordant_bam} "
		"-o {output.vcf} "
		") 2> {log}"

rule Lumpy_convert:
	input:
		vcf = rules.Lumpy_express.output.vcf
	output:
		bed = "results/Lumpy/{sample}.bed"
	log: "logs/Lumpy/vcf_to_bed/{sample}.log"
	threads: 1
	shell:
		"(module load tools bcftools/1.9; "
		"bcftools query -f '%CHROM\\t%POS\\t%INFO/END\\t%ALT\\t%INFO/SU\\n' {input.vcf} | "
		"awk '$4==\"<DEL>\" || $4==\"<DUP>\" {{gsub(\"<DUP>\", \"3\", $4); gsub(\"<DEL>\", \"1\", $4); print $1\"\\t\"$2\"\\t\"$3\"\\t\"$4\"\\t\"$5}}' > {output.bed} "
		") 2> {log}"
