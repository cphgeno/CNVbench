localrules: CNVkit_convert

rule CNVkit:
	input:
		expand("results/{tool}/{sample}.bed", tool = "CNVkit", sample = SAMPLES)

rule CNVkit_calling_WES:
	input:
		bam = expand("results/alignments/{sample}.bam", sample = SAMPLES_WES),
		bai = expand("results/alignments/{sample}.bai", sample = SAMPLES_WES),
		ref = config["human_reference"],
		pon = glob.glob(config["pon"]["WES"] + "/*.bam"),
		bed = config["bed"]["WES_probes"],
		refFlat = config["refFlat"],
		mappable = config["bed"]["WGS_mappable"], #"resources/access-5k-mappable.hg19.bed",
	output:
		expand(["temp/CNVkit/{sample}.antitargetcoverage.cnn",
			"temp/CNVkit/{sample}.cnr",
			"temp/CNVkit/{sample}.cns",
			"temp/CNVkit/{sample}.targetcoverage.cnn",
			"temp/CNVkit/{sample}-diagram.pdf",
			"temp/CNVkit/{sample}-scatter.pdf"], sample = SAMPLES_WES),
		ref_cnn = "temp/CNVkit/WES.ref_cnn"
	params:
		outdir = "temp/CNVkit/"
	resources:
		mem_gb     = 50,
		walltime_h = 8
	threads: 8
	log: "logs/CNVkit_calling/GB-WES.log"
	benchmark: "benchmarks/CNVkit/WES_samples.txt"
	shell:
		"(module load shared tools ngs anaconda3/4.4.0; "
		"cnvkit.py batch {input.bam} --normal {input.pon} "
		"-p {threads} --fasta {input.ref} --annotate {input.refFlat} "
		"--access {input.mappable} --method hybrid --targets {input.bed} "
		"--output-reference {output.ref_cnn} "
		"--output-dir {params.outdir} --scatter --diagram || true; "
		") 2> {log}"

rule CNVkit_calling_WGS:
	input:
		bam = expand("results/alignments/{sample}.bam", sample = SAMPLES_WGS),
		bai = expand("results/alignments/{sample}.bai", sample = SAMPLES_WGS),
		ref = config["human_reference"],
		pon = glob.glob(config["pon"]["WGS"] + "/*.bam"),
		bed = config["bed"]["WGS_mappable"],
		refFlat = config["refFlat"],
		mappable = config["bed"]["WGS_mappable"], #"resources/access-5k-mappable.hg19.bed",
	output:
		expand(["temp/CNVkit/{sample}.antitargetcoverage.cnn",
			"temp/CNVkit/{sample}.cnr",
			"temp/CNVkit/{sample}.cns",
			"temp/CNVkit/{sample}.targetcoverage.cnn",
			"temp/CNVkit/{sample}-diagram.pdf",
			"temp/CNVkit/{sample}-scatter.pdf"], sample = SAMPLES_WGS),
		ref_cnn = "temp/CNVkit/WGS.ref_cnn"
	params:
		outdir = "temp/CNVkit/"
	resources:
		mem_gb     = 50,
		walltime_h = 96
	threads: 9
	log: "logs/CNVkit_calling/GB-WGS.log"
	benchmark: "benchmarks/CNVkit/WGS_samples.txt"
	shell:
		"(module load shared tools ngs anaconda3/4.4.0; "
		"cnvkit.py batch {input.bam} --normal {input.pon} "
		"-p {threads} --fasta {input.ref} --annotate {input.refFlat} "
		"--access {input.mappable} --method wgs --targets {input.bed} "
		"--output-reference {output.ref_cnn} "
		"--output-dir {params.outdir} --scatter --diagram || true; "
		") 2> {log}"

rule CNVkit_convert:
	input:
		cns = "temp/CNVkit/{sample}.cns",
		cnr = "temp/CNVkit/{sample}.cnr",
	output:
		cns = temp("temp/CNVkit/{sample}_tmp.cns"),
		bed = "results/CNVkit/{sample}.bed",
	resources:
		mem_gb     = 1,
		walltime_h = 1
	threads: 1
	log: "logs/CNVkit_convert/{sample}.log"
	shell:
		"(module load shared tools ngs anaconda3/4.4.0; "
		"cnvkit.py call {input.cns} -o {output.cns}; "
		"cnvkit.py segmetrics --sem -o /dev/stdout -s {output.cns} {input.cnr} "
		"| tail -n +2 | cut -f 1,2,3,6,10 > {output.bed} "
		# "cnvkit.py export bed {output.cns} --show all -y -o {output.bed}; "
		") 2> {log}"
