localrules:
	Manta_convert

rule Manta:
	input:
		expand("results/{tool}/{sample}.bed", tool = "Manta", sample = SAMPLES)

rule Manta_calling:
	input:
		bam = "results/alignments/{sample}.bam",
		bai = "results/alignments/{sample}.bai",
		ref = config["human_reference"],
	output:
		vcf_diploid = "temp/Manta/{sample}/results/variants/diploidSV.vcf.gz",
		tbi_diploid = "temp/Manta/{sample}/results/variants/diploidSV.vcf.gz.tbi",
		vcf_small = "temp/Manta/{sample}/results/variants/candidateSmallIndels.vcf.gz",
		tbi_small = "temp/Manta/{sample}/results/variants/candidateSmallIndels.vcf.gz.tbi",
		vcf_cand = "temp/Manta/{sample}/results/variants/candidateSV.vcf.gz",
		tbi_cand = "temp/Manta/{sample}/results/variants/candidateSV.vcf.gz.tbi",
		stats_summary = "temp/Manta/{sample}/results/stats/alignmentStatsSummary.txt",
		stats_cand_tsv = "temp/Manta/{sample}/results/stats/svCandidateGenerationStats.tsv",
		stats_cand_xml = "temp/Manta/{sample}/results/stats/svCandidateGenerationStats.xml",
		stats_graph = "temp/Manta/{sample}/results/stats/svLocusGraphStats.tsv",
	params:
		outdir = "temp/Manta/{sample}",
		wes    = lambda wildcards: "--exome" if "WES" in str(wildcards.sample) else ""
	resources:
		mem_gb     = 120,
		walltime_h = 24
	threads: 28
	log: "logs/Manta_calling/{sample}.log"
	benchmark: "benchmarks/Manta/{sample}.txt"
	shell:
		"(module load shared tools ngs anaconda2/4.4.0 manta/1.5.0; "
		"rm -rf {params.outdir}; "
		# Configure Manta
		"configManta.py "
		"--bam {input.bam} "
		"{params.wes} "
		"--referenceFasta {input.ref} "
		"--runDir {params.outdir}; "
		# Run Manta
		"{params.outdir}/runWorkflow.py --mode local --jobs {threads} --memGb {resources.mem_gb} || true; "
		# Create output folders
		"mkdir -p {params.outdir}/results/variants {params.outdir}/workspace/alignmentStats.xml.tmpdir; "
		# Re-run Manta
		"{params.outdir}/runWorkflow.py --mode local --jobs {threads} --memGb {resources.mem_gb}; "
		") 2> {log}"

rule Manta_convert:
	input:
		vcf = rules.Manta_calling.output.vcf_diploid
	output:
		bed = "results/Manta/{sample}.bed",
	resources:
		mem_gb     = 1,
		walltime_h = 1
	threads: 1
	log: "logs/Manta_convert/{sample}.log"
	shell:
		"(module purge; module load shared tools ngs anaconda2/4.4.0 bcftools/1.9; "
		# Restrict output to DEL, INS, and DUP
		"bcftools filter -i 'FILTER=\"PASS\" && (INFO/SVTYPE=\"DEL\" || INFO/SVTYPE=\"INS\" || INFO/SVTYPE=\"DUP\")' {input.vcf} "
		# Convert to BEDPE
		"| svtools vcftobedpe "
		# Remove comments
		"| fgrep -v '##' "
		# Extract columns
		"| csvcut -t -c '#CHROM_A,START_A,END_B,TYPE,22,QUAL' | csvformat -T | perl -p -e 's/:[^\t]*//' "
		# Format CN column
		"| perl -an -e '$cn=3; $cn=0 if($F[3] eq \"DEL\" && $F[4] eq \"1/1\"); $cn=1 if($F[3] eq \"DEL\" && $F[4] eq \"0/1\"); splice(@F,3,2,$cn); print(join(\"\t\",@F).\"\n\")' "
		# Save output
		"| fgrep -v '#' > {output.bed}; "
		") 2> {log}"
