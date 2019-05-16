localrules:
	# CNVnator_PoN_filtering,
	CNVnator_postprocess

rule CNVnator_calling:
	input:
		cram = "results/alignments/{sample}.cram",
		crai = "results/alignments/{sample}.crai",
		fsa_dir = config["chromosome_sequences_dir"], ## Add note in config.yaml about this
	output:
		root = temp("temp/CNVnator/{sample}/{sample}.root"),
		tsv  = temp("temp/CNVnator/{sample}/{sample}_cnv_calls_bin_100.tsv"),
		vcf  = temp("temp/CNVnator/{sample}/{sample}_cnv_calls_bin_100.vcf"),
	params:
		bin_size = 100,
		outdir   = "temp/CNVnator/{sample}"
	resources:
		mem_gb     = 120,
		walltime_h = 96
	threads: 1
	log: "logs/CNVnator_calling/{sample}.log"
	benchmark: "benchmarks/CNVnator/{sample}.txt"
	shell:
		"("
		"module load shared tools ngs root/6.06.06 yeppp/1.0.0 perl/5.24.0 cnvnator/0.3.3; "
		"mkdir -p {params.outdir}; "
		"touch {params.outdir}/.notempty; "
		"cnvnator -root {output.root} -tree {input.cram}; "
		"cnvnator -root {output.root} -his {params.bin_size} -d {input.fsa_dir}; "
		"cnvnator -root {output.root} -stat {params.bin_size}; "
		"cnvnator -root {output.root} -partition {params.bin_size}; "
		"cnvnator -root {output.root} -call {params.bin_size} > {output.tsv}; "
		"cnvnator2VCF.pl {output.tsv} > {output.vcf}; "
		") 2> {log}"

rule CNVnator_filtering:
	input:
		vcf = rules.CNVnator_calling.output.vcf,
		bed = config["bed"]["WES_strict"],
	output:
		vcf_pval = temp("temp/CNVnator/{sample}/{sample}_filtered_pvalue.vcf"),
		vcf_exon = temp("temp/CNVnator/{sample}/{sample}_filtered_pvalue_exons.vcf"),
	resources:
		mem_gb     = 20,
		walltime_h = 5
	threads: 1
	log: "logs/CNVnator_filtering/{sample}.log"
	shell:
		"("
		"module load shared tools ngs bedtools/2.27.1 vcflib/1.0.0-rc2; "
		"vcffilter -f 'natorP1 < 0.01 & natorP2 < 0.01' {input.vcf} > {output.vcf_pval}; "
		"bedtools intersect -a {output.vcf_pval} -b {input.bed} -wa -wb -header > {output.vcf_exon}; "
		") 2> {log}"

rule CNVnator_postprocess:
	input:
		vcf = rules.CNVnator_calling.output.vcf,
		# vcf = rules.CNVnator_filtering.output.vcf_pval
		# vcf = rules.CNVnator_PoN_filtering.output.vcf
	output:
		bed = "results/CNVnator/{sample}.bed"
	resources:
		mem_gb     = 4,
		walltime_h = 5
	threads: 1
	log: "logs/CNVnator_postprocess/{sample}.log"
	run:
		with open(input.vcf, "r") as i, open(output.bed, "w") as o:
			for line in i:
				line = line.rstrip()
				if line.startswith("#"):
					continue
				cols = line.split("\t")
				## CHR POS
				out = cols[0:2]
				## END
				END = cols[7].split(';', 1)[0]
				assert END.startswith("END=")
				END = END.split('=')[1]
				out.append(END)
				## CN
				CN = cols[-1].split(':', 1)[1]
				out.append(CN)
				pval = cols[7].split(';')[5]
				assert pval.startswith("natorP1=")
				pval = pval.split('=')[1]
				out.append(pval)
				print('\t'.join(out), file = o)

rule CNVnator:
	input:
		expand("results/{tool}/{sample}.bed", tool = "CNVnator", sample = SAMPLES_WGS)
