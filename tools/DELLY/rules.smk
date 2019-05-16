localrules:
	DELLY_convert

rule DELLY:
	input:
		expand("results/{tool}/{sample}.bed", tool = "DELLY", sample = SAMPLES)

rule DELLY_calling:
	input:
		bam = "results/alignments/{sample}.bam"
	output:
		bcf_del = temp("temp/DELLY/{sample}_del.bcf"),
		bcf_ins = temp("temp/DELLY/{sample}_ins.bcf"),
		vcf     = "temp/DELLY/{sample}.vcf"
	params:
		ref = "/home/projects/cu_10047/data/references/hg19/sequence/ref.fa"
	resources:
		walltime_h = 120,
		mem_gb     = 80,
	log: "logs/DELLY_calling/{sample}.log"
	benchmark: "benchmarks/DELLY/{sample}.txt"
	shell:
		"(module load tools bcftools/1.9 delly2/0.7.5 ;\n"
		"delly call -g {params.ref} -t DEL -o {output.bcf_del} {input.bam} ;\n"
		"delly call -g {params.ref} -t DUP -o {output.bcf_ins} {input.bam} ;\n"
		"bcftools concat -a -O v -o {output.vcf} {output.bcf_del} {output.bcf_ins}"
		") 2> {log}"

rule DELLY_convert:
	input:
		vcf = rules.DELLY_calling.output.vcf
	output:
		bed = "results/DELLY/{sample}.bed"
	threads: 1
	log: "logs/vcf_to_bed/{sample}.log"
	run:
		BED = []

		with open(input.vcf, 'r') as vcf:
			for line in vcf:
				if line[0] == "#":
					continue
				line = line.rstrip().split("\t")
				if line[6] != "PASS":
					continue
				chrom = line[0]
				start = line[1]
				end = line[7].split(";")[4]
				number = line[9].split(":")[-5]
				GQ = line[-1].split(":")[2]
				BED.append([chrom,start,end,number,GQ])
		with open(output.new_bed,"w") as out:
			for entry in BED:
				line = "\t".join(entry)
				out.write(line+"\n")

