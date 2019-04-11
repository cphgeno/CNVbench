def get_read_group(read):
	import gzip
	with gzip.open(read) as f:
		(instrument, run_number, flowcell_id, lane, tile, X_pos, Y_pos, read, is_filtered, control_number, sample_barcode) = f.readline().decode()[1:].replace(':',' ').split()
	return '_'.join([instrument, run_number, flowcell_id])

rule bwa_mem:
	input:
		ref   = config["human_reference"],
		reads = expand("reads/{{sample}}_{dir}.fastq.gz", dir = ["R1", "R2"])
	output:
		sam   = temp("temp/align/bwa_mem/{sample}.sam")
	log: "logs/bwa_mem/{sample}.log"
	params:
		rg = lambda wildcards, input: '"@RG\\tID:{id}\\tPL:ILLUMINA\\tPU:{id}\\tLB:{sample}\\tSM:{sample}\\tCN:RH-GM"'.format(id = get_read_group(input.reads[0]), sample = wildcards.sample)
	threads: 14
	resources:
		mem_gb     = 75,
		walltime_h = 15
	version: '0.7.15'
	shell:
		"(module load moab tools samtools/1.9 bwa/{version}; "
		# Run BWA-MEM
		"bwa mem -M -t {threads} -R {params.rg} {input.ref} {input.reads} > {output.sam}; "
		") 2> {log}"

## Note: fix_read_groups is excluded

rule mark_duplicates:
	input:
		sam = rules.bwa_mem.output.sam
	output:
		bam     = temp("temp/align/mark_duplicates/{sample}.bam"),
		bai     = temp("temp/align/mark_duplicates/{sample}.bai"),
		metrics = "stats/mark_duplicates/{sample}.dupmetrics.txt"
	log: "logs/mark_duplicates/{sample}.log"
	params:
		tmpfile = "temp/large_files/{sample}.mark_dup.tmp"
	threads: 7
	resources:
		mem_gb     = 30,
		walltime_h = 10
	version: "2.0.95-release-20190320141403"
	shell:
		"(module load moab tools perl/5.24.0 samtools/1.9 biobambam2/{version}; "
		# Mark Duplicates and Sort reads
		"cat {input.sam} | bamsormadup inputformat=sam threads={threads} tmpfile={params.tmpfile} SO=coordinate M={output.metrics} indexfilename={output.bai} > {output.bam}; "
		# Fix metrics file to be recognized by MultiQC
		"perl -p -i -e 's[^##METRICS$][# MarkDuplicates INPUT={output.bam}\n## METRICS CLASS\tDuplicationMetrics]' {output.metrics}; "
		") 2> {log}"

rule filter_reads:
	input:
		bam = rules.mark_duplicates.output.bam,
		bai = rules.mark_duplicates.output.bai,
		ref = config["human_reference"]
	output:
		bam = temp("temp/align/filter_reads/{sample}.bam"),
		bai = temp("temp/align/filter_reads/{sample}.bai")
	log: "logs/filter_reads/{sample}.log"
	version: "4.0.10.1"
	shell:
		"(module load shared moab tools java/1.8.0 gatk/{version}; "
		# Filter out bad reads (https://software.broadinstitute.org/gatk/documentation/tooldocs/current/)
		"gatk PrintReads -R {input.ref} -I {input.bam} "
		"-RF GoodCigarReadFilter "                          # Keep only reads containing good CIGAR string
		"-RF MappedReadFilter "                             # Filter out unmapped reads
		"-RF MappingQualityAvailableReadFilter "            # Filter out reads without available mapping quality
		"-RF MappingQualityNotZeroReadFilter "              # Filter out reads with mapping quality equal to zero
		"-RF MatchingBasesAndQualsReadFilter "              # Filter out reads where the bases and qualities do not match
		# "-RF MateDifferentStrandReadFilter "                # Keep only reads with mates mapped on the different strand
		# "-RF MateOnSameContigOrNoMappedMateReadFilter "     # Keep only reads whose mate maps to the same contig or is unmapped
		# "-RF NotDuplicateReadFilter "                       # Filter out reads marked as duplicate
		"-RF PassesVendorQualityCheckReadFilter "           # Filter out reads failing platfor/vendor quality checks
		"-RF ReadLengthEqualsCigarLengthReadFilter "        # Filter out reads where the read and CIGAR do not match in length
		"-RF WellformedReadFilter "                         # Keep only reads that are well-formed
		"-O {output.bam}; "
		") 2> {log} "

rule bqsr_create:
	input:
		bam      = rules.filter_reads.output.bam,
		bai      = rules.filter_reads.output.bai,
		ref      = config["human_reference"],
		dbsnp    = config["dbsnp"],
		indels   = config["mills"],
		snps1k   = config["snps1k"],
		indels1k = config["indels1k"]
	output:
		bqsr = "stats/bqsr/{sample}.table",
	log: "logs/bqsr/{sample}.create.log"
	threads: 1
	resources:
		mem_gb     = 10,
		walltime_h = 10
	version: "4.0.10.1"
	shell:
		"(module load moab tools java/1.8.0 samtools/1.9 gatk/{version}; "
		"gatk BaseRecalibrator "
		"--known-sites {input.dbsnp} "
		"--known-sites {input.indels} "
		"--known-sites {input.snps1k} "
		"--known-sites {input.indels1k} "
		"-R {input.ref} "
		"-I {input.bam} "
		"-O {output.bqsr} "
		") 2> {log} "

rule bqsr_apply:
	# https://drive.google.com/file/d/1CV4BGDzt_pD_kvJ-kUuOEQiGsVp4Ybn3/view
	# https://software.broadinstitute.org/gatk/documentation/article?id=7899#ApplyBQSR
	input:
		bam = rules.filter_reads.output.bam,
		bai = rules.filter_reads.output.bai,
		bqsr = rules.bqsr_create.output.bqsr,
		ref  = config["human_reference"]
	output:
		bam = temp("temp/align/bqsr/{sample}.bam"),
		bai = temp("temp/align/bqsr/{sample}.bai")
	log: "logs/bqsr/{sample}.apply.log",
	threads: 1
	resources:
		mem_gb     = 75,
		walltime_h = 15
	version: "4.0.10.1"
	shell:
		"(module load moab tools java/1.8.0 samtools/1.9 gatk/{version}; "
		"gatk ApplyBQSR -R {input.ref} "
		"-bqsr {input.bqsr} "
		"-I {input.bam} "
		"-O {output.bam}; "
		"samtools index -@ {threads} {output.bam} {output.bai}"
		") 2> {log}"

rule archive_cram:
	input:
		bam = rules.bqsr_apply.output.bam,
		bai = rules.bqsr_apply.output.bai,
		ref = config["human_reference"]
	output:
		cram = "results/alignments/{sample}.cram",
		crai = "results/alignments/{sample}.crai",
	log: "logs/archive_cram/{sample}.log"
	threads: 3
	resources:
		mem_gb     = 10,
		walltime_h = 5
	version: "1.9"
	shell:
		"(module load moab tools samtools/{version}; "
		"samtools view -C -T {input.ref} -@ {threads} {input.bam} > {output.cram}; "
		"samtools index -@ {threads} {output.cram} {output.crai}; "
		") 2> {log}"
