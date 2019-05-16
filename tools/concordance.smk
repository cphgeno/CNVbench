localrules: concordance_summarize

rule concordance_filter:
	input:
		rules.cnvs.input
	output:
		"concordance/callset.txt"
	log: "logs/concordance_filter.log"
	shell:
		"(module load intel/redist/2019_update2 intel/compiler/64/2019_update2 R/3.5.0; "
		"Rscript concordance/scripts/01_filter_results.R"
		") 2> {log}"

rule concordance_summarize:
	input:
		rules.concordance_filter.output
	output:
		"concordance/callset.txt"
	log: "logs/concordance_summarize.log"
	shell:
		"(module load intel/redist/2019_update2 intel/compiler/64/2019_update2 R/3.5.0; "
		"Rscript concordance/scripts/02_summarize_results.R"
		") 2> {log}"
