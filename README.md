# CNV tool benchmarking pipeline
This pipeline was set up to test a variety of CNV calling tools on WES and WGS samples sequenced using Illumina short-read sequencing technology.

# Installation and configuration
This pipeline is set up as a Snakemake workflow. More information on how to install and run Snakemake can be found at the [wiki](https://snakemake.readthedocs.io/en/stable/).

After installation of Snakemake, clone this repository and fill in the config.yaml.

This pipeline is set up to work with the [Environment Modules](http://modules.sourceforge.net/). It likely needs customisation or removal of the module statements in the rules (`tools/*/rules.smk`) for compatilibility with a given computational environment.

Some of the tools in this pipeline requires a panel of normals, which the pipeline does not generate automatically. One such tool is [GATK GermlineCNVCaller](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_copynumber_GermlineCNVCaller.php). More info on generating a panel of normals for this tool can be found at the [GATKs tool documentation](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/)

# Running the pipeline
Once the pipeline has been configured, a given tool can be run using the command `snakemake $tool`.

A list of available Snakemake targets can be generated by running `snakemake --list-target-rules`.

By default, Snakemake runs all the tools specified in the config.yaml by just running: `snakemake`.

Many options are available for customizing snakemake behaviour. These are all documented at the [Snakemake Wiki](https://snakemake.readthedocs.io/en/stable/).

# Sample data from project
The alignment file GB-WGS-NA12878.bam used in this project is available at the [NCBI Sequencing Read Archive (SRA)](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA543552).
