#!/bin/sh
### Note: No commands may be executed until after the #PBS lines
### Account information
#PBS -W group_list=cu_10047 -A cu_10047
### Job name (comment out the next line to get the name of the script used as the job name)
#PBS -N Control-FREEC
### Requesting time
#PBS -l walltime=08:00:00
### Requesting nodes and memory
#PBS -l nodes=1:ppn=28
### E-mail settings
#PBS -m abe

cd /home/projects/cu_10047/people/s134891/CNV-project

# Loading modules
module load anaconda3/4.4.0 
module load samtools/1.9
module load rdxplorer/3.2
module load R/3.5.0
module load bedtools/2.27.1
module load sambamba/0.6.7
module load control-freec/11.5


FILEDIR="/home/projects/cu_10047/data/analyzed_samples/project_hierarchy/research_projects/CNV_benchmark/results/alignments/"
SAMPLES=("NA12878" "01" "02" "03" "04" "05" "06" "07" "08")
GENDERS=("XX" "XX" "XX" "XY" "XX" "XX" "XX" "XY" "XX")

# Settings per run
TYPE="WGS"
num=1


# Get one sample like ${SAMPLES[0]}
BAMFILE="${FILEDIR}GB-$TYPE-${SAMPLES[$num]}.bam"
WORKDIR="/home/projects/cu_10047/people/s134891/CNV-project/GB-$TYPE-${SAMPLES[$num]}"
CONFIG="${WORKDIR}/CONFIG_GB-$TYPE-${SAMPLES[$num]}.txt"
gender="${GENDERS[$num]}"

# Make CONFIG file and output dir
mkdir -p ${WORKDIR}


echo -e "[general]

chrLenFile = /home/projects/cu_10047/data/analyzed_samples/project_hierarchy/research_projects/CNV_benchmark/tools/ControlFREEC/hg19.len
chrFiles = /home/projects/cu_10047/pipelines/gatk3/resources/hg19/chromosome_sequences/
ploidy = 2
forceGCcontentNormalization = 0
minCNAlength = 1
maxThreads = 28
outputDir = ${WORKDIR}
sambamba = /services/tools/sambamba/0.6.7/sambamba
SambambaThreads = 2
sex = ${gender}
coefficientOfVariation = 0.05


[sample]

mateFile = ${BAMFILE}
inputFormat = BAM
mateOrientation = FR" > ${CONFIG}



# Run actual command and time it
BEFORE=`date '+%s'`

freec -conf ${CONFIG}

AFTER=`date '+%s'`
TIME=$(($AFTER - $BEFORE))
echo "Processed GB-$TYPE-${SAMPLES[$num]}.bam"
echo "Doing this took $TIME seconds "

# Calculate significance of Control-FREEC predictions with an R script
cat /home/projects/cu_10047/data/analyzed_samples/project_hierarchy/research_projects/CNV_benchmark/tools/ControlFREEC/assess_significance.R | R --slave --args ${WORKDIR}/GB-$TYPE-${SAMPLES[$num]}.bam_CNVs ${WORKDIR}/GB-$TYPE-${SAMPLES[$num]}.bam_ratio.txt


echo "DONE."
