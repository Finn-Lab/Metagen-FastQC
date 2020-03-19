#!/bin/bash

usage()
{
cat << EOF
usage: $0 options

Clean and remove human DNA sequences from FASTQ files.

Needs bsub (-M 15000)

OPTIONS:
   -h      Show help message
   -t      Number of threads (recommended: 16)
   -f      Forward or single-end fastq file (.fastq or *_1.fastq.gz) [REQUIRED]
   -r	   Reverse fastq file (*_2.fastq or *_2.fastq.gz) [OPTIONAL]
EOF
}

THREADS=
FASTQ_R1=
FASTQ_R2=

while getopts “ht:r:f:” OPTION
do
     case $OPTION in
         h)
             usage
             exit 1
             ;;
         t)
             THREADS=${OPTARG}
             ;;
         f)
             FASTQ_R1=$OPTARG
             ;;
	 r)
             FASTQ_R2=$OPTARG
	     ;;
         ?)
             usage
             exit
             ;;
     esac
done

timestamp() {
  date +"%H:%M:%S"
}

echo "$(timestamp) [remove human sequences] Parsing command-line"
# check if all required arguments are supplied
if [[ -z ${FASTQ_R1} ]]
then
     echo "ERROR : Please supply required arguments"
     usage
     exit 1
fi

REF="/nfs/production/interpro/metagenomics/mags-scripts/bwa_hg38/hg38.fa"
if [ ${THREADS} -eq 1 ]
then
    THREADS_SAM=1
else
    THREADS_SAM=$((${THREADS}-1))
fi

if [[ ! -z ${FASTQ_R2} ]]
then
        echo "$(timestamp) [remove human sequences] Cleaning FASTQ files"
        trim_galore --paired ${FASTQ_R1} ${FASTQ_R2} -o $(dirname ${FASTQ_R1})
        name=${FASTQ_R1%%_1.fastq*}
        echo "$(timestamp) [remove human sequences] Mapping files to human genome (hg38)"
        bwa mem -M -t ${THREADS} ${REF} ${name}_1_val_1.fq.gz ${name}_2_val_2.fq.gz | samtools view -@ ${THREADS_SAM} -f 12 -F 256 -uS - -o ${name}_both_unmapped.bam 
	samtools sort -@ ${THREADS_SAM} -n ${name}_both_unmapped.bam -o ${name}_both_unmapped_sorted.bam
	bedtools bamtofastq -i ${name}_both_unmapped_sorted.bam -fq ${name}_clean_1.fastq -fq2 ${name}_clean_2.fastq
	echo "$(timestamp) [remove human sequences] Compressing output files"
	gzip ${name}_clean_1.fastq
	gzip ${name}_clean_2.fastq
        echo "$(timestamp) [remove human sequences] Cleaning tmp files"
        rm -rf ${name}_both_unmapped.bam ${name}_both_unmapped_sorted.bam ${name}_*_trimmed.fq.gz ${name}_*_val_*.fq.gz ${name}_*.fastq.gz_*txt
else
        echo "$(timestamp) [remove human sequences] Cleaning FASTQ files"
        trim_galore ${FASTQ_R1} -o $(dirname ${FASTQ_R1})
        name=${FASTQ_R1%%.fastq*}
        echo "$(timestamp) [remove human sequences] Mapping files to human genome (hg38)"
        bwa mem -M -t ${THREADS} ${REF} ${name}_trimmed.fq.gz | samtools view -@ ${THREADS_SAM} -f 4 -F 256 -uS - -o ${name}_unmapped.bam
        samtools sort -@ ${THREADS_SAM} -n ${name}_unmapped.bam -o ${name}_unmapped_sorted.bam
        bedtools bamtofastq -i ${name}_unmapped_sorted.bam -fq ${name}_clean.fastq
	echo "$(timestamp) [remove human sequences] Compressing output file"
	gzip ${name}_clean.fastq
        echo "$(timestamp) [remove human sequences] Cleaning tmp files"
        rm -rf ${name}_unmapped.bam ${name}_unmapped_sorted.bam ${name}_trimmed.fq.gz ${FASTQ_R1}_*txt
fi
