#!/bin/bash

usage()
{
cat << EOF
usage: $0 options

Clean and remove host DNA sequences from FASTQ files.

OPTIONS:
   -h      Show help message
   -t      Number of threads (recommended: 8) [REQUIRED]
   -f      Forward or single-end fastq file (*.fastq.gz or *_1.fastq.gz) [REQUIRED]
   -c      BWA-indexed host genome (.fa) [REQUIRED]
   -r	   Reverse fastq file (*_2.fastq.gz) [OPTIONAL]
EOF
}

THREADS=
FASTQ_R1=
FASTQ_R2=
REF=

while getopts “ht:r:f:c:” OPTION
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
	 c)
	     REF=$OPTARG
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

echo "$(timestamp) [ metagen-fastqc ] Parsing command-line"
# check if all required arguments are supplied
if [[ -z ${THREADS} ]] || [[ -z ${FASTQ_R1} ]] || [[ -z ${REF} ]]
then
     echo "ERROR : Please supply required arguments"
     usage
     exit 1
fi

if [ ${THREADS} -eq 1 ]
then
    THREADS_SAM=1
else
    THREADS_SAM=$((${THREADS}-1))
fi

if [[ ! -z ${FASTQ_R2} ]]
then
        echo "$(timestamp) [ metagen-fastqc ] Cleaning FASTQ files"
        trim_galore --paired ${FASTQ_R1} ${FASTQ_R2} -o $(dirname ${FASTQ_R1})
        name=${FASTQ_R1%%_1.fastq*}
        echo "$(timestamp) [ metagen-fastqc ] Mapping files to host genome: ${REF}"
        bwa mem -M -t ${THREADS} ${REF} ${name}_1_val_1.fq.gz ${name}_2_val_2.fq.gz | samtools view -@ ${THREADS_SAM} -f 12 -F 256 -uS - -o ${name}_both_unmapped.bam 
	samtools sort -@ ${THREADS_SAM} -n ${name}_both_unmapped.bam -o ${name}_both_unmapped_sorted.bam
	bedtools bamtofastq -i ${name}_both_unmapped_sorted.bam -fq ${name}_clean_1.fastq -fq2 ${name}_clean_2.fastq
	echo "$(timestamp) [ metagen-fastqc ] Compressing output files"
	gzip ${name}_clean_1.fastq
	gzip ${name}_clean_2.fastq
        echo "$(timestamp) [ metagen-fastqc ] Cleaning tmp files"
        rm -rf ${name}_both_unmapped.bam ${name}_both_unmapped_sorted.bam ${name}_*_trimmed.fq.gz ${name}_*_val_*.fq.gz ${name}_*.fastq.gz_*txt
else
        echo "$(timestamp) [ metagen-fastqc ] Cleaning FASTQ files"
        trim_galore ${FASTQ_R1} -o $(dirname ${FASTQ_R1})
        name=${FASTQ_R1%%.fastq*}
        echo "$(timestamp) [ metagen-fastqc ] Mapping files to host genome: ${REF}"
        bwa mem -M -t ${THREADS} ${REF} ${name}_trimmed.fq.gz | samtools view -@ ${THREADS_SAM} -f 4 -F 256 -uS - -o ${name}_unmapped.bam
        samtools sort -@ ${THREADS_SAM} -n ${name}_unmapped.bam -o ${name}_unmapped_sorted.bam
        bedtools bamtofastq -i ${name}_unmapped_sorted.bam -fq ${name}_clean.fastq
	echo "$(timestamp) [ metagen-fastqc ] Compressing output file"
	gzip ${name}_clean.fastq
        echo "$(timestamp) [ metagen-fastqc ] Cleaning tmp files"
        rm -rf ${name}_unmapped.bam ${name}_unmapped_sorted.bam ${name}_trimmed.fq.gz ${FASTQ_R1}_*txt
fi
