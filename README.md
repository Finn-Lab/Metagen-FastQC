Metagen-FastQC
==============
Cleans metagenomic reads to remove adapters, low-quality bases and human contamination:

## Dependencies:
* trim_galore: https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/
* BWA: https://github.com/lh3/bwa/releases
* samtools: http://www.htslib.org/download/
* bedtools: https://github.com/arq5x/bedtools2/releases

## How to run (with 8 threads):

<b>Locally:</b>
```
$ metagen-fastqc.sh -t 8 -f input_1.fastq(gz) -r input_2.fastq(gz)
```

<b>In LSF:</b>
```
$ bsub -M 30000 -n 8 -q production-rh74 -o clean.log "metagen-fastqc.sh -t 8 -f input_1.fastq(gz) -r input_2.fastq(gz)"
```

<b>Notes:</b>
* Cleaned files will be generated in the same directory where the original FASTQ files are located and suffixed with "_clean.fastq.gz".
* <b>-f</b> argument can be either the forward read file or just a single-end FASTQ file (in the latter case -r would be omitted). When using paired-end files, make sure the forward and reverse files end in _1.fastq.gz and _2.fastq.gz, respectively.
