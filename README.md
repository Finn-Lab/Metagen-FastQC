Metagen-FastQC
==============
Cleans metagenomic reads to remove adapters, low-quality bases and human contamination:

## Dependencies:
* trim_galore: https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/
* BWA: https://github.com/lh3/bwa/releases
* samtools: http://www.htslib.org/download/
* bedtools: https://github.com/arq5x/bedtools2/releases

For human decontamination download the hg38 indexed genome here.

## How to run:

Step 1 (If indexing your own host genome):
```
bwa index host_genome.fa
```

Step 2 (or Step 1 if using the above pre-indexed human genome):
```
$ metagen-fastqc.sh -t 8 -f input_1.fastq(gz) -r input_2.fastq(gz) -c host_genome.fa
```

<b>Notes:</b>
* <b>-t</b> controls the number of threads. Going above 8 does not significantly improve performance.
* Cleaned files will be generated in the same directory where the original FASTQ files are located and suffixed with "_clean.fastq.gz".
* <b>-f</b> argument can be either the forward read file or just a single-end FASTQ file (in the latter case -r would be omitted). When using paired-end files, make sure the forward and reverse files end in _1.fastq.gz and _2.fastq.gz, respectively.
