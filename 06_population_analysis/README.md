###1.Call variants of samples using GATK.

note:
1. All input should be full path.
2. The name of ref_genome should be *.fasta.
3. The name of reseq data should be *.R1.fastq.gz *.R2.fastq.gz

``````
usage: bash call_variations.sh ref_genome data_directory thread_number
``````

###2.Merge variants of different samples into a single vcf file.

``````
merge_samples.sh
``````
###3.Filter variants and divide result into SNPs and InDels.
```
filter_variations.sh
```
