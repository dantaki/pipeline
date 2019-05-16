# pipeline
Get your pipe on
--------------------

## Guidelines (just a good rule of thumb)
[Functional equivalence of genome sequencing analysis pipelines enables harmonized variant calling across human genetics projects](https://www.nature.com/articles/s41467-018-06159-4)

-----------

## Pipeline

1. Raw Reads (FASTQ)
  * QC with [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
  * Paired-end sequencing has two FASTQ files `R1.fq` and `R2.fq`
  * You can compress the files with `gzip` to save space
  
2. Align with [bwa mem](http://bio-bwa.sourceforge.net/)

```
bwa mem <REF.fasta> <R1.fq> <R2.fq> -K 100000000 -Y -R 'RG\tID:sequencing_run_id\tSM:sample_name\tLB:library_prep_id\tPL:ILLUMINA' 
```

  * outputs to STDOUT as SAM
  * thread with the `-t THREADS` option
  
3. Convert SAM -> BAM
  * can do with with a pipe with bwa
  * `bwa mem ..... | samtools view -bh >out.bam`
  * thread with the `-@ THREADS` option
  
4. Mark Duplicates
  * Some options here. [Picard tools](https://broadinstitute.github.io/picard/) is the go-to standard
  * [Sambamba](http://lomereiter.github.io/sambamba/) is faster but not as well tested. 
