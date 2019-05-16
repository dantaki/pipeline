# pipeline
Get your pipe on
--------------------

### Guidelines (just a good rule of thumb)

[*Functional equivalence of genome sequencing analysis pipelines enables harmonized variant calling across human genetics projects*](https://www.nature.com/articles/s41467-018-06159-4)

-----------

## Pipeline

1. **Raw Reads (FASTQ)**
  * QC with [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
  * Paired-end sequencing has two FASTQ files `R1.fq` and `R2.fq`
  * You can compress the files with `gzip` to save space
  
2. **Align with [bwa mem](http://bio-bwa.sourceforge.net/)**

```
bwa mem <REF.fasta> <R1.fq> <R2.fq> -K 100000000 -Y -R 'RG\tID:sequencing_run_id\tSM:sample_name\tLB:library_prep_id\tPL:ILLUMINA' 
```

  * outputs to STDOUT as SAM
  * thread with the `-t THREADS` option
  * the `-R` option is important! GATK and other downstream tools require the Read Group information.
    * ID: sequencing run ID
    * SM: sample name 
    * LB: library prep ID
    * PL: Platform (ILLUMINA)
  
3. **Convert SAM -> BAM**
  * can do with with a pipe with bwa
  * `bwa mem ..... | samtools view -bh >out.bam`
  * thread with the `-@ THREADS` option
  
4. **Mark Duplicates**
  * Some options here. [Picard tools](https://broadinstitute.github.io/picard/) is the go-to standard
  * [Sambamba](http://lomereiter.github.io/sambamba/) is faster but not as well tested. 
  * This step identifies and removes (unless needed) PCR duplicates
  
5. **Sort by Coordinate System** 
  * Lots of options here too. 
  * For larger alignment files, use a temporary directory on scratch (local or OASIS)
  * For super-large alignment files, use temporary directory on OASIS scratch and pray the sys-admin does not notice
  
  * `samtools sort` 
  * `sambamba sort`

6. **GATK Indel Realignment**

7. **[GATK](https://software.broadinstitute.org/gatk/) BQSR**
  * Base Quality Score Recalibration
  * This step takes a very long time
  * Two-steps
    * Build the BQSR model
    * Adjust the base quality score for EACH READ!!! (this is the time suck) 


8. **Variant Calling**
  * HaplotypeCaller
  * VQSR: Vartiant Quality Score Recalibration 

.... 

  
