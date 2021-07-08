# Workflows

Collection of Snakefiles for different purposes

These Snakefiles are renamed with a descriptive filename. However, will be necessary to change it back to *Snakefile* again before running the pipeline.

- **bwa.smk**: pipeline for reads realignment (FASTQ), giving a BAM file compatible with GATK.
- **gatkCNV.smk**: pipeline for sensitive CNVs calling from exome sequencing with GATK.
- **GATKsnp.SMK**: pipeline for germline SNP & indel variant calling using GATK best practices. Note this workflow must be accompanied of snpEff.smk, as its last rule will create temp files.
- **haplotype.smk**: pipeline for variant phasing with WhatsHap and HapCUT2.
- **samtoolsVC.smk**: pipeline for exome alignment (BWA-MEM) and variant calling (samtools call).
- **snpEff.smk**: pipeline for VCF annotation.
