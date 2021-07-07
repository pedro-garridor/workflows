# Workflows

Collection of Snakefiles for different purposes

These Snakefiles are renamed with a descriptive filename. However, will be necessary to change it back to *Snakefile* again before running the pipeline.

- **bwa.smk**: *UNDER DEVELOPMENT, NEEDS VALIDATION*. Pipeline intended to realign reads (FASTQ), giving a BAM file compatible with GATK.
- **gatkCNV.smk**: Pipeline intended to sensitively call CNVs from exome sequencing with GATK.
- **GATKsnp.SMK**: *UNDER DEVELOPMENT, NEEDS VALIDATION*. Pipeline intended to call germline SNP & indel variants using GATK best practices.
- **haplotype.smk**: pipeline for variant phasing with WhatsHap and HapCUT2.
- **samtoolsVC.smk**: pipeline for exome alignment (BWA-MEM) and variant calling (samtools call).
- **snpEff.smk**: *UNDER DEVELOPMENT, NEEDS VALIDATION*. Pipeline intended for VCF annotation.
