# Workflow for read realignment using BWA

# https://gatk.broadinstitute.org/hc/en-us/articles/360035535912

# NOTE: activate GATK conda env before running the pipeline

SAMPLES = glob_wildcards('FASTQ/{sample}.fq.gz').sample

rule bwa_mem:
    input:
        ref='~/Documentos/Referencia/hg19/GRCh37/Sequence/BWAIndex/genome.fa',
        L1='FASTQ/{sample}_L1.fq.gz',
        L2='FASTQ/{sample}_L2.fq.gz'
    output:
        temp('BAM/{sample}_ungrouped.bam')
    threads: 8
    shell:
        "mkdir -p BAM; "
        # Illumina/454/IonTorrent < 70bp:
            # bwa aln ref.fa reads.fq > reads.sai; bwa samse ref.fa reads.sai reads.fq > aln-se.sam
        # Illumina/454/IonTorrent > 70bp:
            # bwa mem ref.fa read1.fq read2.fq > aln-pe.sam
        # Illumina paired-end < 70bp:
            # bwa aln ref.fa read1.fq > read1.sai; bwa aln ref.fa read2.fq > read2.sai
            # bwa sampe ref.fa read1.sai read2.sai read1.fq read2.fq > aln-pe.sam
        # PacBio subreads or ONT:
            # bwa mem -x pacbio ref.fa reads.fq > aln.sam
            # bwa mem -x ont2d ref.fa reads.fq > aln.sam
        "bwa mem -t {threads} {input} > {output} | samtools sort -@ {threads} -o {output}"

rule add_read_groups:
    input:
        temp('BAM/{sample}-ungrouped.bam')
    output:
       temp('BAM/{sample}-grouped.bam')
    shell:
        "java -jar ~/Software/picard.jar AddOrReplaceReadGroups "
        "-I {input} "
        "-O {output} "
        "-RGLB lib1 "
        "-RGPL ILLUMINA "
        "-RGPU unit1 "
        "-RGSM: {wildcards.sample}"

'''
rule mark_duplicates:
    # NOTE: DO NOT RUN if sequencing is amplicon-based
    # NOTE: launch this x2 - x3
    # NOTE: if running this rule, modify base_recalibrator & apply_recalibration input
    input:
        'BAM/{sample}-grouped.bam'
    output:
        bam=temp('BAM/{sample}-dedup.bam'),
        metrics=temp('BAM/{sample}.txt')
    threads: 2
    shell:
        "java -jar ~/Software/picard.jar MarkDuplicates "
        "-I {input} "
        "-O {output.bam} "
        "-M {output.metrics}; "
        "samtools sort {output.bam}"
'''

rule base_recalibrator:
    input:
        ref='~/Documentos/Referencia/hg19/GRCh37/Sequence/WholeGenomeFasta/genome.fa',
        bam='BAM/{sample}-grouped.bam',
        dbsnp='~/Documentos/Referencia/hg19/GRCh37/Annotation/Archives/archive-2015-07-17-14-31-42/Variation/Homo_sapiens.vcf'
    output:
        temp('BAM/{sample}.table')
    threads: 2
    shell:
        "~/Software/gatk-4.2.0.0/gatk BaseRecalibrator "
        "-I {input.bam} "
        "-R {input.ref} "
        "--known-sites {input.dbsnp} "
        "-O {output}"

rule apply_recalibration:
    input:
        ref='~/Documentos/Referencia/hg19/GRCh37/Sequence/WholeGenomeFasta/genome.fa',
        bam='BAM/{sample}-grouped.bam',
        table='BAM/{sample}.table'
    output:
        protected('BAM/{sample}.bam')
    shell:
        "~/Software/gatk-4.2.0.0/gatk ApplyBQSR "
        "-R {input.ref} "
        "-I {input.bam} "
        "--bqsr-recal-file {input.table} "
        "-O {output}"

