# Workflow for read realignment using BWA

# https://gatk.broadinstitute.org/hc/en-us/articles/360035535912

# NOTE: activate GATK conda env before running the pipeline

# SAMPLES = glob_wildcards('FASTQ/{sample}_L1.fq.gz').sample

rule bwa_mem:
    input:
        ref='/home/bioinformatica/Documentos/Referencia/hg19/GRCh37/Sequence/BWAIndex/version0.6.0/genome.fa',
        L1='FASTQ/{sample}_1.fastq.gz',
        L2='FASTQ/{sample}_2.fastq.gz'
    output:
        temp('ungrouped/{sample}.bam')
    threads: 8
    priority: 2
    shell:
        "mkdir -p ungrouped; "
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
        "bwa mem -t {threads} {input} | "
        "samtools sort -@ {threads} -o {output}"

rule add_read_groups:
    input:
        'ungrouped/{sample}.bam'
    output:
       temp('grouped/{sample}.bam')
    shell:
        "mkdir -p grouped; "
        "java -jar /home/bioinformatica/Software/picard.jar AddOrReplaceReadGroups "
        "-I {input} "
        "-O {output} "
        "-RGLB lib1 "
        "-RGPL ILLUMINA "
        "-RGPU unit1 "
        "-RGSM {wildcards.sample}"


rule mark_duplicates:
    # NOTE: DO NOT RUN if sequencing is amplicon-based
    # NOTE: launch this x2 - x3
    # NOTE: if not running this rule, modify base_recalibrator & apply_recalibration input
    input:
        'grouped/{sample}.bam'
    output:
        bam=temp('dedup/{sample}.bam'),
        metrics=temp('dedup/{sample}.txt')
    threads: 2
    shell:
        "mkdir -p dedup; "
        "java -jar /home/bioinformatica/Software/picard.jar MarkDuplicates "
        "-I {input} "
        "-O {output.bam} "
        "-M {output.metrics}"


rule base_recalibrator:
    input:
        ref='/home/bioinformatica/Documentos/Referencia/hg19/GRCh37/Sequence/WholeGenomeFasta/genome.fa',
        bam='dedup/{sample}.bam',
        dbsnp='/home/bioinformatica/Documentos/Referencia/hg19/dbSNP/00-All.vcf.gz'
    output:
        temp('grouped/{sample}.table')
    params:
        max_cycle=600
    shell:
        "/home/bioinformatica/Software/gatk-4.2.0.0/gatk BaseRecalibrator "
        "-I {input.bam} "
        "-R {input.ref} "
        "--known-sites {input.dbsnp} "
        "-max-cycle {params} "
        "-O {output}"

rule apply_recalibration:
    input:
        ref='/home/bioinformatica/Documentos/Referencia/hg19/GRCh37/Sequence/WholeGenomeFasta/genome.fa',
        bam='dedup/{sample}.bam',
        table='grouped/{sample}.table'
    output:
        bam=protected('BAM/{sample}.bam'),
        bai='BAM/{sample}.bai'
    shell:
        "mkdir -p BAM; "
        "/home/bioinformatica/Software/gatk-4.2.0.0/gatk ApplyBQSR "
        "-R {input.ref} "
        "-I {input.bam} "
        "--bqsr-recal-file {input.table} "
        "-O {output.bam}"

# rm -rf ungrouped dedup grouped