# Workflow for germline variant calling with GATK

# https://gatk.broadinstitute.org/hc/en-us/articles/360035535932

# NOTE: activate GATK conda env before running the pipeline

# SAMPLES = glob_wildcards('BAM/{sample, [^-]}.bam').sample

# RFE: make GATK conda env available only within this SMK file

rule haplotype_caller:
    # NOTE: launch this using 1/2 the cores
    # NOTE: add -L for interval restriction
    input:
        ref='/home/bioinformatica/Documentos/Referencia/hg19/GRCh37/Sequence/WholeGenomeFasta/genome.fa',
        bam='BAM/{sample}.bam'
    output:
        vcf=temp('Calls/{sample}.vcf.gz'),
        tbi=temp('Calls/{sample}.vcf.gz.tbi')
    threads: 4
    shell:
        "mkdir -p Calls; "
        "java -Xmx8g -jar ~/Software/gatk-4.2.0.0/gatk-package-4.2.0.0-local.jar "
        "HaplotypeCaller "
        "-R {input.ref} "
        "-I {input.bam} "
        # "-L 1 "
        "-O {output.vcf}"

rule cnn_score_variants:
    # NOTE: 2D mode. For 1D, remove BAM & --tensor-type param
    # NOTE: add -L for interval restriction
    input:
        vcf='Calls/{sample}.vcf.gz',
        ref='/home/bioinformatica/Documentos/Referencia/hg19/GRCh37/Sequence/WholeGenomeFasta/genome.fa',
        bam='BAM/{sample}.bam',
        tbi=temp('Calls/{sample}.vcf.gz.tbi')
    output:
        vcf=temp('annotated/{sample}.vcf'),
        index=temp('annotated/{sample}.vcf.idx')
    threads: 8
    priority: 1
    shell:
        "mkdir -p annotated; "
        "~/Software/gatk-4.2.0.0/gatk CNNScoreVariants "
        "-I {input.bam} "
        "-V {input.vcf} "
        "-R {input.ref} "
        # "-L 1 "
        "-O {output.vcf} "
        "--tensor-type read_tensor"

rule filter_variant_tranches:
    input:
        vcf='annotated/{sample}.vcf',
        dbsnp='/home/bioinformatica/Documentos/Referencia/hg19/dbSNP/00-All.vcf.gz'
    output:
        temp('filtered/{sample}.vcf')
    shell:
        "mkdir -p filtered; "
        "~/Software/gatk-4.2.0.0/gatk FilterVariantTranches "
        "-V {input.vcf} "
        "--resource {input.dbsnp} "
        "--info-key CNN_2D "
        "--snp-tranche 99.95 "
        "--indel-tranche 99.4 "
        "-O {output}"

# rm -rf Calls annotated filtered