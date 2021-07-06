# Workflow for germline variant calling with GATK

# https://gatk.broadinstitute.org/hc/en-us/articles/360035535932

# NOTE: activate GATK conda env before running the pipeline

SAMPLES = glob_wildcards('BAM/{sample}.bam').sample 

rule haplotype_caller:
    # NOTE: launch this x2
    input:
        ref='/home/bioinformatica/Documentos/Referencia/hg19/GRCh37/Sequence/WholeGenomeFasta/genome.fa',
        bam='BAM/{sample}.bam'
    output:
        temp('VCF/{sample}-unfiltered.vcf.gz')
    threads: 4
    shell:
        "java -Xmx8g -jar ~/Software/gatk-4.2.0.0/gatk-package-4.2.0.0-local.jar "
        "HaplotypeCaller "
        "-R {input.ref} "
        "-I {input.bam} "
        "-O {output}"

rule cnn_score_variants:
    # NOTE: 2D mode. For 1D, remove BAM & --tensor-type param
    # NOTE: launch 1 by 1
    input:
        vcf='VCF/{sample}-unfiltered.vcf.gz',
        ref='/home/bioinformatica/Documentos/Referencia/hg19/GRCh37/Sequence/WholeGenomeFasta/genome.fa'
        bam='BAM/{sample}.bam'
    output:
        temp('VCF/{sample}-annotated.vcf')
    threads: 8
    shell:
        "~/Software/gatk-4.2.0.0/gatk CNNScoreVariants "
        "-I {input.bam} "
        "-V {input.vcf} "
        "-R {input.ref} "
        "-O {output} "
        "--tensor-type read-tensor"

rule filter_variant_tranches:
    input:
        vcf='VCF/{sample}-annotated.vcf',
        dbsnp='~/Documentos/Referencia/hg19/GRCh37/Annotation/Archives/archive-2015-07-17-14-31-42/Variation/Homo_sapiens.vcf'
    output:
        protected('VCF/{sample}.vcf')
    shell:
        "~/Software/gatk-4.2.0.0/gatk FilterVariantTranches "
        "-V {input.vcf} "
        "--resource {input.dbsnp} "
        "--info-key CNN_1D "
        "--snp-tranche 99.95 "
        "--indel-tranche 99.4 "
        "-O {output}"

# TODO: annotate variants (Funconator or snpEff)