# Snakefile para llamado de variantes de exoma

samples = glob_wildcards("FASTQ/{sample}_L1_1.fq.gz").sample

rule all:
    input:
        expand("variants/{sample}.tsv", sample = samples)
    shell:
        "if [ -d unsorted ]; then rm -Rf unsorted; fi;"
        "if [ -d mpileup ]; then rm -Rf mpileup; fi; "
        "if [ -d vcf_filter ]; then rm -Rf vcf_filter; fi; "
        "if [ -d annotated ]; then rm -Rf annotated; fi"

rule bwa_mem:
    input:
        # ref="/media/genomica/Datos/nanopore/referencia/hg38/hg38.fa",
        # FIXME: CAMBIAR RUTA A hg19
        ref="/media/genomica/Datos/nanopore/referencia/hg38/hg38.fa",
        L1="FASTQ/{sample}_L1_1.fq.gz",
        L2="FASTQ/{sample}_L1_2.fq.gz"
    output:
        temp("unsorted/{sample}.bam")
    threads: 24
    shell:
        # TODO: añadir BED con las posiciones del exoma
        "bwa mem -t {threads} {input} > {output}"

rule sort_bam:
    input:
        "unsorted/{sample}.bam"
    output:
        protected("BAM/{sample}.bam")
    shell:
        "samtools sort {input} > {output}; "
        "samtools index {output}"

rule mpileup:
    input:
        # ref="/media/genomica/Datos/nanopore/referencia/hg38/hg38.fa",
        # FIXME: CAMBIAR RUTA A hg19
        ref="/media/genomica/Datos/nanopore/referencia/hg38/hg38.fa",
        bam="BAM/{sample}.bam"
    output:
        temp("mpileup/{sample}.mpileup")
    params:
        d="8000"
    # NOTE: se lanzan de 3 en 3 para evitar ocupar todo el disco
    threads: 8
    shell:
        # TODO: añadir BED con las posiciones del exoma
        "bcftools mpileup -f {input} --threads {threads} -d {params} -Ou -o {output}"

rule call:
    input:
        "mpileup/{sample}.mpileup"
    output:
        protected("BCF/{sample}.bcf.gz")
    params:
        threads="1"
    threads: 1
    shell:
        "bcftools call -mv --threads {threads} {input} -Ob -o {output}"

rule filter_vcf:
    input:
        "BCF/{sample}.bcf.gz"
    output:
        temp("vcf_filter/{sample}.vcf.gz")
    threads: 1
    shell:
        # TODO: añadir BED con las posiciones del exoma
        "bcftools filter -i 'QUAL>100 && DP>20' {input} --threads {threads} -Oz -o {output}; "
        "tabix -p vcf {output}"

rule snpsift:
    input:
        dbsnp="00-All.vcf.gz",
        vcf="vcf_filter/{sample}.vcf.gz"
    output:
        temp("annotated/{sample}.vcf")
    # NOTE: se lanzan de 8 en 8 para no colapsar la RAM
    threads: 3
    shell:
        "java -jar snpEff/SnpSift.jar annotate {input} > {output}"

rule snpeff:
    input:
        "annotated/{sample}.vcf"
    output:
        protected("VCF/{sample}.vcf")
    params:
        hg="hg19"
    # NOTE: se lanzan de 8 en 8 para no colapsar la RAM
    threads: 3
    shell:
        "java -Xmx8g -jar snpEff/snpEff.jar {params.hg} {input} > {output}"

rule vcf2tsv:
    input:
        "VCF/{sample}.vcf"
    output:
        protected("variants/{sample}.tsv")
    shell:
        "vk vcf2tsv wide --print-header --ANN {input} > {output}"
