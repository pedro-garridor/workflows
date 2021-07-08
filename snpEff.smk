# Workflow for VCF annotation with snpEff

# SAMPLES = glob_wildcards('filtered/{sample}.vcf').sample

rule snpsift:
    # NOTE: launch this using 2 cores
    input:
        dbsnp='/home/bioinformatica/Documentos/Referencia/hg19/dbSNP/00-All.vcf.gz',
        vcf='filtered/{sample}.vcf'
    output:
        temp('rs/{sample}.vcf')
    # threads: 2
    shell:
        "mkdir -p rs; "
        "java -jar ~/Software/snpEff/SnpSift.jar annotate "
        "{input} > {output}"

rule snpeff:
    # NOTE: launch this using 2 cores
    input:
        'rs/{sample}.vcf'
    output:
        protected('VCF/{sample}.vcf')
    params:
        hg='hg19'
    # threads: 2
    shell:
        # RFE: mark snpEff useless output as temp on snpeff output
        "mkdir -p VCF; "
        "java -Xmx8g -jar ~/Software/snpEff/snpEff.jar {params.hg} "
        "{input} > {output}; "
        "rm -f snpEff_genes.txt snpEff_summary.html; "
        "rm -f VCF/{wildcards.sample}-filtered.vcf*"

'''
rule vcf2tsv:
    input:
        'VCF/{sample}.vcf'
    output:
        protected('variants/{sample}.tsv')
    shell:
        # NOTE: VCF-kit (vk) is not compatible with GATK conda env
        # if using this with GATK, 'conda deactivate' it manually
        "mkdir -p variants; "
        "vk vcf2tsv wide --print-header --ANN "
        "{input} > {output}"
'''

# rm -rf rs