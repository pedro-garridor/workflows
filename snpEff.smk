# Workflow for VCF annotation with snpEff

# SAMPLES = glob_wildcards('filtered/{sample}.vcf').sample

rule snpsift:
    # NOTE: launch 2x
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
    # NOTE: launch this x2 - x3
    input:
        'rs/{sample}.vcf'
    output:
        protected('Calls/{sample}.vcf')
    params:
        hg='hg19'
    # threads: 2
    shell:
        # RFE: mark snpEff useless output as temp on snpeff output
        "mkdir -p Calls; "
        "java -Xmx8g -jar ~/Software/snpEff/snpEff.jar {params.hg} "
        "{input} > {output}; "
        "rm -f snpEff_*; "
        "rm -f VCF/{wildcards.sample}-filtered.vcf*"

'''
rule vcf2tsv:
    input:
        'Calls/{sample}.vcf'
    output:
        protected('variants/{sample}.tsv')
    shell:
        # NOTE: VCF-kit (vk) is not compatible with GATK conda env
        # if using this with GATK, add 'conda deactivate' to shell
        "mkdir -p variants; "
        "vk vcf2tsv wide --print-header --ANN "
        "{input} > {output}"
'''