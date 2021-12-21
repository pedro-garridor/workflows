samples = glob_wildcards("/home/bioinf/Documentos/Proyectos/Abiertos/mutacionesRecurrentes/Originales/BAM/{sample}.vcf").sample

callers = ['hapcut2', 'whatshap']

# samples = ["09-4595.IonXpress_088", "10-382.IonXpress_087"]

# samples = "09-4595.IonXpress_088"

# wildcard_constraints:
#     sample = ".IonXpress_\d+$"

rule all:
    input:
        expand("{caller}/{sample}.tsv", caller = callers, sample = samples)

rule snpeff:
    input:
        "/home/bioinf/Documentos/Proyectos/Abiertos/mutacionesRecurrentes/Originales/BAM/{sample}.vcf"
    output:
        "vcf/{sample}.vcf"
    threads: 2
    shell:
        """
        java -Xmx8g -jar /home/bioinf/Software/snpEff/snpEff.jar hg19 {input} > {output}
        """

rule hairs:
    input:
        bam = "/home/bioinf/Documentos/Proyectos/Abiertos/mutacionesRecurrentes/Originales/BAM/{sample}.bam",
        vcf = "vcf/{sample}.vcf"
    output:
        temp('hapcut2/hairs/{sample}')
    shell:
        """
        /home/bioinf/Software/HapCUT2/build/extractHAIRS \
        --bam {input.bam} \
        --VCF {input.vcf} \
        --out {output}
        """

rule hapcut2:
    input:
        vcf = "vcf/{sample}.vcf",
        fragments = "hapcut2/hairs/{sample}"
    output:
        "hapcut2/{sample}.vcf"
    shell:
        """
        /home/bioinf/Software/HapCUT2/build/HAPCUT2 \
        --fragments {input.fragments} \
        --VCF {input.vcf} \
        --out {wildcards.sample}; 
        mv {wildcards.sample}.phased.VCF {output}; 
        rm {wildcards.sample}
        """

rule whatshap:
    input:
        ref = "/home/bioinf/Documentos/Referencia/hg19/hg19.fa",
        vcf = "vcf/{sample}.vcf",
        bam = "/home/bioinf/Documentos/Proyectos/Abiertos/mutacionesRecurrentes/Originales/BAM/{sample}.bam"
    output:
        'whatshap/{sample}.vcf'
    shell:
        "whatshap phase -o {output} --reference {input}"

rule tabulator:
    input:
        "{caller}/{sample}.vcf"
    output:
        '{caller}/{sample}.tsv'
    shell:
        """
        python3 tabulador.py \
        {input} > {output}; 
        sed -i 1d {output}
        """
