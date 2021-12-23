rule longshot:
    input:
        bam='BAM/{sample}.bam',
        ref='/home/pedro/Documentos/Referencia/hg38/hg38.fa'
    output:
        'minAlleleQual{qual}_AF{af}/{sample}.vcf'
    params:
        quality='{qual}',
        af='{af}'
    conda: 'longshot.yaml'
    # threads: 2
    shell:
        '''
        mkdir -p minAlleleQual{wildcards.qual}_AF{wildcards.af}; 
        longshot --bam {input.bam} \
        --ref {input.ref} \
        --out {output} \
        -a {params.quality} \
        -E {params.af}
        '''
