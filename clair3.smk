rule clair3:
    input:
        bam='/home/pedro/Documentos/{sample}.bam',
        ref='/home/pedro/Documentos/Referencia/hg38/hg38.fa'
    output:
        '{sample}/merge_output.vcf.gz'
    threads: 8
    shell:
        '''
        mkdir -p {wildcards.sample}; 
        cp {input.bam} /home/pedro/Documentos/{wildcards.sample}/input.bam; 
        cp {input.bam}.bai /home/pedro/Documentos/{wildcards.sample}/input.bam.bai; 
        cp {input.ref} /home/pedro/Documentos/{wildcards.sample}/ref.fa; 
        cp {input.ref}.fai /home/pedro/Documentos/{wildcards.sample}/ref.fa.fai; 
        singularity exec /media/bioinformatica/Bioinf/Software/clair3_latest.sif \
        /opt/bin/run_clair3.sh \
        --bam_fn=/home/pedro/Documentos/{wildcards.sample}/input.bam \
        --ref_fn=/home/pedro/Documentos/{wildcards.sample}/ref.fa \
        --threads={threads} \
        --platform="ont" \
        --model_path="/opt/models/r941_prom_sup_g506" \
        --output=/home/pedro/Documentos/{wildcards.sample} \
        --snp_min_af=0.4 \
        --indel_min_af=0.4 \
        --call_snp_only \
        --fast_mode \
        --remove_intermediate_dir; 
        rm /home/pedro/Documentos/{wildcards.sample}/input.bam; 
        rm /home/pedro/Documentos/{wildcards.sample}/input.bam.bai; 
        rm /home/pedro/Documentos/{wildcards.sample}/ref.fa; 
        rm /home/pedro/Documentos/{wildcards.sample}/ref.fa.fai
        '''
