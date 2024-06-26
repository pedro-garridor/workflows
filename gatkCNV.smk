# Snakefile for sentitive CNV calling with GATK
# Part I: https://gatk.broadinstitute.org/hc/en-us/articles/360035531092
# Part II: https://gatk.broadinstitute.org/hc/en-us/articles/360035890011

CONTROLS = glob_wildcards('PON/{control}.bam').control
SAMPLES = glob_wildcards("Samples/{sample, [^/^\.]+}.bam").sample

rule all:
    input:
        expand('GATK/Plots/denoised/{sample}.denoised.png', sample = SAMPLES),
        expand('GATK/Plots/segments/{sample}.modeled.png', sample = SAMPLES)

rule REF_collect_read_counts:
    input:
        bam='PON/{control}.bam',
        intervals='GATK/targets.preprocessed.interval_list'
    output:
        'GATK/REF/readCounts/{control}.hdf5'
    shell:
        """
        mkdir -p GATK/REF/readCounts/; 
        gatk-4.2.0.0/gatk CollectReadCounts \
        -I {input.bam} \
        -L {input.intervals} \
        -imr OVERLAPPING_ONLY \
        -O {output}
        """

rule REF_pon:
    input:
        controls=expand('GATK/REF/readCounts/{control}.tsv', control = CONTROLS)
    output:
        'GATK/REF/pon.hdf5'
    params:
        files = lambda wildcards, input: ' -I '.join(input.controls)
    shell:
        """
        jre1.8.0_291/bin/java -Xmx6500m \
        -jar gatk-4.2.0.0/gatk-package-4.2.0.0-local.jar CreateReadCountPanelOfNormals \
        -I {params.files} \
        -O {output}
        """

rule read_counts:
    input:
        bam='Samples/{sample}.bam',
        intervals='GATK/targets.preprocessed.interval_list'
    output:
        protected('GATK/readCounts/{sample}.hdf5')
    shell:
        """
        mkdir -p GATK/readCounts; 
        gatk-4.2.0.0/gatk CollectReadCounts \
        -I {input.bam} \
        -L {input.intervals} \
        -imr OVERLAPPING_ONLY \
        -O {output}
        """

rule denoise_read_counts:
    input:
        sample='GATK/readCounts/{sample}.hdf5',
        pon='GATK/REF/pon.hdf5'
    output:
        standardized=protected('GATK/copyRatios/{sample}-standardized.tsv'),
        denoised=protected('GATK/copyRatios/{sample}-denoised.tsv')
    shell:
        """
        mkdir -p GATK/copyRatios; 
        gatk-4.2.0.0/gatk --java-options "-Xmx12g" DenoiseReadCounts \
        -I {input.sample} \
        --count-panel-of-normals {input.pon} \
        --standardized-copy-ratios {output.standardized} \
        --denoised-copy-ratios {output.denoised}
        """

rule plot_denoised_copy_ratios:
    input:
        standarized='GATK/copyRatios/{sample}-standardized.tsv',
        denoised='GATK/copyRatios/{sample}-denoised.tsv',
        dictionary='Homo_sapiens/Ensembl/GRCh37/Sequence/WholeGenomeFasta/genome.dict'
    output:
        protected('GATK/Plots/denoised/{sample}.denoised.png')
    shell:
        """
        mkdir -p GATK/Plots/denoised; 
        gatk-4.2.0.0/gatk PlotDenoisedCopyRatios \
        --standardized-copy-ratios {input.standarized} \
        --denoised-copy-ratios {input.denoised} \
        --sequence-dictionary {input.dictionary} \
        --minimum-contig-length 46709983 \
        --output Plots/denoised/ \
        --output-prefix {wildcards.sample}
        """

rule collect_allelic_counts:
    input:
        sample='Samples/{sample}.bam',
        intervals='GATK/targets.preprocessed.interval_list',
        ref='Homo_sapiens/Ensembl/GRCh37/Sequence/WholeGenomeFasta/genome.fa'
    output:
        'GATK/allelicCounts/{sample}.tsv'
    threads: 4
    shell:
        """
        mkdir -p allelicCounts; 
        gatk-4.2.0.0/gatk --java-options "-Xmx16g" CollectAllelicCounts \
        -L {input.intervals} \
        -I {input.sample} \
        -R {input.ref} \
        -O {output}
        """

rule model_segments:
    input:
        denoised='GATK/copyRatios/{sample}-denoised.tsv',
        allelic='GATK/allelicCounts/{sample}.tsv'
    output:
        segments=protected('GATK/modelSegments/{sample}.cr.seg'),
        final=protected('GATK/modelSegments/{sample}.modelFinal.seg'),
        allelic=protected('GATK/modelSegments/{sample}.hets.tsv)'
    threads: 12
    shell:
        """
        mkdir -p modelSegments; 
        gatk-4.2.0.0/gatk --java-options "-Xmx46g" ModelSegments \
        --denoised-copy-ratios {input.denoised} \
        --allelic-counts {input.allelic} \
        --output modelSegments \
        --output-prefix {wildcards.sample}
        """

rule copy_ratio_segments:
    input:
        'GATK/modelSegments/{sample}.cr.seg'
    output:
        'GATK/crSegments/{sample}.called.seg'
    shell:
        """
        mkdir -p GATK/crSegments; 
        gatk-4.2.0.0/gatk CallCopyRatioSegments \
        -I {input} \
        -O {output}
        """

rule plot_modeled_segments:
    input:
        denoised='GATK/copyRatios/{sample}-denoised.tsv',
        allelic='GATK/modelSegments/{sample}.hets.tsv',
        segments='GATK/modelSegments/{sample}.modelFinal.seg',
        dictionary='Homo_sapiens/Ensembl/GRCh37/Sequence/WholeGenomeFasta/genome.dict'
    output:
        protected('GATK/Plots/segments/{sample}.modeled.png')
    shell:
        """
        mkdir -p GATK/Plots/segments/; 
        gatk-4.2.0.0/gatk PlotModeledSegments \
        --denoised-copy-ratios {input.denoised} \
        --allelic-counts {input.allelic} \
        --segments {input.segments} \
        --sequence-dictionary {input.dictionary} \
        --minimum-contig-length 46709983 \
        --output GATK/Plots/segments/ \
        --output-prefix {wildcards.sample}
        """
