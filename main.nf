#!/usr/bin/env nextflow

params.input_dir = "$baseDir/data"
params.genome_dir = "$baseDir/genome"
params.output_dir = "$baseDir/results"

process fastqc {
    label 'light'  

    input:
    file fastq_files from file("${params.input_dir}/*.fastq.gz")

    output:
    file "FastQC" into fastqc_results

    script:
    """
    mkdir -p FastQC
    fastqc ${fastq_files} --outdir FastQC
    """
}

process trimming {
    label 'moderate'  

    input:
    file fastq_files from file("${params.input_dir}/*_R1_001.fastq.gz")
    
    output:
    file "FastqTrimmed" into trimmed_fastq

    script:
    """
    mkdir -p FastqTrimmed
    for file in ${fastq_files}; do
        base=\$(basename "\$file" _R1_001.fastq.gz)
        java -jar Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 4 -phred33 \
        "\$base"_R1_001.fastq.gz "\$base"_R2_001.fastq.gz \
        FastqTrimmed/"\$base"_R1_trimmed.fastq.gz FastqTrimmed/"\$base"_R1_unpaired.fastq.gz \
        FastqTrimmed/"\$base"_R2_trimmed.fastq.gz FastqTrimmed/"\$base"_R2_unpaired.fastq.gz \
        ILLUMINACLIP:Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    done
    """
}

process genome_indexing {
    label 'memory_intensive'  

    input:
    file genome_fasta from file("${params.genome_dir}/genome.fa")
    file genome_gtf from file("${params.genome_dir}/annotation.gtf")

    output:
    file "STARIndex" into star_index

    script:
    """
    mkdir -p STARIndex
    STAR --runMode genomeGenerate --genomeDir STARIndex --genomeFastaFiles ${genome_fasta} \
    --sjdbGTFfile ${genome_gtf} --runThreadN 4
    """
}

process first_pass_alignment {
    label 'heavy'  

    input:
    file fastq_trimmed from trimmed_fastq
    file star_idx from star_index

    output:
    file "STARAlignmentsPass1" into first_pass_bam

    script:
    """
    mkdir -p STARAlignmentsPass1
    for file in ${fastq_trimmed}; do
        base=\$(basename "\$file" _R1_trimmed.fastq.gz)
        STAR --runThreadN 8 --genomeDir ${star_idx} \
        --readFilesIn \$file FastqTrimmed/"\$base"_R2_trimmed.fastq.gz \
        --outFileNamePrefix STARAlignmentsPass1/"\$base"_firstpass_ --outSAMtype BAM SortedByCoordinate
    done
    """
}

process second_pass_alignment {
    label 'intensive'  

    input:
    file bam_files from first_pass_bam
    file star_idx from star_index

    output:
    file "STARAlignmentsPass2" into second_pass_bam

    script:
    """
    mkdir -p STARAlignmentsPass2
    for file in FastqTrimmed/*.fastq.gz; do
        base=\$(basename "\$file" _R1_trimmed.fastq.gz)
        STAR --runThreadN 8 --genomeDir ${star_idx} \
        --readFilesIn \$file FastqTrimmed/"\$base"_R2_trimmed.fastq.gz \
        --sjdbFileChrStartEnd STARAlignmentsPass1/*.tab \
        --outFileNamePrefix STARAlignmentsPass2/"\$base"_ --outSAMtype BAM SortedByCoordinate
    done
    """
}

process mappability_stats {
    label 'light'  

    input:
    file bam_logs from second_pass_bam

    output:
    file "AlignmentStats" into alignment_stats

    script:
    """
    python3 alignment_analysis.py STARAlignmentsPass2 AlignmentStats
    """
}

process bam_processing {
    label 'moderate'  

    input:
    file bam_files from second_pass_bam

    output:
    file "BAMDuplicatesRemoved" into processed_bam
    file "AlignmentStats" into alignment_stats

    script:
    """
    mkdir -p BAMDuplicatesRemoved
    for file in ${bam_files}; do
        base=\$(basename "\$file" .bam)
        java -jar picard.jar MarkDuplicates INPUT=\$file OUTPUT=BAMDuplicatesRemoved/\$base"_marked.bam" \
        METRICS_FILE=BAMDuplicatesRemoved/\$base"_metrics.txt" REMOVE_DUPLICATES=true
        samtools index BAMDuplicatesRemoved/\$base"_marked.bam"
    done
    
    # Plot duplicate metrics
    python3 duplicate_metrics_plot.py BAMDuplicatesRemoved AlignmentStats
    """
}

process merge_bam {
    label 'heavy'  

    input:
    file bam_files from processed_bam

    output:
    file "MergedBam" into merged_bam

    script:
    """
    mkdir -p MergedBam
    samtools merge MergedBam/merged.bam BAMDuplicatesRemoved/*.bam
    samtools index MergedBam/*.bam
    """
}

process read_coverage {
    label 'moderate'  

    input:
    file merged_bam from merged_bam

    output:
    file "CoverageStats" into coverage_results

    script:
    """
    plotCoverage -b MergedBam/*.bam -o CoverageStats/read_coverage.png --outRawCounts read_coverage.tab

    # Generate read coverage plot
    python3 calculate_coverage.py MergedBam CoverageStats 4
    """
}

process feature_counts {
    label 'moderate'  

    input:
    file bam_files from processed_bam

    output:
    file "FeatureCounts" into count_matrix

    script:
    """
    mkdir -p FeatureCounts
    featureCounts -T 4 -a ${params.genome_dir}/annotation.gtf -o FeatureCounts/gene_counts.txt BAMDuplicatesRemoved/*.bam
    """
}

process differential_expression {
    label 'moderate'

    input:
    file counts from count_matrix

    output:
    file "DifferentialExpression" into diff_expr_results

    script:
    """
    mkdir -p DifferentialExpression
    Rscript scriptDESeq2.R
    """
}

process functional_enrichment {
    label 'light'

    input:
    file gene_list from diff_expr_results

    output:
    file "FunctionalEnrichment" into enrichment_results

    script:
    """
    mkdir -p FunctionalEnrichment
    Rscript scriptEnrichR.R
    """
}

workflow {
    fastqc()
    trimming()
    genome_indexing()
    first_pass_alignment()
    second_pass_alignment()
    mappability_stats()
    bam_processing()
    merge_bam()
    read_coverage()
    feature_counts()
    differential_expression()
    functional_enrichment()
}

