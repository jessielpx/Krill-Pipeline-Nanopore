/*
 * Krill - BAM indexing and alignment QC
 *
 * Validates coordinate-sorted BAM files, generates `.bai` index files, and computes alignment summary metrics via samtools flagstat and idxstats.
 * 
 * Input:
 *    - tuple val(meta), path(bam): Sample metadata & coordinate-sorted BAM (expected filename: ${sample}.sorted.bam)
 *
 * Output:
 *    - tuple val(meta), path("${sample}.sorted.bam"), path("${sample}.sorted.bam.bai")
 *    - tuple val(meta), path("${sample}.flagstat.txt"), path("${sample}.idxstats.txt")
 */

process SAMTOOLS_INDEX {

    tag "${meta.sample}"

    cpus 2
    memory '8 GB'
    time '2h'

    publishDir "${params.outdir}/01_preprocessing/alignment/${meta.sample}",
        mode: 'copy'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta),
          path("${meta.sample}.sorted.bam"),
          path("${meta.sample}.sorted.bam.bai"),
          emit: indexed_bam

    tuple val(meta),
          path("${meta.sample}.flagstat.txt"),
          path("${meta.sample}.idxstats.txt"),
          emit: qc

    script:
    """
    set -euo pipefail

    module load StdEnv/2023 samtools

    # Validate BAM
    test -s "${bam}"
    samtools quickcheck -v "${bam}"

    # Generate BAM index
    samtools index \\
        -@ ${task.cpus} \\
        -o "${meta.sample}.sorted.bam.bai" \\
        "${bam}"

    # Alignment QC
    samtools flagstat \\
        -@ ${task.cpus} \\
        "${bam}" \\
        > "${meta.sample}.flagstat.txt"

    samtools idxstats \\
        "${bam}" \\
        > "${meta.sample}.idxstats.txt"

    # Validate outputs
    test -s "${meta.sample}.sorted.bam.bai"
    test -s "${meta.sample}.flagstat.txt"
    test -s "${meta.sample}.idxstats.txt"

    echo "BAM indexing and QC completed: ${meta.sample}"
    """
}
