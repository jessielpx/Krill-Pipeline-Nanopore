/*
 * Krill - Minimap2 alignment
 *
 * Aligns merged Nanopore cDNA reads against the pre-built minimap2 genome index (.mmi) in splice-aware mode and streams directly into coordinate-sorted BAM format.
 *
 * Input:
 *    - tuple val(meta), path(reads): Sample metadata & merged FASTQ (.fastq.gz)
 *    - path minimap_index: Pre-built minimap2 splice index (.mmi)
 *
 * Output:
 *    - tuple val(meta), path("${sample}.sorted.bam"): Coordinate-sorted BAM
 *    - tuple val(meta), path("${sample}.minimap2.log"): Alignment stderr log
 */

process MINIMAP2_ALIGN {

    tag "${meta.sample}"

    cpus 8
    memory '32 GB'
    time '12h'

    publishDir "${params.outdir}/01_preprocessing/alignment/${meta.sample}",
        pattern: "*.minimap2.log",
        mode: 'copy'

    input:
    tuple val(meta), path(reads)
    path minimap_index

    output:
    tuple val(meta),
          path("${meta.sample}.sorted.bam"),
          emit: bam

    tuple val(meta),
          path("${meta.sample}.minimap2.log"),
          emit: log

    script:

    def sample = meta.sample.toString()
    def readType = params.read_type?.toString()?.trim()?.toLowerCase()

    if (!(sample ==~ /[A-Za-z0-9._-]+/)) {
        error "Invalid sample ID: ${sample}"
    }

    if (readType != 'cdna') {
        error "Krill Renew currently supports --read_type cdna"
    }

    """
    set -euo pipefail

    module load StdEnv/2023 minimap2 samtools

    test -s "${reads}"
    test -s "${minimap_index}"

    echo "minimap2 version:"
    minimap2 --version

    echo "samtools version:"
    samtools --version | head -n 1

    # Splice-aware alignment for Nanopore cDNA.
    # -ub searches both orientations for canonical splice sites.

    minimap2 \\
        -ax splice \\
        -ub \\
        -t 6 \\
        "${minimap_index}" \\
        "${reads}" \\
        2> "${sample}.minimap2.log" \\
    | samtools sort \\
        -@ 2 \\
        -m 2G \\
        -O BAM \\
        -o "${sample}.sorted.bam" -

    test -s "${sample}.sorted.bam"

    samtools quickcheck -v "${sample}.sorted.bam"

    echo "Genome alignment completed: ${sample}"
    """
}
