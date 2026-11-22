/*
 * Krill - Build minimap2 genome index
 *
 * Reference: GRCh38
 * Sequencing: Nanopore cDNA RNA-seq
 * Builds a splice-aware minimap2 genome index (.mmi) from the reference FASTA for downstream long-read spliced alignment.
 *
 * Input:
 *    - path reference_fasta: Reference genome FASTA (e.g., GRCh38)
 *
 * Output:
 *    - path "reference.mmi": Pre-built minimap2 splice index file
 */

process BUILD_MINIMAP2_INDEX {

    tag "GRCh38"

    cpus 4
    memory '32 GB'
    time '4h'

    publishDir "${params.outdir}/01_preprocessing/reference/cdna",
        mode: 'copy'

    input:
    path reference_fasta

    output:
    path "reference.mmi"

    script:
    """
    set -euo pipefail

    module load StdEnv/2023 minimap2

    test -s "${reference_fasta}"

    echo "Building minimap2 genome index"
    minimap2 --version

    minimap2 \\
        -x splice \\
        -t ${task.cpus} \\
        -d reference.mmi \\
        "${reference_fasta}"

    test -s reference.mmi

    echo "Genome index completed"
    """
}
