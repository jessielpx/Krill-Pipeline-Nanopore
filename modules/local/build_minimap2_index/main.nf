/*
 * Krill Renew - Build minimap2 genome index
 *
 * Reference: GRCh38
 * Sequencing: Nanopore cDNA RNA-seq
 *
 * Native execution on Rorqual
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
