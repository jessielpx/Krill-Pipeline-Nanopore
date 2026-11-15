/*
 * Krill Renew - 08 ISOFORM_STRUCTURE_ALL
 *
 * Annotation of all Bambu transcripts.
 *
 * No dependency on ISOFORM_USAGE or DTU.
 */

process ISOFORM_STRUCTURE_ALL {

    tag "all_transcripts"

    cpus 2
    memory '16 GB'
    time '8h'

    publishDir "${params.outdir}/08_isoform_structure_all",
        mode: 'copy'

    input:
    path bambu_results
    path reference_gtf
    path structure_script

    output:
    path 'isoform_structure_results', emit: results

    shell:
    '''
    set -euo pipefail

    module --force purge
    module load StdEnv/2023
    module load python/3.11.5

    unset PYTHONHOME PYTHONPATH
    hash -r

    BAMBU_GTF="!{bambu_results}/extended_annotations.gtf"
    REF_GTF="!{reference_gtf}"
    SCRIPT="!{structure_script}"

    for required in \
        "$BAMBU_GTF" \
        "$REF_GTF" \
        "$SCRIPT"
    do
        test -s "$required" || {
            echo "Missing input: $required" >&2
            exit 1
        }
    done

    python3 "$SCRIPT" \
        --bambu-gtf "$BAMBU_GTF" \
        --reference-gtf "$REF_GTF" \
        --outdir isoform_structure_results

    test -s \
        isoform_structure_results/isoform_structure_annotations.tsv

    test -s \
        isoform_structure_results/isoform_structure_qc.tsv

    echo "ISOFORM_STRUCTURE_ALL completed"
    '''
}
