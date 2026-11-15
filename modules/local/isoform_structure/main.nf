/*
 * Krill Renew - Isoform Structure Annotation
 *
 * Uses:
 *   - Bambu extended transcript annotation
 *   - Isoform usage results
 *   - Ensembl 104 reference GTF
 *
 * Outputs:
 *   - Exon chain comparison
 *   - CDS and UTR annotation
 *   - Structural annotation QC
 *
 * Native execution on Rorqual.
 */

process ISOFORM_STRUCTURE {

    tag "cohort"

    cpus 2
    memory '16 GB'
    time '8h'

    publishDir "${params.outdir}/08_isoform_structure",
        mode: 'copy'

    input:
    path bambu_results
    path usage_results
    path reference_gtf
    path structure_script

    output:
    path 'isoform_structure_results', emit: results

    shell:
    '''
    set -euo pipefail

    # Load consistent Python environment
    module --force purge
    module load StdEnv/2023
    module load python/3.11.5

    unset PYTHONHOME PYTHONPATH
    hash -r

    python3 -c 'import csv, re, argparse; print("Python OK")'

    # Validate inputs
    BAMBU_GTF="!{bambu_results}/extended_annotations.gtf"
    REF_GTF="!{reference_gtf}"

    CANDIDATES="!{usage_results}/isoform_usage_gene_summary.tsv"
    USAGE="!{usage_results}/isoform_usage_supported.tsv"

    SCRIPT="!{structure_script}"

    for file in \
        "$BAMBU_GTF" \
        "$REF_GTF" \
        "$CANDIDATES" \
        "$USAGE" \
        "$SCRIPT"
    do
        test -s "$file" || {
            echo "ERROR: Missing input: $file" >&2
            exit 1
        }
    done

    # Verify candidate summary schema
    python3 - "$CANDIDATES" <<'PY'
import csv
import sys

with open(sys.argv[1], newline="") as f:
    reader = csv.DictReader(f, delimiter=chr(9))
    fields = set(reader.fieldnames or [])

required = {"gene_id", "top_transcript"}

if not required.issubset(fields):
    raise SystemExit(
        "ERROR: Candidate summary must contain "
        "gene_id and top_transcript"
    )

print("Candidate summary schema OK")
PY

    echo "===== KRILL ISOFORM STRUCTURE ====="

    python3 "$SCRIPT" \
        --candidates "$CANDIDATES" \
        --usage "$USAGE" \
        --reference-gtf "$REF_GTF" \
        --bambu-gtf "$BAMBU_GTF" \
        --outdir isoform_structure_results

    # Validate outputs
    test -s isoform_structure_results/isoform_structure_annotations.tsv
    test -s isoform_structure_results/isoform_structure_qc.tsv

    echo "===== STRUCTURE QC ====="
    cat isoform_structure_results/isoform_structure_qc.tsv

    echo "Krill ISOFORM_STRUCTURE completed."
    '''
}
