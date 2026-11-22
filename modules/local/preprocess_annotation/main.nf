/*
 * Krill - Reference annotation preparation
 *
 * Filters a gene annotation GTF file to retain only records on contigs present in the target reference genome FASTA
 * Preserve all original gene/transcript IDs and attributes.
 *
 * Input:
 *    - path annotation_gtf: Uncompressed reference annotation file (.gtf)
 *    - path reference_fasta: Uncompressed reference genome FASTA (.fa/.fasta)
 *
 * Output:
 *    - annotation.cleaned.gtf: Contig-synchronized GTF annotation
 *    - annotation_qc/reference_contigs.txt: Whitelist of contigs from FASTA
 *    - annotation_qc/summary.tsv: Counts of retained/removed features
 *    - annotation_qc/removed_contigs.tsv: Breakdown of dropped contigs
 */

process PREPROCESS_ANNOTATION {

    tag "reference_annotation"

    cpus 1
    memory '8 GB'
    time '2h'

    publishDir "${params.outdir}/01_preprocessing/reference/annotation",
        mode: 'copy'

    input:
    path annotation_gtf
    path reference_fasta

    output:
    path "annotation.cleaned.gtf", emit: cleaned_gtf
    path "annotation_qc", emit: qc

    script:
    """
    set -euo pipefail

    mkdir -p annotation_qc

    test -s "${annotation_gtf}"
    test -s "${reference_fasta}"

    # Extract contig names from the genome FASTA.
    awk '
        /^>/ {
            sub(/^>/, "", \$1)
            print \$1
        }
    ' "${reference_fasta}" \
        > annotation_qc/reference_contigs.txt

    test -s annotation_qc/reference_contigs.txt

    # Preserve GTF records only for contigs
    # present in the genome reference.
    awk -F '\\t' '
        NR == FNR {
            allowed[\$1] = 1
            next
        }

        /^#/ {
            print
            next
        }

        NF != 9 {
            print "ERROR: Invalid GTF at line " FNR > "/dev/stderr"
            bad = 1
            exit 2
        }

        {
            total++

            if (\$1 in allowed) {
                print
                retained++

                if (\$3 == "gene") genes++
                if (\$3 == "transcript") transcripts++
                if (\$3 == "exon") exons++
            } else {
                removed++
                removed_contigs[\$1]++
            }
        }

        END {
            if (bad) exit 2

            print "metric\\tvalue" > "annotation_qc/summary.tsv"
            print "total_records\\t" total+0 >> "annotation_qc/summary.tsv"
            print "retained_records\\t" retained+0 >> "annotation_qc/summary.tsv"
            print "removed_records\\t" removed+0 >> "annotation_qc/summary.tsv"
            print "retained_genes\\t" genes+0 >> "annotation_qc/summary.tsv"
            print "retained_transcripts\\t" transcripts+0 >> "annotation_qc/summary.tsv"
            print "retained_exons\\t" exons+0 >> "annotation_qc/summary.tsv"

            print "contig\\tremoved_records" > "annotation_qc/removed_contigs.tsv"

            for (contig in removed_contigs) {
                print contig "\\t" removed_contigs[contig] \
                    >> "annotation_qc/removed_contigs.tsv"
            }

            if (retained == 0 || exons == 0) {
                print "ERROR: No usable annotation retained" \
                    > "/dev/stderr"
                exit 3
            }
        }
    ' annotation_qc/reference_contigs.txt \
      "${annotation_gtf}" \
      > annotation.cleaned.gtf

    test -s annotation.cleaned.gtf
    test -s annotation_qc/summary.tsv

    echo "Reference annotation preparation completed."
    """
}
