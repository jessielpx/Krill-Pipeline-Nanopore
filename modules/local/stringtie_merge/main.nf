process STRINGTIE_MERGE {

    tag "all_samples"

    container '/project/def-lefranco/peixiliu/NXF_SINGULARITY_CACHEDIR/ontresearch-wf-transcriptomes-shaaaf20a5a0e76f9e18bad21af639a6b69e4a31a2f.sif'

    cpus 4
    memory '120 GB'
    time '12h'

    publishDir "${params.outdir}/transcripts/merged", mode: 'copy'

    input:
    path assembled_gtfs
    path ref_annotation
    path ref_genome

    output:
    path "merged_transcriptome.gtf", emit: merged_gtf
    path "final_non_redundant_transcriptome.fasta", emit: transcriptome_fasta

    script:
    """
    set -euo pipefail

    mkdir -p query_annotations

    cp ${assembled_gtfs} query_annotations/

    stringtie --merge \
        -G ${ref_annotation} \
        -p ${task.cpus} \
        -o merged_transcriptome.gtf \
        query_annotations/*

    gffread \
        -g ${ref_genome} \
        -w final_non_redundant_transcriptome.raw.fasta \
        merged_transcriptome.gtf

    awk '
    function write_record() {
        if (header != "" && sequence != "" && !seen[id]++) {
            print header
            print sequence
        }
    }

    /^>/ {
        write_record()

        header = \$0
        id = \$1
        sub(/^>/, "", id)
        sequence = ""
        next
    }

    {
        line = \$0
        gsub(/[[:space:]]/, "", line)
        sequence = sequence line
    }

    END {
        write_record()
    }
    ' final_non_redundant_transcriptome.raw.fasta \
      > final_non_redundant_transcriptome.fasta

    rm -f final_non_redundant_transcriptome.raw.fasta
    """
}
