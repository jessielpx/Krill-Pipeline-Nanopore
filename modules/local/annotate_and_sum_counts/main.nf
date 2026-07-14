process ANNOTATE_AND_SUM_COUNTS {

    tag "all_samples"

    container '/project/def-lefranco/peixiliu/NXF_SINGULARITY_CACHEDIR/ontresearch-wf-common-sha72f3517dd994984e0e2da0b97cb3f23f8540be4b.sif'

    cpus 1
    memory '8 GB'
    time '2h'

    publishDir "${params.outdir}/quantification/merged", mode: 'copy'

    input:
    path transcript_counts
    path merged_gtf
    path annotation_script

    output:
    path "unfiltered_transcript_counts_with_genes.tsv", emit: transcript_counts_with_genes
    path "all_gene_counts.tsv", emit: gene_counts

    script:
    """
    python ${annotation_script} \
        --gtf ${merged_gtf} \
        --transcript-counts ${transcript_counts} \
        --transcript-output unfiltered_transcript_counts_with_genes.tsv \
        --gene-output all_gene_counts.tsv
    """
}
