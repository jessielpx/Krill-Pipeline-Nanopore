process FASTCAT {

    tag "${meta.sample}"

    container '/project/def-lefranco/peixiliu/NXF_SINGULARITY_CACHEDIR/ontresearch-wf-common-sha72f3517dd994984e0e2da0b97cb3f23f8540be4b.sif'

    cpus 2
    memory '32 GB'
    time '6h'

    publishDir "${params.outdir}/fastcat/${meta.sample}", mode: 'copy'

    input:
    tuple val(meta), path(read_dir)

    output:
    tuple val(meta),
          path("merged_fastq/*.fastq.gz"),
          path("fastcat_stats"),
          emit: results

    script:
    """
    set -euo pipefail

    mkdir -p merged_fastq fastcat_stats

    fastcat \
        -s "${meta.sample}" \
        -f fastcat_stats/per-file-stats.tsv \
        -i fastcat_stats/per-file-runids.tsv \
        -l fastcat_stats/per-file-basecallers.tsv \
        --histograms histograms \
        -x \
        ${read_dir} \
    | bgzip -@ ${task.cpus} \
      > merged_fastq/${meta.sample}.fastq.gz

    mv histograms/* fastcat_stats/ 2>/dev/null || true

    awk '
        NR == 1 {
            for (i = 1; i <= NF; i++) {
                ix[\$i] = i
            }
        }
        NR > 1 {
            count += \$ix["n_seqs"]
        }
        END {
            print count + 0
        }
    ' fastcat_stats/per-file-stats.tsv \
    > fastcat_stats/n_seqs

    awk -F "\\t" '
        NR == 1 {
            for (i = 1; i <= NF; i++) {
                ix[\$i] = i
            }
        }
        NR > 1 && \$ix["run_id"] != "" {
            print \$ix["run_id"]
        }
    ' fastcat_stats/per-file-runids.tsv \
    | sort -u \
    > fastcat_stats/run_ids

    awk -F "\\t" '
        NR == 1 {
            for (i = 1; i <= NF; i++) {
                ix[\$i] = i
            }
        }
        NR > 1 && \$ix["basecaller"] != "" {
            print \$ix["basecaller"]
        }
    ' fastcat_stats/per-file-basecallers.tsv \
    | sort -u \
    > fastcat_stats/basecallers
    """
}
