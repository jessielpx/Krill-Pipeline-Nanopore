/*
 * Krill - Joint bambu analysis
 *
 * Aggregates all cohort BAM/BAI files, constructs a deterministic sample manifest.
 * Executes joint transcript discovery and gene/isoform quantification via Bioconductor `bambu`.
 * 
 * Offline Bioconductor configuration.
 *
 *  Input:
 *    - path bams, stageAs: 'bams/*': All cohort sorted BAMs (${sample}.sorted.bam)
 *    - path bais, stageAs: 'bams/*': Matching BAM indices (${sample}.sorted.bam.bai)
 *    - path reference_gtf: Contig-matched reference GTF
 *    - path reference_fasta: Reference genome FASTA
 *    - path run_script: Driver R script invoking `bambu::bambu()`
 *
 *  Output (published to `02_transcriptome/bambu_results`):
 *    - bambu_transcripts.rds: Serialized SummarizedExperiment object
 *    - extended_annotations.gtf: Reference + newly discovered isoform GTF
 *    - counts_transcript.txt: Isoform-level read count matrix
 *    - counts_gene.txt: Gene-level read count matrix
 *    - bambu_sample_qc.tsv: Per-sample bambu QC metrics
 *    - sample_manifest.tsv: Verified cohort sample manifest
 */

process BAMBU_JOINT {

    tag "cohort"

    cpus 8
    memory '128 GB'
    time '48h'

    publishDir "${params.outdir}/02_transcriptome",
        mode: 'copy'

    input:
    path bams, stageAs: 'bams/*'
    path bais, stageAs: 'bams/*'
    path reference_gtf
    path reference_fasta
    path run_script

    output:
    path "bambu_results", emit: results

    shell:
    '''
    set -euo pipefail

    module load StdEnv/2023 r/4.6.1

    export R_LIBS_USER="!{params.r_libs_user ?: ''}"

    # Use a locally downloaded Bioconductor version map.
    BIOC_CONFIG="!{params.bioconductor_config_file ?: ''}"

    export BIOCONDUCTOR_CONFIG_FILE="$BIOC_CONFIG"

    if [[ -z "$R_LIBS_USER" ]]; then
        echo "ERROR: Set r_libs_user in params.yaml" >&2
        exit 1
    fi

    if [[ ! -d "$R_LIBS_USER" ]]; then
        echo "ERROR: R library directory not found: $R_LIBS_USER" >&2
        exit 1
    fi

    if [[ -z "${BIOCONDUCTOR_CONFIG_FILE:-}" ]]; then
        echo "ERROR: Set bioconductor_config_file in params.yaml" >&2
        exit 1
    fi

    BIOC_CONFIG_PATH="${BIOCONDUCTOR_CONFIG_FILE#file://}"

    if [[ ! -s "$BIOC_CONFIG_PATH" ]]; then
        echo "ERROR: Bioconductor config not found: $BIOC_CONFIG_PATH" >&2
        exit 1
    fi

    export BIOCONDUCTOR_ONLINE_VERSION_DIAGNOSIS=FALSE

    test -d "$R_LIBS_USER"
    test -s "!{reference_gtf}"
    test -s "!{reference_fasta}"
    test -s "!{run_script}"

    # Verify local Bioconductor configuration and packages.
    Rscript --vanilla -e '
        cat(
            "Bioconductor:",
            as.character(BiocManager::version()),
            "\\n"
        )

        for (pkg in c(
            "bambu",
            "Rsamtools",
            "SummarizedExperiment"
        )) {
            stopifnot(requireNamespace(pkg, quietly=TRUE))
            cat(
                pkg,
                as.character(packageVersion(pkg)),
                "\\n"
            )
        }
    '

    # Collect BAM files.
    shopt -s nullglob

    bam_files=(bams/*.sorted.bam)

    if (( ${#bam_files[@]} == 0 )); then
        echo "ERROR: No BAM files found" >&2
        exit 1
    fi

    printf 'sample\\tbam\\n' > bam_manifest.tsv

    for bam in "${bam_files[@]}"; do

        test -s "$bam"

        if [[ ! -s "${bam}.bai" ]]; then
            echo "ERROR: Missing BAI for $bam" >&2
            exit 1
        fi

        sample=$(basename "$bam" .sorted.bam)

        printf '%s\\t%s\\n' \
            "$sample" \
            "$(pwd)/$bam" \
            >> bam_manifest.tsv

    done

    # Deterministic sample ordering.
    {
        head -n 1 bam_manifest.tsv
        tail -n +2 bam_manifest.tsv | LC_ALL=C sort -k1,1
    } > bam_manifest.sorted.tsv

    echo "===== BAM MANIFEST ====="
    cat bam_manifest.sorted.tsv

    # Joint discovery and quantification.
    Rscript --vanilla "!{run_script}" \
        bam_manifest.sorted.tsv \
        "!{reference_gtf}" \
        "!{reference_fasta}" \
        bambu_results \
        !{task.cpus}

    # Validate outputs.
    test -s bambu_results/bambu_transcripts.rds
    test -s bambu_results/bambu_sample_qc.tsv
    test -s bambu_results/sample_manifest.tsv

    test -s bambu_results/extended_annotations.gtf
    test -s bambu_results/counts_transcript.txt
    test -s bambu_results/counts_gene.txt

    echo "Joint bambu analysis completed."
    '''
}
