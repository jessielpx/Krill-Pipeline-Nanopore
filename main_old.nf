nextflow.enable.dsl = 2

include { PYCOQC }                    from './modules/local/pycoqc/main'
include { FASTCAT }                   from './modules/local/fastcat/main'
include { BUILD_MINIMAP2_INDEX }      from './modules/local/build_minimap2_index/main'
include { PREPROCESS_ANNOTATION }     from './modules/local/preprocess_annotation/main'
include { MINIMAP2_ALIGN }            from './modules/local/minimap2_align/main'
include { SAMTOOLS_SORT }             from './modules/local/samtools_sort/main'
include { SAMTOOLS_INDEX }            from './modules/local/samtools_index/main'
include { STRINGTIE_ASSEMBLY }        from './modules/local/stringtie_assembly/main'
include { STRINGTIE_MERGE }           from './modules/local/stringtie_merge/main'
include { BUILD_TRANSCRIPTOME_INDEX } from './modules/local/build_transcriptome_index/main'
include { MAP_TRANSCRIPTOME }         from './modules/local/map_transcriptome/main'
include { SALMON_QUANT }              from './modules/local/salmon_quant/main'
include { MERGE_TRANSCRIPT_COUNTS }   from './modules/local/merge_transcript_counts/main'
include { ANNOTATE_AND_SUM_COUNTS }   from './modules/local/annotate_and_sum_counts/main'
include { DIFFERENTIAL_EXPRESSION }   from './modules/local/differential_expression/main'
include { FILTER_UNSTRANDED_ANNOTATION } from './modules/local/filter_unstranded_annotation/main'
include { GFFCOMPARE } from './modules/local/gffcompare/main'
include { PARSE_GFFCOMPARE } from './modules/local/parse_gffcompare/main'

workflow {

    /*
     * Check required inputs
     */
    if (!params.summary) {
        error "Please provide --summary /path/to/sequencing_summary.txt"
    }

    if (!params.fastq_dir) {
        error "Please provide --fastq_dir /path/to/barcode_directories"
    }

    if (!params.ref_genome) {
        error "Please provide --ref_genome /path/to/reference.fa"
    }

    if (!params.ref_annotation) {
        error "Please provide --ref_annotation /path/to/annotation.gtf"
    }


    /*
     * Run-level quality control
     */
    summary_ch = Channel.value(
        tuple(
            [id: params.run_id],
            file(params.summary, checkIfExists: true)
        )
    )

    PYCOQC(summary_ch)


    /*
     * Read sample metadata and locate barcode directories
     */
    samples_ch = Channel
        .fromPath(
            "${projectDir}/assets/samples.csv",
            checkIfExists: true
        )
        .splitCsv(header: true)
        .map { row ->

            def meta = [
                barcode  : row.barcode,
                sample   : row.sample,
                alias    : row.alias,
                condition: row.condition
            ]

            def barcode_dir = file(
                "${params.fastq_dir}/${row.barcode}",
                checkIfExists: true
            )

            tuple(meta, barcode_dir)
        }


    /*
     * Merge FASTQ files and generate Fastcat statistics
     */
    FASTCAT(samples_ch)


    /*
     * Reusable reference inputs
     */
    reference_fasta_ch = Channel.value(
        file(
            params.ref_genome,
            checkIfExists: true
        )
    )

    annotation_gtf_ch = Channel.value(
        file(
            params.ref_annotation,
            checkIfExists: true
        )
    )


    /*
     * Build genome Minimap2 index
     */
    BUILD_MINIMAP2_INDEX(
        reference_fasta_ch
    )


    /*
     * Clean reference annotation
     */
    PREPROCESS_ANNOTATION(
        annotation_gtf_ch,
        reference_fasta_ch
    )


    /*
     * Extract merged FASTQ files from FASTCAT output
     */
    merged_fastq_ch = FASTCAT.out.results.map {
        meta,
        fastq_files,
        stats_dir ->

        tuple(meta, fastq_files)
    }


    /*
     * Align reads to the genome
     */
    genome_index_ch = BUILD_MINIMAP2_INDEX.out

    MINIMAP2_ALIGN(
        merged_fastq_ch,
        genome_index_ch
    )


    /*
     * Sort genome-aligned SAM files
     */
    SAMTOOLS_SORT(
        MINIMAP2_ALIGN.out.sam
    )


    /*
     * Index genome-aligned BAM files
     */
    SAMTOOLS_INDEX(
        SAMTOOLS_SORT.out.bam
    )


    /*
     * Assemble transcripts separately for each sample
     */
    cleaned_gtf_ch = PREPROCESS_ANNOTATION.out.cleaned_gtf

    STRINGTIE_ASSEMBLY(
        SAMTOOLS_INDEX.out.indexed_bam,
        cleaned_gtf_ch
    )


    /*
     * Collect all sample-level assembled GTF files
     */
    assembled_gtfs_ch = STRINGTIE_ASSEMBLY.out.assembled_gtf
        .map { meta, assembled_gtf ->
            assembled_gtf
        }
        .collect()


    /*
     * Merge all sample transcriptomes
     */
    STRINGTIE_MERGE(
        assembled_gtfs_ch,
        cleaned_gtf_ch,
        reference_fasta_ch
    )

    /*
     * Compare merged transcriptome against the reference annotation
     */
    GFFCOMPARE(
        STRINGTIE_MERGE.out.merged_gtf,
        cleaned_gtf_ch
    )

    parse_gffcompare_script_ch = Channel.value(
        file(
            "${projectDir}/bin/parse_gffcompare_annotation.py",
            checkIfExists: true
        )
    )

    PARSE_GFFCOMPARE(
        GFFCOMPARE.out.annotated_gtf,
        parse_gffcompare_script_ch
    )


    /*
     * Build Minimap2 index for the merged transcriptome
     */
    transcriptome_fasta_ch =
        STRINGTIE_MERGE.out.transcriptome_fasta

    BUILD_TRANSCRIPTOME_INDEX(
        transcriptome_fasta_ch
    )


    /*
     * Align reads to the merged transcriptome
     */
    transcriptome_index_ch =
        BUILD_TRANSCRIPTOME_INDEX.out.index

    MAP_TRANSCRIPTOME(
        merged_fastq_ch,
        transcriptome_index_ch
    )


    /*
     * Quantify transcript abundance with Salmon
     */
    SALMON_QUANT(
        MAP_TRANSCRIPTOME.out.bam,
        transcriptome_fasta_ch
    )


    /*
     * Collect all sample-level Salmon outputs
     */
    transcript_count_files_ch = SALMON_QUANT.out.counts
        .map { meta, count_file ->
            count_file
        }
        .collect()


    /*
     * Load the local Salmon count-merging script
     */
    merge_script_ch = Channel.value(
        file(
            "${projectDir}/bin/merge_salmon_counts.py",
            checkIfExists: true
        )
    )


    /*
     * Merge transcript count and TPM matrices
     */
    MERGE_TRANSCRIPT_COUNTS(
        transcript_count_files_ch,
        merge_script_ch
    )


    /*
     * Load the transcript-to-gene annotation script
     */
    annotation_script_ch = Channel.value(
        file(
            "${projectDir}/bin/add_gene_ids_and_sum_counts.py",
            checkIfExists: true
        )
    )


    /*
     * Add gene IDs and generate a gene-level count matrix
     */
    ANNOTATE_AND_SUM_COUNTS(
        MERGE_TRANSCRIPT_COUNTS.out.counts,
        STRINGTIE_MERGE.out.merged_gtf,
        annotation_script_ch
    )


    /*
     * Load inputs required for differential expression analysis
     */
    de_sample_sheet_ch = Channel.value(
        file(
            "${projectDir}/assets/samples.csv",
            checkIfExists: true
        )
    )

    de_script_ch = Channel.value(
        file(
            "${projectDir}/bin/de_analysis.R",
            checkIfExists: true
        )
    )

    /*
     * Remove unstranded transcript records before DGE/DTU
     */
    FILTER_UNSTRANDED_ANNOTATION(
        STRINGTIE_MERGE.out.merged_gtf
    )

    /*
     * Run transcript filtering, DGE and DTU analysis
     */
    DIFFERENTIAL_EXPRESSION(
        de_sample_sheet_ch,
        MERGE_TRANSCRIPT_COUNTS.out.counts,
        FILTER_UNSTRANDED_ANNOTATION.out.stranded_gtf,
        de_script_ch
    )
}
