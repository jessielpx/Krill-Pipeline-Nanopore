nextflow.enable.dsl = 2

include { PYCOQC }                    from './modules/local/pycoqc/main'
include { FASTCAT }                    from './modules/local/fastcat/main'
include { BUILD_MINIMAP2_INDEX }       from './modules/local/build_minimap2_index/main'
include { PREPROCESS_ANNOTATION }      from './modules/local/preprocess_annotation/main'
include { MINIMAP2_ALIGN }             from './modules/local/minimap2_align/main'
include { SAMTOOLS_SORT }              from './modules/local/samtools_sort/main'
include { SAMTOOLS_INDEX }             from './modules/local/samtools_index/main'
include { STRINGTIE_ASSEMBLY }         from './modules/local/stringtie_assembly/main'
include { STRINGTIE_MERGE }            from './modules/local/stringtie_merge/main'
include { BUILD_TRANSCRIPTOME_INDEX }  from './modules/local/build_transcriptome_index/main'
include { MAP_TRANSCRIPTOME }          from './modules/local/map_transcriptome/main'
include { SALMON_QUANT }               from './modules/local/salmon_quant/main'

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
        file(params.ref_genome, checkIfExists: true)
    )

    annotation_gtf_ch = Channel.value(
        file(params.ref_annotation, checkIfExists: true)
    )

    /*
     * Build genome Minimap2 index
     */
    BUILD_MINIMAP2_INDEX(reference_fasta_ch)

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
     * Reuse genome index for every sample
     */
    genome_index_ch = BUILD_MINIMAP2_INDEX.out.first()

    /*
     * Align reads to the genome
     */
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
     * Reuse cleaned annotation for every sample
     */
    cleaned_gtf_ch = PREPROCESS_ANNOTATION.out.cleaned_gtf.first()

    /*
     * Assemble transcripts separately for each sample
     */
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
     * Reuse merged transcriptome FASTA
     */
    transcriptome_fasta_ch =
        STRINGTIE_MERGE.out.transcriptome_fasta.first()

    /*
     * Build Minimap2 index for the merged transcriptome
     */
    BUILD_TRANSCRIPTOME_INDEX(
        transcriptome_fasta_ch
    )

    /*
     * Reuse transcriptome index for every sample
     */
    transcriptome_index_ch =
        BUILD_TRANSCRIPTOME_INDEX.out.index.first()

    /*
     * Align reads to the merged transcriptome
     */
    MAP_TRANSCRIPTOME(
        merged_fastq_ch,
        transcriptome_index_ch
    )

    /*
     * Quantify transcript abundance
     */
    SALMON_QUANT(
        MAP_TRANSCRIPTOME.out.bam,
        transcriptome_fasta_ch
    )
}
