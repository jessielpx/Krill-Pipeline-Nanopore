/*
 * Krill Renew - Unified CORE Workflow
 *
 * Nanopore cDNA long-read RNA-seq
 *
 * 1. FASTQ preparation
 * 2. minimap2 genome alignment
 * 3. samtools indexing and QC
 * 4. bambu joint transcript discovery and quantification
 * 5. SQANTI3 structural QC
 * 6. Unified QC
 *
 * Native execution on Rorqual.
 * No containers.
 */

include { PREPROCESSING } from './preprocessing'
include { BAMBU_JOINT } from '../modules/local/bambu_joint/main'
include { SQANTI3_QC } from '../modules/local/sqanti3_qc/main'
include { UNIFIED_QC } from '../modules/local/unified_qc/main'


workflow CORE {

    take:
    samples_ch
    reference_fasta_ch
    annotation_gtf_ch

    main:

    /*
     * STEP 1: FASTQ preparation and genome alignment
     */

    PREPROCESSING(
        samples_ch,
        reference_fasta_ch,
        annotation_gtf_ch
    )


    /*
     * STEP 2: Collect indexed BAMs for joint bambu analysis
     */

    sorted_bams_ch = PREPROCESSING.out.indexed_bam
        .toSortedList { a, b ->
            a[0].sample.toString() <=> b[0].sample.toString()
        }
        .map { rows ->

            if (rows.isEmpty()) {
                error "No genome-aligned BAMs available for bambu"
            }

            def sample_names = rows.collect {
                it[0].sample.toString()
            }

            if (sample_names.size() != sample_names.toSet().size()) {
                error "Duplicate sample IDs detected"
            }

            rows
        }

    cohort_bams_ch = sorted_bams_ch.map { rows ->
        rows.collect { it[1] }
    }

    cohort_bais_ch = sorted_bams_ch.map { rows ->
        rows.collect { it[2] }
    }


    /*
     * STEP 3: Joint bambu discovery and quantification
     */

    bambu_script_ch = Channel.value(
        file(
            "${projectDir}/bin/run_bambu.R",
            checkIfExists: true
        )
    )

    BAMBU_JOINT(
        cohort_bams_ch,
        cohort_bais_ch,
        PREPROCESSING.out.cleaned_gtf,
        reference_fasta_ch,
        bambu_script_ch
    )


    /*
     * STEP 4: SQANTI3 structural quality control
     */

    SQANTI3_QC(
        BAMBU_JOINT.out.results,
        PREPROCESSING.out.cleaned_gtf,
        reference_fasta_ch
    )


    /*
     * STEP 5: Collect FASTQ QC directories
     *
     * FASTCAT output:
     * tuple(meta, merged_fastq, fastcat_stats_dir)
     */

    fastcat_qc_ch = PREPROCESSING.out.fastcat_stats
        .toSortedList { a, b ->
            a[0].sample.toString() <=> b[0].sample.toString()
        }
        .map { rows ->
            rows.collect { it[2] }
        }


    /*
     * STEP 6: Collect alignment flagstat files
     *
     * Alignment QC output:
     * tuple(meta, flagstat, idxstats)
     */

    alignment_flagstats_ch = PREPROCESSING.out.alignment_qc
        .toSortedList { a, b ->
            a[0].sample.toString() <=> b[0].sample.toString()
        }
        .map { rows ->
            rows.collect { it[1] }
        }


    /*
     * STEP 7: Unified QC
     *
     * Combine:
     *   FASTQ statistics
     *   samtools flagstat
     *   bambu sample QC
     *   SQANTI3 structural classification
     */

    qc_script_ch = Channel.value(
        file(
            "${projectDir}/bin/build_unified_qc.py",
            checkIfExists: true
        )
    )

    UNIFIED_QC(
        fastcat_qc_ch,
        alignment_flagstats_ch,
        BAMBU_JOINT.out.results,
        SQANTI3_QC.out.results,
        qc_script_ch
    )


    /*
     * CORE OUTPUTS
     */

    emit:

    // FASTQ
    merged_fastq = PREPROCESSING.out.merged_fastq
    fastcat_stats = PREPROCESSING.out.fastcat_stats

    // Genome alignment
    indexed_bam = PREPROCESSING.out.indexed_bam
    genome_alignment_qc = PREPROCESSING.out.alignment_qc

    // Reference
    cleaned_gtf = PREPROCESSING.out.cleaned_gtf
    genome_index = PREPROCESSING.out.genome_index
    annotation_qc = PREPROCESSING.out.annotation_qc

    // Joint transcriptome and quantification
    bambu_results = BAMBU_JOINT.out.results

    // Structural transcript QC
    sqanti3_results = SQANTI3_QC.out.results

    // Unified QC
    unified_qc = UNIFIED_QC.out.results
}
