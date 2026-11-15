/*
 * Krill Renew - SQANTI3 Quality Control
 *
 * Structural classification of bambu transcripts:
 * FSM, ISM, NIC, NNC, etc.
 *
 * ORF prediction enabled using TD2.
 *
 * SQANTI3 v6.0.2
 * Native execution on Rorqual.
 * No containers.
 */

process SQANTI3_QC {

    tag "cohort"

    cpus 8
    memory '64 GB'
    time '48h'

    publishDir "${params.outdir}/03_sqanti3",
        mode: 'copy'

    input:
    path bambu_results
    path reference_gtf
    path reference_fasta

    output:
    path "sqanti3_results", emit: results

    shell:
    '''
    set -euo pipefail

    # ==========================================
    # 1. SQANTI3 installation
    # ==========================================

    SQANTI_DIR="!{params.sqanti_dir ?: ''}"

    SQANTI_ENV="!{params.sqanti_env ?: ''}"

    if [[ -z "$SQANTI_DIR" || -z "$SQANTI_ENV" ]]; then
        echo "ERROR: Configure sqanti_dir and sqanti_env" >&2
        exit 1
    fi

    # ==========================================
    # 2. Clean Python and R environment
    # ==========================================

    unset PYTHONHOME PYTHONPATH PYTHONUSERBASE
    unset R_LIBS_USER R_HOME

    export PATH="$SQANTI_ENV/bin:$PATH"
    export PYTHONNOUSERSITE=1

    # Prefer micromamba C++ runtime over system
    # /lib64/libstdc++.so.6.
    #
    # Required for scipy.optimize._highs_wrapper
    # and TD2 ORF prediction.

    export LD_LIBRARY_PATH="$SQANTI_ENV/lib:${LD_LIBRARY_PATH:-}"

    # ==========================================
    # 3. Validate inputs and software
    # ==========================================

    test -s "$SQANTI_DIR/sqanti3_qc.py"
    test -x "$SQANTI_ENV/bin/python"
    test -e "$SQANTI_ENV/lib/libstdc++.so.6"

    test -s "bambu_results/extended_annotations.gtf"
    test -s "!{reference_gtf}"
    test -s "!{reference_fasta}"

    echo "===== SQANTI3 ENVIRONMENT ====="

    "$SQANTI_ENV/bin/python" \
        "$SQANTI_DIR/sqanti3_qc.py" \
        --version

    echo "===== SCIPY DEPENDENCY CHECK ====="

    "$SQANTI_ENV/bin/python" -c \
        'from scipy import stats, optimize; print("SciPy dependencies OK")'

    # ==========================================
    # 4. Run SQANTI3 with ORF prediction
    # ==========================================

    echo "===== SQANTI3 QC + ORF ====="

    mkdir -p sqanti3_results

    "$SQANTI_ENV/bin/python" \
        "$SQANTI_DIR/sqanti3_qc.py" \
        --isoforms bambu_results/extended_annotations.gtf \
        --refGTF "!{reference_gtf}" \
        --refFasta "!{reference_fasta}" \
        --include_ORF \
        --output bambu_sqanti3 \
        --dir sqanti3_results \
        --cpus !{task.cpus} \
        --chunks 1 \
        --report html

    # ==========================================
    # 5. Validate essential outputs
    # ==========================================

    echo "===== VALIDATING OUTPUTS ====="

    test -s sqanti3_results/bambu_sqanti3_classification.txt
    test -s sqanti3_results/bambu_sqanti3_junctions.txt

    echo "SQANTI3 structural classification completed."
    echo "SQANTI3 ORF prediction requested."
    echo "===== SQANTI3 QC COMPLETED ====="
    '''
}
