/*
 * Krill Renew 12 - JAFFAL
 *
 * Native execution on Rorqual.
 *
 * Each task creates an isolated JAFFA runtime:
 *   - original JAFFA installation remains unchanged
 *   - existing GENCODE49 reference is reused
 *   - JAFFA_stages.groovy is patched locally
 *
 * Input:
 *   tuple(meta, merged FASTQ)
 *
 * Output:
 *   tuple(meta, sample.jaffa_results.csv)
 */

process JAFFAL {

    tag "${meta.sample}"

    cpus 16
    memory '64 GB'
    time '24h'

    publishDir "${params.outdir}/12_fusion_detection/per_sample/${meta.sample}",
        mode: 'copy'

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta),
          path("${meta.sample}.jaffa_results.csv"),
          emit: calls

    shell:
    '''
    set -euo pipefail

    module --force purge
    module load StdEnv/2023
    module load java/17.0.21
    module load python/3.11.5
    module load minimap2
    module load samtools

    unset PYTHONHOME PYTHONPATH
    hash -r

    JAFFA_HOME="!{params.jaffal_home ?: ''}"

    JAFFA_REF="!{params.jaffal_ref_dir ?: ''}"

    if [[ -z "$JAFFA_HOME" || -z "$JAFFA_REF" ]]; then
        echo "ERROR: Configure jaffal_home and jaffal_ref_dir" >&2
        exit 1
    fi

    SAMPLE="!{meta.sample}"
    FASTQ="!{fastq}"

    BPIPE="$JAFFA_HOME/tools/bin/bpipe"

    [[ "$SAMPLE" =~ ^[A-Za-z0-9._-]+$ ]] || {
        echo "ERROR: Invalid sample ID: $SAMPLE" >&2
        exit 1
    }

    echo "===== JAFFAL CONFIGURATION ====="
    echo "Sample: $SAMPLE"
    echo "JAFFA home: $JAFFA_HOME"
    echo "Reference: $JAFFA_REF"
    echo "FASTQ: $FASTQ"

    test -s "$FASTQ"
    test -s "$JAFFA_HOME/JAFFAL.groovy"
    test -s "$JAFFA_HOME/JAFFA_stages.groovy"
    test -s "$JAFFA_HOME/tools.groovy"
    test -x "$BPIPE"

    # Check the essential existing reference files.
    for file in \
        "$JAFFA_REF/hg38.fa" \
        "$JAFFA_REF/hg38_gencode49.fa" \
        "$JAFFA_REF/hg38_gencode49.tab" \
        "$JAFFA_REF/hg38_gencode49.bed" \
        "$JAFFA_REF/hg38_gencode49.1.bt2" \
        "$JAFFA_REF/Masked_hg38.1.bt2"
    do
        test -s "$file" || {
            echo "ERROR: Missing JAFFA reference: $file" >&2
            exit 1
        }
    done

    # Create an isolated copy of only the configuration
    # and entrypoint. Other JAFFA assets are symlinked.
    #
    # The installed JAFFA directory is never modified.

    python3 - "$JAFFA_HOME" "$JAFFA_REF" <<'PY'
import json
import re
import shutil
import sys
from pathlib import Path

home = Path(sys.argv[1]).resolve()
reference = Path(sys.argv[2]).resolve()

runtime = Path("jaffa_runtime")
runtime.mkdir(exist_ok=True)

for source in home.iterdir():

    if source.name in {
        "JAFFA_stages.groovy",
        "JAFFAL.groovy",
    }:
        continue

    target = runtime / source.name

    if not target.exists() and not target.is_symlink():
        target.symlink_to(source)

shutil.copy2(
    home / "JAFFAL.groovy",
    runtime / "JAFFAL.groovy"
)

source_text = (
    home / "JAFFA_stages.groovy"
).read_text()

pattern = r'(?m)^refBase\s*=\s*codeBase\s*$'

replacement = (
    "refBase = "
    + json.dumps(str(reference))
)

patched, count = re.subn(
    pattern,
    lambda _: replacement,
    source_text
)

if count != 1:
    raise SystemExit(
        "ERROR: Expected exactly one "
        "refBase = codeBase assignment"
    )

(runtime / "JAFFA_stages.groovy").write_text(
    patched
)

print("===== LOCAL JAFFA RUNTIME =====")
print("Runtime:", runtime.resolve())
print("Reference:", reference)
print("Original installation unchanged")
PY

    test -s jaffa_runtime/JAFFAL.groovy
    test -s jaffa_runtime/JAFFA_stages.groovy

    echo "===== JAFFAL RUN ====="

    "$BPIPE" run \
        -p threads=!{task.cpus} \
        jaffa_runtime/JAFFAL.groovy \
        "$FASTQ"

    # Validate the completed JAFFAL output.
    if [[ ! -f jaffa_results.csv ]]; then
        echo "ERROR: JAFFAL did not produce jaffa_results.csv" >&2
        exit 1
    fi

    # A header-only CSV is valid when no fusion is found.
    python3 - jaffa_results.csv <<'PY'
import csv
import sys

with open(sys.argv[1], newline="") as f:
    header = csv.DictReader(f).fieldnames or []

required = {
    "fusion genes",
    "classification",
    "spanning reads",
}

missing = required - set(header)

if missing:
    raise SystemExit(
        "ERROR: Invalid JAFFAL CSV columns: "
        + str(sorted(missing))
    )

print("JAFFAL output schema validated")
PY

    cp jaffa_results.csv \
       "$SAMPLE.jaffa_results.csv"

    echo "===== JAFFAL FINISHED ====="
    echo "Sample: $SAMPLE"
    '''
}
