#!/usr/bin/env bash

set -euo pipefail

###############################################################################
# Download and prepare the JAFFAL hg38 + GENCODE v49 reference.
#
# Usage:
#   bash scripts/download_jaffal_reference.sh /path/to/install/directory
#
# Example:
#   bash scripts/download_jaffal_reference.sh \
#     /project/def-lefranco/peixiliu/references/JAFFAL
###############################################################################

INSTALL_ROOT="${1:-}"

if [[ -z "${INSTALL_ROOT}" ]]; then
    echo "ERROR: Please provide an installation directory."
    echo
    echo "Usage:"
    echo "  bash scripts/download_jaffal_reference.sh /path/to/install/directory"
    exit 1
fi

for command_name in curl python3 tar find; do
    if ! command -v "${command_name}" >/dev/null 2>&1; then
        echo "ERROR: Required command not found: ${command_name}"
        exit 1
    fi
done

mkdir -p "${INSTALL_ROOT}"
INSTALL_ROOT="$(cd "${INSTALL_ROOT}" && pwd)"

DOWNLOAD_DIR="${INSTALL_ROOT}/downloads"
REFERENCE_DIR="${INSTALL_ROOT}/hg38_gencode49"

mkdir -p "${DOWNLOAD_DIR}"
mkdir -p "${REFERENCE_DIR}"

echo "============================================================"
echo "JAFFAL reference setup"
echo "============================================================"
echo "Installation root: ${INSTALL_ROOT}"
echo

###############################################################################
# Find the reference archive from the latest JAFFA GitHub release
###############################################################################

echo "Checking the latest JAFFA GitHub release..."

ASSET_INFO="$(
python3 <<'PY'
import json
import sys
import urllib.request

api_url = "https://api.github.com/repos/Oshlack/JAFFA/releases/latest"

request = urllib.request.Request(
    api_url,
    headers={
        "Accept": "application/vnd.github+json",
        "User-Agent": "Krill-Pipeline-Nanopore"
    }
)

try:
    with urllib.request.urlopen(request) as response:
        release = json.load(response)
except Exception as exc:
    print(f"ERROR\tCould not query GitHub: {exc}")
    sys.exit(1)

assets = release.get("assets", [])

matches = []

for asset in assets:
    name = asset.get("name", "")
    lower = name.lower()

    is_archive = lower.endswith((
        ".tar.gz",
        ".tgz",
        ".tar",
        ".zip"
    ))

    is_reference = (
        "hg38" in lower
        and (
            "gencode49" in lower
            or "gencode_49" in lower
            or "gencode-v49" in lower
            or "gencodev49" in lower
        )
    )

    if is_archive and is_reference:
        matches.append({
            "name": name,
            "url": asset.get("browser_download_url", "")
        })

if not matches:
    print("ERROR\tNo hg38 + GENCODE v49 archive was found.")
    print("AVAILABLE\t" + "|".join(
        asset.get("name", "unnamed")
        for asset in assets
    ))
    sys.exit(1)

selected = matches[0]

print(f"OK\t{selected['name']}\t{selected['url']}")
PY
)"

if [[ "${ASSET_INFO}" == ERROR* ]]; then
    echo "${ASSET_INFO}" | tr '\t' ' '
    exit 1
fi

if [[ "${ASSET_INFO}" == AVAILABLE* ]]; then
    echo "${ASSET_INFO}" | tr '|' '\n'
    exit 1
fi

ASSET_NAME="$(printf '%s\n' "${ASSET_INFO}" | cut -f2)"
ASSET_URL="$(printf '%s\n' "${ASSET_INFO}" | cut -f3-)"

if [[ -z "${ASSET_NAME}" || -z "${ASSET_URL}" ]]; then
    echo "ERROR: Could not determine the reference archive."
    exit 1
fi

ARCHIVE_PATH="${DOWNLOAD_DIR}/${ASSET_NAME}"

echo "Reference archive: ${ASSET_NAME}"
echo

###############################################################################
# Download
###############################################################################

if [[ -s "${ARCHIVE_PATH}" ]]; then
    echo "Archive already exists. Download skipped:"
    echo "  ${ARCHIVE_PATH}"
else
    echo "Downloading reference archive..."

    curl \
        --fail \
        --location \
        --retry 5 \
        --retry-delay 10 \
        --continue-at - \
        "${ASSET_URL}" \
        --output "${ARCHIVE_PATH}"
fi

echo

###############################################################################
# Extract
###############################################################################

MARKER_FILE="${REFERENCE_DIR}/.extraction_complete"

if [[ -f "${MARKER_FILE}" ]]; then
    echo "Reference has already been extracted. Extraction skipped."
else
    echo "Extracting reference archive..."

    case "${ARCHIVE_PATH}" in
        *.tar.gz|*.tgz)
            tar -xzf "${ARCHIVE_PATH}" -C "${REFERENCE_DIR}"
            ;;

        *.tar)
            tar -xf "${ARCHIVE_PATH}" -C "${REFERENCE_DIR}"
            ;;

        *.zip)
            if ! command -v unzip >/dev/null 2>&1; then
                echo "ERROR: unzip is required for this archive."
                exit 1
            fi

            unzip -q "${ARCHIVE_PATH}" -d "${REFERENCE_DIR}"
            ;;

        *)
            echo "ERROR: Unsupported archive type:"
            echo "  ${ARCHIVE_PATH}"
            exit 1
            ;;
    esac

    touch "${MARKER_FILE}"
fi

echo

###############################################################################
# Find required JAFFAL files
###############################################################################

find_first() {
    find "${REFERENCE_DIR}" \
        -type f \
        "$@" \
        -print \
        -quit
}

GENOME_FASTA="$(
    find_first -iname "hg38.fa"
)"

TRANSCRIPTOME_FASTA="$(
    find_first \
        \( \
            -iname "hg38_gencode49.fa" \
            -o -iname "hg38_gencode_49.fa" \
            -o -iname "*gencode*49*.fa" \
        \)
)"

ANNOTATION_TAB="$(
    find_first \
        \( \
            -iname "hg38_gencode49.tab" \
            -o -iname "hg38_gencode_49.tab" \
            -o -iname "*gencode*49*.tab" \
        \)
)"

ANNOTATION_BED="$(
    find_first \
        \( \
            -iname "hg38_gencode49.bed" \
            -o -iname "hg38_gencode_49.bed" \
            -o -iname "*gencode*49*.bed" \
        \)
)"

check_file() {
    local label="$1"
    local path="$2"

    if [[ -n "${path}" && -s "${path}" ]]; then
        echo "Found ${label}:"
        echo "  ${path}"
    else
        echo "ERROR: Missing ${label}"
        return 1
    fi
}

echo "Checking required files..."

FAILED=0

check_file "genome FASTA" "${GENOME_FASTA}" || FAILED=1
check_file "transcriptome FASTA" "${TRANSCRIPTOME_FASTA}" || FAILED=1
check_file "annotation TAB" "${ANNOTATION_TAB}" || FAILED=1
check_file "annotation BED" "${ANNOTATION_BED}" || FAILED=1

if [[ "${FAILED}" -ne 0 ]]; then
    echo
    echo "Reference setup did not find all expected files."
    echo "Inspect:"
    echo "  ${REFERENCE_DIR}"
    exit 1
fi

ACTUAL_REFERENCE_DIR="$(dirname "${GENOME_FASTA}")"

###############################################################################
# Check Bowtie2 indexes
###############################################################################

has_bowtie2_index() {
    local prefix="$1"

    [[ -s "${prefix}.1.bt2" ]] \
        || [[ -s "${prefix}.1.bt2l" ]]
}

GENOME_INDEX_PREFIX="${GENOME_FASTA%.fa}"
TRANSCRIPTOME_INDEX_PREFIX="${TRANSCRIPTOME_FASTA%.fa}"

echo
echo "Checking Bowtie2 indexes..."

if has_bowtie2_index "${GENOME_INDEX_PREFIX}"; then
    echo "Genome Bowtie2 index already exists."
else
    echo "Genome Bowtie2 index is missing."

    if command -v bowtie2-build >/dev/null 2>&1; then
        echo "Building genome Bowtie2 index..."

        bowtie2-build \
            --threads "${SLURM_CPUS_PER_TASK:-4}" \
            "${GENOME_FASTA}" \
            "${GENOME_INDEX_PREFIX}"
    else
        echo
        echo "ERROR: bowtie2-build is not available."
        echo
        echo "Load Bowtie2 and rerun the script:"
        echo "  module spider bowtie2"
        echo "  module load bowtie2"
        exit 1
    fi
fi

if has_bowtie2_index "${TRANSCRIPTOME_INDEX_PREFIX}"; then
    echo "Transcriptome Bowtie2 index already exists."
else
    echo "Transcriptome Bowtie2 index is missing."

    if command -v bowtie2-build >/dev/null 2>&1; then
        echo "Building transcriptome Bowtie2 index..."

        bowtie2-build \
            --threads "${SLURM_CPUS_PER_TASK:-4}" \
            "${TRANSCRIPTOME_FASTA}" \
            "${TRANSCRIPTOME_INDEX_PREFIX}"
    else
        echo
        echo "ERROR: bowtie2-build is not available."
        exit 1
    fi
fi

###############################################################################
# Final result
###############################################################################

echo
echo "============================================================"
echo "JAFFAL reference setup completed"
echo "============================================================"
echo
echo "JAFFAL reference directory:"
echo
echo "  ${ACTUAL_REFERENCE_DIR}"
echo
echo
echo "Run the pipeline with:"
echo
echo "  --jaffal_ref_dir ${ACTUAL_REFERENCE_DIR}"
echo
echo "or add this to nextflow.config:"
echo
echo "  params.jaffal_ref_dir = '${ACTUAL_REFERENCE_DIR}'"
echo
