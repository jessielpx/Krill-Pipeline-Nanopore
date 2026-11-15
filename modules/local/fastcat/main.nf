/*
 * Krill Renew - FASTQ preparation
 *
 * Merge FASTQ files belonging to one sample.
 * Generate basic read QC.
 *
 * Native execution on Rorqual.
 * No FASTCAT binary or container required.
 */

process FASTCAT {

    tag "${meta.sample}"

    cpus 2
    memory '8 GB'
    time '8h'

    publishDir "${params.outdir}/01_preprocessing/fastcat/${meta.sample}",
        mode: 'copy'

    input:
    tuple val(meta), path(read_dir)

    output:
    tuple val(meta),
          path("merged_fastq/${meta.sample}.fastq.gz"),
          path("fastcat_stats"),
          emit: results

    script:

    def sample = meta.sample.toString()

    if (!(sample ==~ /[A-Za-z0-9._-]+/)) {
        error "Invalid sample ID: ${sample}"
    }

    """
    set -euo pipefail

    mkdir -p merged_fastq fastcat_stats

    env -u PYTHONPATH -u PYTHONHOME python3 - \\
        "${sample}" "${read_dir}" <<'PY'

import csv
import gzip
import sys
from pathlib import Path

sample = sys.argv[1]
read_dir = Path(sys.argv[2])

extensions = (
    ".fastq.gz",
    ".fq.gz",
    ".fastq",
    ".fq",
)

fastqs = sorted(
    p for p in read_dir.iterdir()
    if p.is_file() and p.name.endswith(extensions)
)

if not fastqs:
    raise RuntimeError(
        f"No FASTQ files found for {sample}"
    )

output = Path("merged_fastq") / f"{sample}.fastq.gz"
stats_path = Path("fastcat_stats/per-file-stats.tsv")
validation_path = Path("fastcat_stats/validation.tsv")

total_reads = 0

with gzip.open(output, "wb", compresslevel=1) as merged:
    with stats_path.open("w", newline="") as stats_file:

        writer = csv.writer(
            stats_file,
            delimiter="\\t"
        )

        writer.writerow([
            "filename",
            "sample_name",
            "n_seqs",
            "n_bases",
            "min_length",
            "max_length",
            "mean_quality",
        ])

        for fastq in fastqs:

            opener = (
                gzip.open
                if fastq.name.endswith(".gz")
                else open
            )

            n_reads = 0
            n_bases = 0
            min_length = None
            max_length = 0
            quality_sum = 0

            with opener(fastq, "rb") as source:

                while True:
                    header = source.readline()

                    if not header:
                        break

                    sequence = source.readline()
                    plus = source.readline()
                    quality = source.readline()

                    if (
                        not header.startswith(b"@")
                        or not plus.startswith(b"+")
                        or not sequence
                        or not quality
                    ):
                        raise RuntimeError(
                            f"Invalid FASTQ record: {fastq}"
                        )

                    seq = sequence.rstrip(b"\\r\\n")
                    qual = quality.rstrip(b"\\r\\n")

                    if len(seq) != len(qual):
                        raise RuntimeError(
                            f"Sequence/quality mismatch: {fastq}"
                        )

                    length = len(seq)

                    n_reads += 1
                    n_bases += length
                    quality_sum += sum(qual) - 33 * length

                    if min_length is None:
                        min_length = length
                    else:
                        min_length = min(min_length, length)

                    max_length = max(max_length, length)

                    merged.write(header)
                    merged.write(sequence)
                    merged.write(plus)
                    merged.write(quality)

            if n_reads == 0:
                raise RuntimeError(
                    f"Empty FASTQ file: {fastq}"
                )

            mean_quality = (
                quality_sum / n_bases
                if n_bases else 0
            )

            writer.writerow([
                str(fastq),
                sample,
                n_reads,
                n_bases,
                min_length,
                max_length,
                f"{mean_quality:.3f}",
            ])

            total_reads += n_reads

with validation_path.open("w") as handle:
    handle.write(
        "sample\\tinput_files\\ttotal_reads\\n"
    )
    handle.write(
        f"{sample}\\t{len(fastqs)}\\t{total_reads}\\n"
    )

print(f"Sample: {sample}")
print(f"Input FASTQ files: {len(fastqs)}")
print(f"Total reads: {total_reads}")
print(f"Output: {output}")

PY

    test -s "merged_fastq/${sample}.fastq.gz"
    gzip -t "merged_fastq/${sample}.fastq.gz"

    echo "FASTQ preparation completed: ${sample}"
    """
}
