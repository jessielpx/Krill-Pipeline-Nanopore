import os
#!/usr/bin/env python3

"""
Krill Renew - Unified Core QC

Combine:
  FASTQ read statistics
  Genome alignment statistics
  bambu transcript quantification QC
  SQANTI3 structural classification

No third-party Python packages required.
"""

import argparse
import csv
import re
from collections import Counter
from pathlib import Path


def read_table(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def require_file(path):
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Missing or empty file: {path}")
    return path


def collect_fastq_qc(directory):
    results = {}

    validation_paths = []

    for root, dirs, files in os.walk(directory, followlinks=True):
        dirs.sort()
        if "validation.tsv" in files:
            validation_paths.append(Path(root) / "validation.tsv")

    for path in sorted(validation_paths):
        for row in read_table(path):
            sample = row["sample"]

            if sample in results:
                raise ValueError(f"Duplicate FASTQ QC: {sample}")

            results[sample] = {
                "fastq_reads": int(row["total_reads"]),
                "fastq_input_files": int(row["input_files"]),
            }

    if not results:
        raise ValueError("No FASTQ validation.tsv files found")

    return results


def collect_alignment_qc(directory):
    results = {}

    for path in sorted(Path(directory).rglob("*.flagstat.txt")):
        sample = path.name.removesuffix(".flagstat.txt")

        primary = None
        primary_mapped = None

        for line in path.read_text().splitlines():

            match = re.match(
                r"^\s*(\d+)\s+\+\s+\d+\s+primary(?:\s|$)",
                line,
            )
            if match:
                primary = int(match.group(1))

            match = re.match(
                r"^\s*(\d+)\s+\+\s+\d+\s+primary mapped(?:\s|$)",
                line,
            )
            if match:
                primary_mapped = int(match.group(1))

        if primary is None or primary_mapped is None:
            raise ValueError(f"Cannot parse flagstat: {path}")

        if sample in results:
            raise ValueError(f"Duplicate alignment QC: {sample}")

        results[sample] = {
            "primary_alignments": primary,
            "primary_mapped": primary_mapped,
        }

    if not results:
        raise ValueError("No samtools flagstat files found")

    return results


def collect_bambu_qc(directory):
    path = require_file(
        Path(directory) / "bambu_sample_qc.tsv"
    )

    results = {}

    for row in read_table(path):
        sample = row["sample"]

        if sample in results:
            raise ValueError(f"Duplicate bambu sample: {sample}")

        results[sample] = {
            "total_transcript_features": int(
                row["reference_transcripts"]
            ),
            "detected_transcripts": int(
                row["detected_transcripts"]
            ),
            "transcript_estimated_counts": round(
                float(row["transcript_estimated_counts"]), 2
            ),
        }

    return results


def collect_sqanti3_qc(directory):
    path = require_file(
        Path(directory) / "bambu_sqanti3_classification.txt"
    )

    categories = Counter()

    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        if "structural_category" not in (reader.fieldnames or []):
            raise ValueError(
                "SQANTI3 classification lacks structural_category"
            )

        for row in reader:
            category = row["structural_category"].strip()
            categories[category or "Unclassified"] += 1

    if not categories:
        raise ValueError("SQANTI3 classification is empty")

    return categories


def write_tsv(path, rows, columns):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--fastcat", required=True)
    parser.add_argument("--alignment", required=True)
    parser.add_argument("--bambu", required=True)
    parser.add_argument("--sqanti", required=True)
    parser.add_argument("--out", required=True)

    args = parser.parse_args()

    output = Path(args.out)
    output.mkdir(parents=True, exist_ok=True)

    fastq = collect_fastq_qc(args.fastcat)
    alignment = collect_alignment_qc(args.alignment)
    bambu = collect_bambu_qc(args.bambu)
    categories = collect_sqanti3_qc(args.sqanti)

    samples = sorted(bambu)

    for source_name, source in [
        ("FASTQ", fastq),
        ("alignment", alignment),
    ]:
        if set(source) != set(samples):
            raise ValueError(
                f"{source_name} sample IDs do not match bambu: "
                f"missing={sorted(set(samples) - set(source))}, "
                f"extra={sorted(set(source) - set(samples))}"
            )

    sample_rows = []

    for sample in samples:
        fastq_reads = fastq[sample]["fastq_reads"]
        primary_mapped = alignment[sample]["primary_mapped"]

        if fastq_reads <= 0:
            raise ValueError(f"No FASTQ reads for {sample}")

        sample_rows.append({
            "sample": sample,
            **fastq[sample],
            **alignment[sample],
            **bambu[sample],
            "primary_mapping_pct": round(
                primary_mapped / fastq_reads * 100, 2
            ),
        })

    sample_columns = [
        "sample",
        "fastq_reads",
        "fastq_input_files",
        "primary_alignments",
        "primary_mapped",
        "primary_mapping_pct",
        "total_transcript_features",
        "detected_transcripts",
        "transcript_estimated_counts",
    ]

    write_tsv(
        output / "sample_qc_summary.tsv",
        sample_rows,
        sample_columns,
    )

    total_classified = sum(categories.values())

    category_rows = [
        {
            "structural_category": category,
            "transcripts": count,
            "percentage": round(
                count / total_classified * 100, 2
            ),
        }
        for category, count in categories.most_common()
    ]

    write_tsv(
        output / "sqanti3_category_summary.tsv",
        category_rows,
        [
            "structural_category",
            "transcripts",
            "percentage",
        ],
    )

    with open(output / "qc_summary.txt", "w") as handle:
        handle.write("KRILL RENEW - UNIFIED CORE QC\n")
        handle.write("============================\n\n")
        handle.write(f"Samples: {len(samples)}\n")
        handle.write(
            f"SQANTI3 classified transcripts: {total_classified}\n\n"
        )

        handle.write("Sample-level QC\n")
        for row in sample_rows:
            handle.write(
                f"{row['sample']}: "
                f"FASTQ reads={row['fastq_reads']}, "
                f"primary mapping={row['primary_mapping_pct']}%, "
                f"detected transcripts={row['detected_transcripts']}, "
                f"estimated counts={row['transcript_estimated_counts']}\n"
            )

        handle.write("\nSQANTI3 structural categories\n")
        for row in category_rows:
            handle.write(
                f"{row['structural_category']}: "
                f"{row['transcripts']} "
                f"({row['percentage']}%)\n"
            )

    print("Unified QC completed")
    print(f"Samples: {len(samples)}")
    print(f"SQANTI3 classified transcripts: {total_classified}")
    print(f"Output: {output}")


if __name__ == "__main__":
    main()
