#!/usr/bin/env python3

"""
Krill Renew 12 - JAFFAL fusion aggregation.

Inputs:
  manifest TSV with sample, csv
  samplesheet CSV/TSV

Outputs:
  all_fusions.tsv
  high_confidence_fusions.tsv
  filtered_fusions.tsv
  recurrent_fusions.tsv
  fusion_qc.tsv

Preserves original JAFFAL columns.
Does not infer functional validation.
"""

import argparse
import csv
import math
from collections import Counter, defaultdict
from pathlib import Path


def read_table(path):
    with open(path, newline="") as handle:
        first = handle.readline()
        handle.seek(0)

        delimiter = "\t" if "\t" in first else ","

        reader = csv.DictReader(
            handle,
            delimiter=delimiter
        )

        return list(reader), reader.fieldnames or []


def write_table(path, rows, columns):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            delimiter="\t",
            extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)


def as_number(value):
    try:
        n = float(value)
        return n if math.isfinite(n) else None
    except (ValueError, TypeError):
        return None


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--manifest", required=True
    )
    parser.add_argument(
        "--samplesheet", required=True
    )
    parser.add_argument(
        "--outdir", required=True
    )
    parser.add_argument(
        "--min-spanning-reads",
        type=int,
        default=2
    )
    parser.add_argument(
        "--min-recurrent-samples",
        type=int,
        default=2
    )

    args = parser.parse_args()

    if (
        args.min_spanning_reads < 0
        or args.min_recurrent_samples < 2
    ):
        raise ValueError("Invalid filtering thresholds")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    metadata, metadata_cols = read_table(
        args.samplesheet
    )

    if not {"sample", "condition"}.issubset(
        metadata_cols
    ):
        raise ValueError(
            "Samplesheet requires sample and condition"
        )

    sample_info = {}

    for row in metadata:
        sample = row["sample"].strip()

        if not sample or sample in sample_info:
            raise ValueError(
                "Invalid/duplicate sample: " + sample
            )

        sample_info[sample] = row

    manifest, manifest_cols = read_table(
        args.manifest
    )

    if not {"sample", "csv"}.issubset(
        manifest_cols
    ):
        raise ValueError(
            "Manifest requires sample and csv"
        )

    manifest_map = {}

    for entry in manifest:
        sample = entry["sample"].strip()

        if sample in manifest_map:
            raise ValueError(
                "Duplicate JAFFAL sample: " + sample
            )

        manifest_map[sample] = Path(entry["csv"])

    if set(manifest_map) != set(sample_info):
        raise ValueError(
            "JAFFAL/sample mismatch. Missing: "
            + str(sorted(set(sample_info) - set(manifest_map)))
            + "; unexpected: "
            + str(sorted(set(manifest_map) - set(sample_info)))
        )

    records = []
    jaffa_columns = []
    sample_qc = []

    required = {
        "fusion genes",
        "classification",
        "spanning reads",
    }

    for sample in sample_info:

        path = manifest_map[sample]

        if not path.is_file():
            raise FileNotFoundError(path)

        rows, columns = read_table(path)

        if not required.issubset(columns):
            raise ValueError(
                f"Missing JAFFAL columns in {path}: "
                f"{sorted(required - set(columns))}"
            )

        for column in columns:
            if column not in jaffa_columns:
                jaffa_columns.append(column)

        sample_total = 0
        sample_high = 0
        sample_filtered = 0

        for row in rows:

            fusion = row["fusion genes"].strip()
            classification = row[
                "classification"
            ].strip()

            reads = as_number(
                row["spanning reads"]
            )

            if not fusion:
                continue

            high = (
                classification == "HighConfidence"
            )

            read_pass = (
                reads is not None
                and reads >= args.min_spanning_reads
            )

            filtered = high and read_pass

            record = {
                "sample_name": sample,
                "condition": sample_info[sample][
                    "condition"
                ],
                **row,
                "is_high_confidence": high,
                "read_support_pass": read_pass,
                "fusion_filter_pass": filtered,
            }

            records.append(record)
            sample_total += 1
            sample_high += int(high)
            sample_filtered += int(filtered)

        sample_qc.append({
            "sample": sample,
            "condition": sample_info[sample][
                "condition"
            ],
            "all_calls": sample_total,
            "high_confidence_calls": sample_high,
            "filtered_calls": sample_filtered,
        })

    # Preserve the order of the original JAFFAL fields.
    fields = (
        ["sample_name", "condition"]
        + [
            c for c in jaffa_columns
            if c not in ("sample_name", "condition")
        ]
        + [
            "is_high_confidence",
            "read_support_pass",
            "fusion_filter_pass",
        ]
    )

    high_confidence = [
        r for r in records
        if r["is_high_confidence"]
    ]

    filtered = [
        r for r in records
        if r["fusion_filter_pass"]
    ]

    write_table(
        outdir / "all_fusions.tsv",
        records,
        fields
    )

    write_table(
        outdir / "high_confidence_fusions.tsv",
        high_confidence,
        fields
    )

    write_table(
        outdir / "filtered_fusions.tsv",
        filtered,
        fields
    )

    # Recurrence refers to distinct samples,
    # not the number of JAFFAL rows or read counts.
    # Use reported gene pair identity; retain
    # partner ordering rather than silently reversing.
    by_fusion = defaultdict(list)

    for record in filtered:
        by_fusion[
            record["fusion genes"].strip()
        ].append(record)

    recurrent = []

    for fusion, calls in sorted(by_fusion.items()):

        samples = sorted(set(
            r["sample_name"] for r in calls
        ))

        if len(samples) < args.min_recurrent_samples:
            continue

        recurrent.append({
            "fusion_genes": fusion,
            "n_samples": len(samples),
            "samples": ";".join(samples),
            "conditions": ";".join(sorted(set(
                r["condition"] for r in calls
            ))),
            "n_calls": len(calls),
        })

    write_table(
        outdir / "recurrent_fusions.tsv",
        recurrent,
        [
            "fusion_genes",
            "n_samples",
            "samples",
            "conditions",
            "n_calls",
        ]
    )

    write_table(
        outdir / "fusion_sample_qc.tsv",
        sample_qc,
        [
            "sample",
            "condition",
            "all_calls",
            "high_confidence_calls",
            "filtered_calls",
        ]
    )

    metrics = [
        ("samples", len(sample_info)),
        ("all_fusion_calls", len(records)),
        (
            "high_confidence_calls",
            len(high_confidence)
        ),
        ("filtered_calls", len(filtered)),
        (
            "unique_filtered_gene_pairs",
            len(by_fusion)
        ),
        (
            "recurrent_gene_pairs",
            len(recurrent)
        ),
        (
            "min_spanning_reads",
            args.min_spanning_reads
        ),
        (
            "min_recurrent_samples",
            args.min_recurrent_samples
        ),
    ]

    write_table(
        outdir / "fusion_qc.tsv",
        [
            {"metric": key, "value": value}
            for key, value in metrics
        ],
        ["metric", "value"]
    )

    print("===== FUSION QC =====")

    for key, value in metrics:
        print(f"{key}\t{value}")

    print(
        "Note: JAFFAL high-confidence calls "
        "are candidates, not validated fusions."
    )


if __name__ == "__main__":
    main()
