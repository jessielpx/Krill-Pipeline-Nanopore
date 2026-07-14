#!/usr/bin/env python3

import argparse
import re
from pathlib import Path

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Add gene IDs to a transcript count matrix using a GTF file "
            "and sum transcript counts into gene-level counts."
        )
    )
    parser.add_argument(
        "--gtf",
        required=True,
        help="Merged transcriptome GTF file."
    )
    parser.add_argument(
        "--transcript-counts",
        required=True,
        help="Merged transcript count matrix."
    )
    parser.add_argument(
        "--transcript-output",
        required=True,
        help="Output transcript count matrix containing gene_id."
    )
    parser.add_argument(
        "--gene-output",
        required=True,
        help="Output gene-level count matrix."
    )
    return parser.parse_args()


def extract_attribute(attributes: str, key: str):
    match = re.search(rf'{key}\s+"([^"]+)"', attributes)
    return match.group(1) if match else None


def read_transcript_gene_map(gtf_path: str) -> pd.DataFrame:
    records = []

    with open(gtf_path, "r", encoding="utf-8") as gtf:
        for line in gtf:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) != 9:
                continue

            feature = fields[2]
            attributes = fields[8]

            if feature != "transcript":
                continue

            transcript_id = extract_attribute(attributes, "transcript_id")
            gene_id = extract_attribute(attributes, "gene_id")

            if transcript_id and gene_id:
                records.append(
                    {
                        "transcript_id": transcript_id,
                        "gene_id": gene_id
                    }
                )

    mapping = pd.DataFrame(records)

    if mapping.empty:
        raise ValueError(
            f"No transcript_id and gene_id pairs were found in {gtf_path}"
        )

    return mapping.drop_duplicates(subset=["transcript_id"])


def main():
    args = parse_args()

    mapping = read_transcript_gene_map(args.gtf)

    counts = pd.read_csv(args.transcript_counts, sep="\t")

    if "Name" not in counts.columns:
        raise ValueError(
            "Transcript count matrix must contain a column named 'Name'."
        )

    counts = counts.rename(columns={"Name": "transcript_id"})

    annotated = mapping.merge(
        counts,
        on="transcript_id",
        how="right"
    )

    missing_gene_ids = annotated["gene_id"].isna().sum()

    if missing_gene_ids:
        print(
            f"Warning: {missing_gene_ids} transcripts did not have a gene_id "
            "in the GTF and will be excluded from gene-level aggregation."
        )

    annotated.to_csv(
        args.transcript_output,
        sep="\t",
        index=False
    )

    sample_columns = [
        column
        for column in counts.columns
        if column != "transcript_id"
    ]

    gene_counts = (
        annotated
        .dropna(subset=["gene_id"])
        .groupby("gene_id", as_index=False)[sample_columns]
        .sum()
    )

    gene_counts.to_csv(
        args.gene_output,
        sep="\t",
        index=False
    )


if __name__ == "__main__":
    main()
