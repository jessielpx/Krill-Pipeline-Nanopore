#!/usr/bin/env python3

import argparse
import re

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Add gene annotations to a transcript count matrix using a GTF "
            "file and sum transcript counts into gene-level counts."
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
        help=(
            "Output transcript count matrix containing transcript_id, "
            "gene_id, ref_gene_id, gene_name, and aggregation_gene_id."
        )
    )

    parser.add_argument(
        "--gene-output",
        required=True,
        help="Output gene-level count matrix."
    )

    return parser.parse_args()


def extract_attribute(attributes: str, key: str):
    match = re.search(
        rf'(?:^|;\s*){re.escape(key)}\s+"([^"]+)"',
        attributes
    )

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

            transcript_id = extract_attribute(
                attributes,
                "transcript_id"
            )

            stringtie_gene_id = extract_attribute(
                attributes,
                "gene_id"
            )

            ref_gene_id = extract_attribute(
                attributes,
                "ref_gene_id"
            )

            gene_name = extract_attribute(
                attributes,
                "gene_name"
            )

            if not transcript_id:
                continue

            # Use the stable reference gene ID whenever available.
            # For transcripts without ref_gene_id, use gene_name if present.
            # Truly novel transcripts fall back to the StringTie MSTRG ID.
            aggregation_gene_id = (
                ref_gene_id
                or gene_name
                or stringtie_gene_id
            )

            if not aggregation_gene_id:
                continue

            records.append(
                {
                    "transcript_id": transcript_id,
                    "stringtie_gene_id": stringtie_gene_id,
                    "ref_gene_id": ref_gene_id,
                    "gene_name": gene_name,
                    "aggregation_gene_id": aggregation_gene_id
                }
            )

    mapping = pd.DataFrame(records)

    if mapping.empty:
        raise ValueError(
            "No usable transcript annotations were found in "
            f"{gtf_path}"
        )

    duplicate_transcripts = mapping[
        mapping.duplicated(
            subset=["transcript_id"],
            keep=False
        )
    ]

    if not duplicate_transcripts.empty:
        conflicting = (
            duplicate_transcripts
            .groupby("transcript_id")["aggregation_gene_id"]
            .nunique()
        )

        conflicting = conflicting[conflicting > 1]

        if not conflicting.empty:
            examples = ", ".join(
                conflicting.index.astype(str)[:10]
            )

            raise ValueError(
                "Some transcript IDs map to more than one gene. "
                f"Examples: {examples}"
            )

    mapping = mapping.drop_duplicates(
        subset=["transcript_id"],
        keep="first"
    )

    return mapping


def main():
    args = parse_args()

    mapping = read_transcript_gene_map(args.gtf)

    counts = pd.read_csv(
        args.transcript_counts,
        sep="\t"
    )

    if "Name" not in counts.columns:
        raise ValueError(
            "Transcript count matrix must contain a column named 'Name'."
        )

    counts = counts.rename(
        columns={"Name": "transcript_id"}
    )

    if counts["transcript_id"].duplicated().any():
        duplicated_ids = (
            counts.loc[
                counts["transcript_id"].duplicated(
                    keep=False
                ),
                "transcript_id"
            ]
            .drop_duplicates()
            .astype(str)
            .tolist()
        )

        examples = ", ".join(duplicated_ids[:10])

        raise ValueError(
            "Transcript count matrix contains duplicated transcript IDs. "
            f"Examples: {examples}"
        )

    sample_columns = [
        column
        for column in counts.columns
        if column != "transcript_id"
    ]

    if not sample_columns:
        raise ValueError(
            "No sample count columns were found in the transcript "
            "count matrix."
        )

    for column in sample_columns:
        counts[column] = pd.to_numeric(
            counts[column],
            errors="raise"
        )

    annotated = mapping.merge(
        counts,
        on="transcript_id",
        how="right",
        validate="one_to_one"
    )

    missing_gene_ids = (
        annotated["aggregation_gene_id"]
        .isna()
        .sum()
    )

    if missing_gene_ids:
        print(
            f"Warning: {missing_gene_ids} transcripts did not have a "
            "matching gene annotation in the GTF and will be excluded "
            "from gene-level aggregation."
        )

    annotated = annotated[
        [
            "transcript_id",
            "stringtie_gene_id",
            "ref_gene_id",
            "gene_name",
            "aggregation_gene_id",
            *sample_columns
        ]
    ]

    annotated.to_csv(
        args.transcript_output,
        sep="\t",
        index=False
    )

    gene_counts = (
        annotated
        .dropna(subset=["aggregation_gene_id"])
        .groupby(
            "aggregation_gene_id",
            as_index=False
        )[sample_columns]
        .sum()
        .rename(
            columns={
                "aggregation_gene_id": "gene_id"
            }
        )
    )

    gene_counts.to_csv(
        args.gene_output,
        sep="\t",
        index=False
    )

    reference_gene_count = (
        gene_counts["gene_id"]
        .astype(str)
        .str.startswith("ENSG")
        .sum()
    )

    novel_gene_count = (
        gene_counts["gene_id"]
        .astype(str)
        .str.startswith("MSTRG.")
        .sum()
    )

    print(
        f"Transcript annotations read: {len(mapping)}"
    )

    print(
        f"Transcripts in count matrix: {len(counts)}"
    )

    print(
        f"Gene-level rows written: {len(gene_counts)}"
    )

    print(
        f"Reference ENSG rows: {reference_gene_count}"
    )

    print(
        f"Novel MSTRG rows: {novel_gene_count}"
    )


if __name__ == "__main__":
    main()
