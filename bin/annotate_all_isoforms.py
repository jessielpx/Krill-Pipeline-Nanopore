#!/usr/bin/env python3

"""
Krill Renew - reference-based annotation of ALL Bambu isoforms.

Input:
    Bambu extended_annotations.gtf
    Ensembl reference GTF

Output:
    isoform_structure_annotations.tsv
    isoform_structure_qc.tsv

No HR/control comparison.
No gene expression filter.
No ORF prediction.
"""

import argparse
import csv
import gzip
import re
from collections import Counter, defaultdict
from pathlib import Path


ATTR = re.compile(
    r'([A-Za-z_][A-Za-z0-9_]*) "([^"]*)"'
)

KINDS = {
    "transcript",
    "exon",
    "CDS",
    "five_prime_utr",
    "three_prime_utr",
    "stop_codon",
}


def open_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "rt", encoding="utf-8")


def parse_attributes(value):
    return dict(ATTR.findall(value))


def new_record():
    return {
        "gene_id": "",
        "gene_name": "",
        "biotype": "",
        "protein_id": "",
        "chrom": "",
        "strand": "",
        "exon": [],
        "CDS": [],
        "five_prime_utr": [],
        "three_prime_utr": [],
        "stop_codon": [],
    }


def load_gtf(path, wanted=None):
    records = defaultdict(new_record)

    with open_text(path) as handle:
        for line in handle:

            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 9:
                continue

            kind = fields[2]

            if kind not in KINDS:
                continue

            attributes = parse_attributes(fields[8])
            transcript = attributes.get("transcript_id")

            if not transcript:
                continue

            if wanted is not None and transcript not in wanted:
                continue

            record = records[transcript]

            for key in ("gene_id", "gene_name", "protein_id"):
                value = attributes.get(key)
                if value:
                    record[key] = value

            biotype = (
                attributes.get("transcript_biotype")
                or attributes.get("transcript_type")
            )

            if biotype:
                record["biotype"] = biotype

            record["chrom"] = fields[0]
            record["strand"] = fields[6]

            if kind == "transcript":
                continue

            start = int(fields[3])
            end = int(fields[4])

            if start > end:
                raise ValueError(
                    f"Invalid coordinates: {transcript}"
                )

            phase = (
                fields[7] if kind == "CDS" else "."
            )

            record[kind].append((start, end, phase))

    return dict(records)


def coordinates(features, include_phase=False):
    if include_phase:
        return tuple(sorted(features))

    return tuple(
        sorted(
            (start, end)
            for start, end, _ in features
        )
    )


def total_length(features):
    return sum(
        end - start + 1
        for start, end, _ in features
    )


def write_tsv(path, fieldnames, rows):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--bambu-gtf",
        required=True
    )
    parser.add_argument(
        "--reference-gtf",
        required=True
    )
    parser.add_argument(
        "--outdir",
        required=True
    )

    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print("Loading Bambu transcript annotations", flush=True)

    bambu = load_gtf(args.bambu_gtf)

    wanted = set(bambu)

    if not wanted:
        raise ValueError(
            "No transcripts found in Bambu GTF"
        )

    print(
        "Bambu transcript IDs:",
        len(wanted),
        flush=True
    )

    print("Loading reference annotations", flush=True)

    reference = load_gtf(
        args.reference_gtf,
        wanted=wanted
    )

    columns = [
        "gene_id",
        "transcript_id",
        "gene_name",
        "biotype",
        "annotation_status",
        "bambu_ref_exon_match",
        "chrom",
        "strand",
        "bambu_exons",
        "reference_exons",
        "reference_spliced_length",
        "cds_segments",
        "cds_length_excluding_stop",
        "stop_codon_length",
        "utr5_length",
        "utr3_length",
        "protein_id",
        "comparison_to_top",
        "same_annotated_cds_as_top",
    ]

    counts = Counter()
    genes = set()
    mismatches = 0

    output_file = (
        outdir / "isoform_structure_annotations.tsv"
    )

    with open(output_file, "w", newline="") as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            delimiter="\t"
        )
        writer.writeheader()

        for transcript in sorted(bambu):

            b = bambu[transcript]
            r = reference.get(transcript)

            gene_id = b["gene_id"] or (
                r["gene_id"] if r else ""
            )

            if not gene_id:
                raise ValueError(
                    "Missing gene ID for " + transcript
                )

            genes.add(gene_id)

            b_exons = coordinates(b["exon"])

            r_exons = (
                coordinates(r["exon"])
                if r else ()
            )

            # Compare chromosome, strand,
            # and the complete exon chain.
            exon_match = bool(
                b_exons
                and r_exons
                and b_exons == r_exons
                and b["chrom"] == r["chrom"]
                and b["strand"] == r["strand"]
            )

            if b_exons and r_exons and not exon_match:
                mismatches += 1

            if r is None:
                status = "not_in_reference"
            elif not b_exons:
                status = "no_bambu_exons"
            elif not exon_match:
                status = "exon_chain_mismatch"
            elif not r["CDS"]:
                status = "matched_no_CDS"
            else:
                status = "matched_with_CDS"

            counts[status] += 1

            ref_ok = bool(r and exon_match)

            cds = r["CDS"] if ref_ok else []

            row = {
                "gene_id": gene_id,
                "transcript_id": transcript,
                "gene_name": (
                    r["gene_name"] if r and r["gene_name"]
                    else b["gene_name"]
                ),
                "biotype": (
                    r["biotype"] if r else b["biotype"]
                ),
                "annotation_status": status,
                "bambu_ref_exon_match": (
                    "TRUE" if exon_match else "FALSE"
                ),
                "chrom": b["chrom"],
                "strand": b["strand"],
                "bambu_exons": len(b["exon"]),
                "reference_exons": (
                    len(r["exon"]) if r else 0
                ),
                "reference_spliced_length": (
                    total_length(r["exon"])
                    if ref_ok else "NA"
                ),
                "cds_segments": len(cds),
                "cds_length_excluding_stop": (
                    total_length(cds)
                    if cds else "NA"
                ),
                "stop_codon_length": (
                    total_length(r["stop_codon"])
                    if ref_ok and r["stop_codon"]
                    else "NA"
                ),
                "utr5_length": (
                    total_length(r["five_prime_utr"])
                    if ref_ok and r["five_prime_utr"]
                    else "NA"
                ),
                "utr3_length": (
                    total_length(r["three_prime_utr"])
                    if ref_ok and r["three_prime_utr"]
                    else "NA"
                ),
                "protein_id": (
                    r["protein_id"] if ref_ok else ""
                ),
                # There is no differential-usage
                # top transcript in all-isoform mode.
                "comparison_to_top": "NA",
                "same_annotated_cds_as_top": "NA",
            }

            writer.writerow(row)

    metrics = [
        ("mode", "all_transcripts"),
        ("candidate_genes", len(genes)),
        ("candidate_transcripts", len(bambu)),
        (
            "reference_bambu_exon_mismatches",
            mismatches
        ),
    ]

    for status, number in sorted(counts.items()):
        metrics.append((status, number))

    write_tsv(
        outdir / "isoform_structure_qc.tsv",
        ["metric", "value"],
        [
            {"metric": key, "value": value}
            for key, value in metrics
        ]
    )

    print("===== STRUCTURE ANNOTATION QC =====")

    for key, value in metrics:
        print(f"{key}\t{value}")

    print(
        "CDS coordinates are reference annotations, "
        "not translated protein predictions."
    )


if __name__ == "__main__":
    main()
