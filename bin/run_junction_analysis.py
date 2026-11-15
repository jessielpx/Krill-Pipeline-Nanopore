#!/usr/bin/env python3

"""
Krill Renew 11 - Long-read splice junction quantification.

Inputs:
  BAM manifest TSV: sample, bam
  Ensembl GTF

Outputs:
  junction_counts.tsv
  junction_annotation.tsv
  junction_sample_qc.tsv

Coordinates are 1-based, inclusive intron coordinates.

Counts are junction-supporting primary alignments,
NOT unique molecules or independent AS events.

No strand inferred from cDNA read orientation.
"""

import argparse
import csv
import re
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path


CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
ATTR_RE = re.compile(
    r'([A-Za-z_][A-Za-z0-9_]*) "([^"]*)"'
)


def read_tsv(path):
    with open(path, newline="") as handle:
        reader = csv.DictReader(
            handle, delimiter="\t"
        )
        return list(reader), reader.fieldnames or []


def write_tsv(path, columns, rows):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            delimiter="\t",
            extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)


def parse_gtf(gtf_path):
    """
    Build reference introns from adjacent exons
    of individual annotated transcripts.
    """

    exons = defaultdict(list)
    tx_info = {}

    opener = open

    if str(gtf_path).endswith(".gz"):
        import gzip
        opener = gzip.open

    with opener(gtf_path, "rt") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 9 or fields[2] != "exon":
                continue

            attr = dict(ATTR_RE.findall(fields[8]))

            transcript = attr.get("transcript_id")
            gene = attr.get("gene_id")

            if not transcript or not gene:
                continue

            chrom = fields[0]
            strand = fields[6]

            start = int(fields[3])
            end = int(fields[4])

            if transcript in tx_info:
                if tx_info[transcript][:3] != (
                    gene, chrom, strand
                ):
                    raise ValueError(
                        "Inconsistent transcript coordinates: "
                        + transcript
                    )
            else:
                tx_info[transcript] = (
                    gene,
                    chrom,
                    strand,
                    attr.get("gene_name", "")
                )

            exons[transcript].append(
                (start, end)
            )

    junction_info = defaultdict(
        lambda: {
            "genes": set(),
            "names": set(),
            "strands": set(),
            "transcripts": set(),
        }
    )

    for transcript, blocks in exons.items():

        if len(blocks) < 2:
            continue

        gene, chrom, strand, name = tx_info[transcript]

        ordered = sorted(set(blocks))

        for left, right in zip(
            ordered[:-1], ordered[1:]
        ):
            intron_start = left[1] + 1
            intron_end = right[0] - 1

            if intron_end < intron_start:
                continue

            key = (
                chrom,
                intron_start,
                intron_end
            )

            info = junction_info[key]

            info["genes"].add(gene)
            if name:
                info["names"].add(name)
            info["strands"].add(strand)
            info["transcripts"].add(transcript)

    return junction_info


def extract_cigar_junctions(
    chrom, start, cigar, min_intron
):
    """
    SAM POS is 1-based.
    The CIGAR N operation advances the reference.
    """

    position = start
    observed = []

    operations = CIGAR_RE.findall(cigar)

    if not operations:
        return observed

    reconstructed = "".join(
        n + op for n, op in operations
    )

    if reconstructed != cigar:
        raise ValueError(
            "Invalid CIGAR: " + cigar
        )

    for length_text, op in operations:

        length = int(length_text)

        if op == "N":

            if length >= min_intron:
                observed.append(
                    (
                        chrom,
                        position,
                        position + length - 1
                    )
                )

            position += length

        elif op in ("M", "D", "=", "X"):
            position += length

        # I, S, H, P do not advance
        # the reference coordinate.

    return set(observed)


def count_bam(
    bam_path,
    samtools,
    min_mapq,
    min_intron
):
    """
    Exclude:
      0x4   unmapped
      0x100 secondary
      0x800 supplementary

    Keep primary alignments only.
    """

    command = [
        samtools,
        "view",
        "-F", "2308",
        str(bam_path)
    ]

    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=1
    )

    counts = Counter()

    qc = Counter()

    try:
        for line in process.stdout:

            qc["primary_alignments"] += 1

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 11:
                raise ValueError(
                    "Malformed SAM alignment"
                )

            mapq = int(fields[4])

            if mapq < min_mapq:
                qc["low_mapq_alignments"] += 1
                continue

            qc["eligible_alignments"] += 1

            chrom = fields[2]
            position = int(fields[3])
            cigar = fields[5]

            junctions = extract_cigar_junctions(
                chrom,
                position,
                cigar,
                min_intron
            )

            if junctions:
                qc["spliced_alignments"] += 1

            for junction in junctions:
                counts[junction] += 1
                qc["junction_observations"] += 1

    except Exception:
        process.kill()
        process.communicate()
        raise

    stderr = process.stderr.read()
    returncode = process.wait()

    if returncode != 0:
        raise RuntimeError(
            f"samtools view failed for {bam_path}\n"
            + stderr
        )

    qc["distinct_junctions"] = len(counts)

    return counts, qc


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--manifest", required=True
    )
    parser.add_argument(
        "--reference-gtf", required=True
    )
    parser.add_argument(
        "--outdir", required=True
    )
    parser.add_argument(
        "--samtools", default="samtools"
    )
    parser.add_argument(
        "--min-mapq", type=int, default=10
    )
    parser.add_argument(
        "--min-intron", type=int, default=20
    )

    args = parser.parse_args()

    if args.min_mapq < 0 or args.min_intron < 1:
        raise ValueError(
            "Invalid MAPQ or intron threshold"
        )

    outdir = Path(args.outdir)
    outdir.mkdir(
        parents=True, exist_ok=True
    )

    manifest, cols = read_tsv(args.manifest)

    if not {"sample", "bam"}.issubset(cols):
        raise ValueError(
            "Manifest requires sample and bam columns"
        )

    if not manifest:
        raise ValueError("Empty BAM manifest")

    samples = [
        row["sample"].strip()
        for row in manifest
    ]

    if (
        any(not s for s in samples)
        or len(samples) != len(set(samples))
    ):
        raise ValueError(
            "Blank or duplicate sample IDs"
        )

    all_counts = {}
    qc_rows = []
    observed = set()

    for entry in manifest:

        sample = entry["sample"].strip()
        bam = Path(entry["bam"].strip())

        if not bam.is_file():
            raise FileNotFoundError(bam)

        print(
            "Counting junctions:",
            sample,
            flush=True
        )

        counts, qc = count_bam(
            bam,
            args.samtools,
            args.min_mapq,
            args.min_intron
        )

        all_counts[sample] = counts
        observed.update(counts)

        qc_rows.append({
            "sample": sample,
            **qc
        })

    print(
        "Loading Ensembl junction annotation",
        flush=True
    )

    reference = parse_gtf(
        args.reference_gtf
    )

    annotation_rows = []
    matrix_rows = []

    for chrom, start, end in sorted(observed):

        junction_id = (
            f"{chrom}:{start}-{end}"
        )

        ref = reference.get(
            (chrom, start, end)
        )

        if ref:

            genes = sorted(ref["genes"])
            gene_names = sorted(ref["names"])
            strands = sorted(ref["strands"])

            known = "TRUE"

            uniquely_assigned = (
                len(genes) == 1
            )

            gene_id = (
                genes[0]
                if uniquely_assigned else ""
            )

        else:
            genes = []
            gene_names = []
            strands = []
            known = "FALSE"
            gene_id = ""

        annotation_rows.append({
            "junction_id": junction_id,
            "chrom": chrom,
            "intron_start": start,
            "intron_end": end,
            "reference_known": known,
            "gene_id": gene_id,
            "all_gene_ids": ";".join(genes),
            "gene_names": ";".join(gene_names),
            "reference_strands": ";".join(strands),
            "n_reference_genes": len(genes)
        })

        row = {
            "junction_id": junction_id,
            "gene_id": gene_id,
            "reference_known": known
        }

        for sample in samples:
            row[sample] = all_counts[sample].get(
                (chrom, start, end),
                0
            )

        matrix_rows.append(row)

    write_tsv(
        outdir / "junction_counts.tsv",
        [
            "junction_id",
            "gene_id",
            "reference_known"
        ] + samples,
        matrix_rows
    )

    write_tsv(
        outdir / "junction_annotation.tsv",
        [
            "junction_id",
            "chrom",
            "intron_start",
            "intron_end",
            "reference_known",
            "gene_id",
            "all_gene_ids",
            "gene_names",
            "reference_strands",
            "n_reference_genes"
        ],
        annotation_rows
    )

    qc_columns = [
        "sample",
        "primary_alignments",
        "low_mapq_alignments",
        "eligible_alignments",
        "spliced_alignments",
        "junction_observations",
        "distinct_junctions"
    ]

    write_tsv(
        outdir / "junction_sample_qc.tsv",
        qc_columns,
        qc_rows
    )

    summary = [
        ("samples", len(samples)),
        ("observed_junctions", len(observed)),
        (
            "reference_known_junctions",
            sum(
                x["reference_known"] == "TRUE"
                for x in annotation_rows
            )
        ),
        (
            "unique_gene_assigned_junctions",
            sum(
                bool(x["gene_id"])
                for x in annotation_rows
            )
        ),
        (
            "unannotated_junctions",
            sum(
                x["reference_known"] == "FALSE"
                for x in annotation_rows
            )
        ),
    ]

    write_tsv(
        outdir / "junction_qc_summary.tsv",
        ["metric", "value"],
        [
            {"metric": k, "value": v}
            for k, v in summary
        ]
    )

    for key, value in summary:
        print(f"{key}\t{value}")

    print(
        "Junction quantification completed",
        flush=True
    )


if __name__ == "__main__":
    main()
