#!/usr/bin/env python3

"""
Krill Renew - Gene Expression Export

Inputs:
  Bambu counts_gene.txt
  Bambu counts_transcript.txt
  Bambu extended_annotations.gtf
  Ensembl reference GTF

Outputs:
  Gene-level counts and TPM
  Ensembl ID and gene-symbol versions
  Gene mapping and export QC

TPM method:
  transcript estimated counts / spliced transcript length
  -> sample-wise TPM normalization
  -> sum transcript TPM by gene

Standard Python library only.
"""

import csv
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path


def fail(message):
    raise SystemExit("ERROR: " + message)


def attributes(text):
    return dict(
        re.findall(r'([\w.-]+)\s+"([^"]*)"', text)
    )


def read_reference_symbols(gtf):
    symbols = {}

    with gtf.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")

            if len(parts) < 9:
                continue

            if parts[2] != "gene":
                continue

            attr = attributes(parts[8])
            gene = attr.get("gene_id", "")
            symbol = attr.get("gene_name", "")

            if gene and symbol:
                symbols[gene] = symbol
                symbols.setdefault(
                    gene.split(".")[0], symbol
                )

    return symbols


def transcript_annotations(gtf):
    tx_gene = {}
    exons = defaultdict(list)

    with gtf.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")

            if len(parts) < 9:
                continue

            feature = parts[2]

            if feature not in {"transcript", "exon"}:
                continue

            attr = attributes(parts[8])
            tx = attr.get("transcript_id")
            gene = attr.get("gene_id")

            if not tx or not gene:
                continue

            previous = tx_gene.get(tx)

            if previous is not None and previous != gene:
                fail("Conflicting GTF gene for " + tx)

            tx_gene[tx] = gene

            if feature == "exon":
                start = int(parts[3])
                end = int(parts[4])

                if end < start:
                    fail("Invalid exon coordinates: " + tx)

                exons[tx].append((start, end))

    lengths = {}

    for tx, intervals in exons.items():
        intervals.sort()

        start, end = intervals[0]
        total = 0

        for left, right in intervals[1:]:
            if left <= end + 1:
                end = max(end, right)
            else:
                total += end - start + 1
                start, end = left, right

        total += end - start + 1
        lengths[tx] = total

    return tx_gene, lengths


def read_gene_counts(path):
    genes = {}
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)

        if header[0] != "GENEID":
            fail("Invalid gene counts header")

        samples = header[1:]

        if not samples or len(samples) != len(set(samples)):
            fail("Invalid sample columns")

        for row in reader:
            if not row:
                continue

            gene = row[0]

            if gene in genes:
                fail("Duplicate gene ID: " + gene)

            if len(row) != len(header):
                fail("Invalid gene counts row: " + gene)

            genes[gene] = [float(v) for v in row[1:]]

    return samples, genes


def write_matrix(path, samples, matrix, key_name):
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([key_name] + samples)

        for key in sorted(matrix):
            writer.writerow(
                [key] + [
                    f"{v:.6f}" for v in matrix[key]
                ]
            )


def main():
    if len(sys.argv) != 4:
        fail(
            "Usage: export_gene_expression.py "
            "<bambu_results> <reference.gtf> <outdir>"
        )

    bambu = Path(sys.argv[1])
    reference = Path(sys.argv[2])
    outdir = Path(sys.argv[3])
    outdir.mkdir(parents=True, exist_ok=True)

    print("===== READING GENE COUNTS =====", flush=True)

    samples, gene_counts = read_gene_counts(
        bambu / "counts_gene.txt"
    )

    print("Samples:", samples, flush=True)
    print("Genes:", len(gene_counts), flush=True)

    print("===== READING REFERENCE SYMBOLS =====", flush=True)

    symbols = read_reference_symbols(reference)

    print("Reference mappings:", len(symbols), flush=True)

    print("===== READING BAMBU TRANSCRIPT STRUCTURE =====",
          flush=True)

    tx_gene, tx_lengths = transcript_annotations(
        bambu / "extended_annotations.gtf"
    )

    print("Transcript lengths:", len(tx_lengths), flush=True)

    print("===== COMPUTING TRANSCRIPT-DERIVED TPM =====",
          flush=True)

    n = len(samples)
    gene_rates = {}
    total_rates = [0.0] * n
    processed = 0
    missing_tx = []
    gene_mismatches = []

    path = bambu / "counts_transcript.txt"

    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        expected = ["TXNAME", "GENEID"] + samples

        if reader.fieldnames != expected:
            fail("Gene and transcript sample columns differ")

        seen_tx = set()

        for row in reader:
            tx = row["TXNAME"]
            gene = row["GENEID"]

            if tx in seen_tx:
                fail("Duplicate transcript count ID: " + tx)
            seen_tx.add(tx)

            if gene not in gene_counts:
                fail("Transcript gene missing in gene counts: " + gene)

            if tx not in tx_lengths:
                missing_tx.append(tx)
                continue

            if tx_gene[tx] != gene:
                gene_mismatches.append(tx)
                continue

            length_kb = tx_lengths[tx] / 1000.0

            if length_kb <= 0:
                fail("Invalid transcript length: " + tx)

            values = [
                float(row[s]) / length_kb
                for s in samples
            ]

            if gene not in gene_rates:
                gene_rates[gene] = [0.0] * n

            for i, value in enumerate(values):
                gene_rates[gene][i] += value
                total_rates[i] += value

            processed += 1

    if missing_tx or gene_mismatches:
        fail(
            f"Missing transcript lengths: {len(missing_tx)}; "
            f"gene mismatches: {len(gene_mismatches)}. "
            f"First missing: {missing_tx[:5]}; "
            f"first mismatches: {gene_mismatches[:5]}"
        )

    if len(seen_tx) != len(tx_gene):
        fail(
            "GTF transcript count differs from counts_transcript: "
            f"{len(tx_gene)} vs {len(seen_tx)}"
        )

    if any(v <= 0 for v in total_rates):
        fail("A sample has zero total transcript expression")

    gene_tpm = {}

    for gene in gene_counts:
        rates = gene_rates.get(gene, [0.0] * n)

        gene_tpm[gene] = [
            rates[i] / total_rates[i] * 1e6
            for i in range(n)
        ]

    print("Transcripts processed:", processed, flush=True)

    print("===== MAPPING GENE SYMBOLS =====", flush=True)

    symbol_counts = {}
    symbol_tpm = {}
    symbol_sources = defaultdict(list)
    unmapped = []

    for gene in gene_counts:
        symbol = symbols.get(gene)

        if symbol is None:
            symbol = symbols.get(gene.split(".")[0])

        if not symbol:
            unmapped.append(gene)
            continue

        symbol_sources[symbol].append(gene)

        if symbol not in symbol_counts:
            symbol_counts[symbol] = [0.0] * n
            symbol_tpm[symbol] = [0.0] * n

        for i in range(n):
            symbol_counts[symbol][i] += gene_counts[gene][i]
            symbol_tpm[symbol][i] += gene_tpm[gene][i]

    print("===== WRITING MATRICES =====", flush=True)

    write_matrix(
        outdir / "gene_counts_ENSG.tsv",
        samples, gene_counts, "GENEID"
    )

    write_matrix(
        outdir / "gene_counts_symbols.tsv",
        samples, symbol_counts, "GeneSymbol"
    )

    write_matrix(
        outdir / "gene_TPM_ENSG.tsv",
        samples, gene_tpm, "GENEID"
    )

    write_matrix(
        outdir / "gene_TPM_symbols.tsv",
        samples, symbol_tpm, "GeneSymbol"
    )

    print("===== WRITING MAPPING QC =====", flush=True)

    with (outdir / "gene_id_mapping.tsv").open(
        "w", newline=""
    ) as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "GENEID", "GeneSymbol", "mapping_status",
            "symbol_gene_id_count"
        ])

        for gene in sorted(gene_counts):
            symbol = (
                symbols.get(gene)
                or symbols.get(gene.split(".")[0])
                or ""
            )

            writer.writerow([
                gene,
                symbol,
                "mapped" if symbol else "unmapped",
                len(symbol_sources[symbol]) if symbol else 0
            ])

    duplicated_symbols = sum(
        len(ids) > 1
        for ids in symbol_sources.values()
    )

    with (outdir / "export_qc_summary.tsv").open(
        "w", newline=""
    ) as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])

        metrics = {
            "samples": n,
            "bambu_genes": len(gene_counts),
            "transcripts_processed": processed,
            "genes_with_symbols": (
                len(gene_counts) - len(unmapped)
            ),
            "genes_without_symbols": len(unmapped),
            "unique_gene_symbols": len(symbol_counts),
            "symbols_mapping_to_multiple_gene_ids": (
                duplicated_symbols
            ),
            "gene_tpm_method": (
                "transcript_estimated_counts_divided_"
                "by_spliced_length_then_gene_sum"
            ),
        }

        for key, value in metrics.items():
            writer.writerow([key, value])

        for i, sample in enumerate(samples):
            writer.writerow([
                f"{sample}_all_gene_TPM_sum",
                f"{sum(v[i] for v in gene_tpm.values()):.4f}"
            ])
            writer.writerow([
                f"{sample}_mapped_symbol_TPM_sum",
                f"{sum(v[i] for v in symbol_tpm.values()):.4f}"
            ])

    print("===== EXPORT COMPLETE =====", flush=True)
    print("Mapped gene IDs:", len(gene_counts) - len(unmapped))
    print("Unique symbols:", len(symbol_counts))
    print("Unmapped gene IDs:", len(unmapped))
    print("Duplicate symbols:", duplicated_symbols)
    print("Output:", outdir)


if __name__ == "__main__":
    main()
