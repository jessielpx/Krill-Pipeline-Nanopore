#!/usr/bin/env python3

"""
Krill Renew 09 - Functional consequence prioritization.

Integrates:
  - Isoform usage
  - All-transcript structural annotation
  - Reference GTF CDS/exon/UTR coordinates
  - Gene-level DTU results
  - Sample metadata

This is structural consequence prioritization,
NOT experimental functional validation or
de novo ORF/protein prediction.
"""

import argparse
import csv
import gzip
import math
import re
from collections import Counter, defaultdict
from pathlib import Path


ATTR_RE = re.compile(
    r'([A-Za-z_][A-Za-z0-9_]*) "([^"]*)"'
)


def read_table(path):
    with open(path, newline="") as handle:
        first = handle.readline()
        handle.seek(0)

        delimiter = "\t" if "\t" in first else ","

        reader = csv.DictReader(
            handle,
            delimiter=delimiter
        )

        rows = list(reader)
        columns = reader.fieldnames or []

    return rows, columns


def require(columns, expected, name):
    missing = set(expected) - set(columns)

    if missing:
        raise ValueError(
            f"{name}: missing columns {sorted(missing)}"
        )


def numeric(value):
    try:
        x = float(value)
        return x if math.isfinite(x) else None
    except (TypeError, ValueError):
        return None


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


def open_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def load_reference_features(path, wanted):
    """
    Load feature coordinates for the transcripts
    participating in candidate pairs only.
    """

    features = defaultdict(
        lambda: {
            "chrom": "",
            "strand": "",
            "exon": [],
            "CDS": [],
            "five_prime_utr": [],
            "three_prime_utr": [],
        }
    )

    accepted = {
        "exon",
        "CDS",
        "five_prime_utr",
        "three_prime_utr",
    }

    with open_text(path) as handle:
        for line in handle:

            if not line.strip() or line.startswith("#"):
                continue

            f = line.rstrip("\n").split("\t")

            if len(f) < 9 or f[2] not in accepted:
                continue

            attrs = dict(ATTR_RE.findall(f[8]))
            transcript = attrs.get("transcript_id")

            if transcript not in wanted:
                continue

            item = features[transcript]

            item["chrom"] = f[0]
            item["strand"] = f[6]

            phase = f[7] if f[2] == "CDS" else "."

            item[f[2]].append(
                (int(f[3]), int(f[4]), phase)
            )

    return dict(features)


def signature(record, kind, phase=False):
    if record is None:
        return None

    coordinates = record[kind]

    if not coordinates:
        return None

    if phase:
        segments = tuple(sorted(coordinates))
    else:
        segments = tuple(
            sorted((a, b) for a, b, _ in coordinates)
        )

    return (
        record["chrom"],
        record["strand"],
        segments
    )


def compare_feature(a, b, kind, phase=False):
    """
    Return same / different / not_comparable.
    """

    sa = signature(a, kind, phase)
    sb = signature(b, kind, phase)

    if sa is None or sb is None:
        return "not_comparable"

    return "same" if sa == sb else "different"


def mean(values):
    return sum(values) / len(values) if values else None


def range_of(values):
    return (
        max(values) - min(values)
        if values else None
    )


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--samplesheet", required=True)
    parser.add_argument("--usage", required=True)
    parser.add_argument("--structure", required=True)
    parser.add_argument("--dtu", required=True)
    parser.add_argument("--reference-gtf", required=True)
    parser.add_argument("--outdir", required=True)

    parser.add_argument(
        "--min-reads",
        type=float,
        default=3
    )

    parser.add_argument(
        "--min-supported-samples",
        type=int,
        default=2
    )

    parser.add_argument(
        "--group-fraction",
        type=float,
        default=0.25
    )

    args = parser.parse_args()

    if args.min_reads < 0:
        raise ValueError("--min-reads must be >= 0")

    if args.min_supported_samples < 1:
        raise ValueError(
            "--min-supported-samples must be >= 1"
        )

    if not 0 < args.group_fraction <= 1:
        raise ValueError(
            "--group-fraction must be in (0,1]"
        )

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------
    # Metadata
    # ------------------------------------------

    metadata, mcols = read_table(args.samplesheet)

    require(
        mcols,
        ["sample", "condition"],
        "samplesheet"
    )

    samples = [
        row["sample"].strip()
        for row in metadata
    ]

    if not samples:
        raise ValueError("Empty samplesheet")

    if any(not s for s in samples):
        raise ValueError("Blank sample ID")

    if len(samples) != len(set(samples)):
        raise ValueError("Duplicate sample IDs")

    groups = defaultdict(list)

    for row in metadata:
        condition = row["condition"].strip()

        if not condition:
            raise ValueError("Blank condition")

        groups[condition].append(
            row["sample"].strip()
        )

    if len(groups) != 2:
        raise ValueError(
            "09 requires exactly two comparison groups "
            "because upstream delta_IF is a two-group "
            "contrast. Observed: "
            + ", ".join(sorted(groups))
        )

    if min(map(len, groups.values())) < 2:
        raise ValueError(
            "At least two samples per group required"
        )

    group_names = sorted(groups)

    # ------------------------------------------
    # Upstream tables
    # ------------------------------------------

    usage, ucols = read_table(args.usage)
    structure, scols = read_table(args.structure)
    dtu, dcols = read_table(args.dtu)

    require(
        ucols,
        ["gene_id", "feature_id", "delta_IF"]
        + samples
        + [s + "_IF" for s in samples],
        "isoform usage"
    )

    require(
        scols,
        [
            "gene_id",
            "transcript_id",
            "annotation_status",
            "bambu_ref_exon_match"
        ],
        "structure"
    )

    require(
        dcols,
        ["gene_id", "pvalue", "adj_pvalue"],
        "DTU"
    )

    annotations = {}

    for row in structure:
        key = (
            row["gene_id"],
            row["transcript_id"]
        )

        if key in annotations:
            raise ValueError(
                f"Duplicate structure annotation: {key}"
            )

        annotations[key] = row

    dtu_by_gene = {}

    for row in dtu:
        gene = row["gene_id"]

        if gene in dtu_by_gene:
            raise ValueError(
                f"Duplicate DTU gene: {gene}"
            )

        dtu_by_gene[gene] = row

    # ------------------------------------------
    # Identify opposite-direction isoform pairs
    # ------------------------------------------

    by_gene = defaultdict(list)
    seen_usage = set()

    for row in usage:

        key = (
            row["gene_id"],
            row["feature_id"]
        )

        if key in seen_usage:
            raise ValueError(
                f"Duplicate usage transcript: {key}"
            )

        seen_usage.add(key)

        delta = numeric(row.get("delta_IF"))

        if delta is not None and delta != 0:
            by_gene[row["gene_id"]].append(
                (abs(delta), delta, row)
            )

    selected = []

    for gene, transcripts in by_gene.items():

        transcripts.sort(
            key=lambda item: (
                -item[0],
                item[2]["feature_id"]
            )
        )

        _, major_delta, major = transcripts[0]

        opposite = [
            item
            for item in transcripts[1:]
            if item[1] * major_delta < 0
        ]

        if not opposite:
            continue

        opposite.sort(
            key=lambda item: (
                -item[0],
                item[2]["feature_id"]
            )
        )

        _, partner_delta, partner = opposite[0]

        selected.append(
            (
                gene,
                major_delta,
                major,
                partner_delta,
                partner
            )
        )

    wanted = {
        row["feature_id"]
        for _, _, major, _, partner in selected
        for row in (major, partner)
    }

    reference = load_reference_features(
        args.reference_gtf,
        wanted
    )

    results = []
    missing_annotations = 0

    # ------------------------------------------
    # Pair annotation
    # ------------------------------------------

    for (
        gene,
        major_delta,
        major,
        partner_delta,
        partner
    ) in selected:

        major_id = major["feature_id"]
        partner_id = partner["feature_id"]

        a = annotations.get((gene, major_id), {})
        b = annotations.get((gene, partner_id), {})

        missing_annotations += int(not a) + int(not b)

        a_status = a.get(
            "annotation_status",
            "missing"
        )
        b_status = b.get(
            "annotation_status",
            "missing"
        )

        a_ref = (
            reference.get(major_id)
            if a.get("bambu_ref_exon_match") == "TRUE"
            else None
        )

        b_ref = (
            reference.get(partner_id)
            if b.get("bambu_ref_exon_match") == "TRUE"
            else None
        )

        # CDS coordinates may only be compared when
        # both transcript exon chains match reference.
        if (
            a_status == "matched_with_CDS"
            and b_status == "matched_with_CDS"
            and a_ref is not None
            and b_ref is not None
        ):
            cds_cmp = compare_feature(
                a_ref,
                b_ref,
                "CDS",
                phase=True
            )

            cds_relation = {
                "same": "same_CDS",
                "different": "different_CDS"
            }.get(cds_cmp, "not_comparable")

        else:
            cds_relation = "not_comparable"

        exon_relation = compare_feature(
            a_ref,
            b_ref,
            "exon"
        )

        utr5_relation = compare_feature(
            a_ref,
            b_ref,
            "five_prime_utr"
        )

        utr3_relation = compare_feature(
            a_ref,
            b_ref,
            "three_prime_utr"
        )

        row = {
            "gene_id": gene,
            "gene_name": (
                a.get("gene_name")
                or b.get("gene_name", "")
            ),
            "major_transcript": major_id,
            "partner_transcript": partner_id,
            "major_delta_IF": major_delta,
            "partner_delta_IF": partner_delta,
            "pair_effect": min(
                abs(major_delta),
                abs(partner_delta)
            ),
            "CDS_relation": cds_relation,
            "exon_chain_relation": exon_relation,
            "utr5_relation": utr5_relation,
            "utr3_relation": utr3_relation,
            "major_annotation_status": a_status,
            "partner_annotation_status": b_status,
            "major_biotype": a.get("biotype", ""),
            "partner_biotype": b.get("biotype", ""),
            "major_utr5": a.get("utr5_length", ""),
            "partner_utr5": b.get("utr5_length", ""),
            "major_utr3": a.get("utr3_length", ""),
            "partner_utr3": b.get("utr3_length", ""),
        }

        pair_support = []
        pair_consistency = []

        for prefix, transcript in (
            ("major", major),
            ("partner", partner)
        ):

            count_map = {}
            fraction_map = {}

            for sample in samples:

                count = numeric(
                    transcript.get(sample)
                )

                fraction = numeric(
                    transcript.get(sample + "_IF")
                )

                if count is None or count < 0:
                    raise ValueError(
                        f"Invalid count for "
                        f"{gene}/{transcript['feature_id']}/"
                        f"{sample}"
                    )

                if fraction is not None and not (
                    0 <= fraction <= 1 + 1e-6
                ):
                    raise ValueError(
                        f"Invalid isoform fraction for "
                        f"{gene}/{transcript['feature_id']}/"
                        f"{sample}"
                    )

                count_map[sample] = count
                fraction_map[sample] = fraction

                row[
                    f"{prefix}_{sample}_count"
                ] = count

                row[
                    f"{prefix}_{sample}_IF"
                ] = (
                    fraction
                    if fraction is not None
                    else "NA"
                )

            supported_groups = []

            ranges = {}

            for group, members in groups.items():

                minimum = max(
                    args.min_supported_samples,
                    math.ceil(
                        args.group_fraction * len(members)
                    )
                )

                supported_count = sum(
                    count_map[sample] >= args.min_reads
                    for sample in members
                )

                supported_groups.append(
                    supported_count >= minimum
                )

                defined_fractions = [
                    fraction_map[s]
                    for s in members
                    if fraction_map[s] is not None
                ]

                ranges[group] = range_of(
                    defined_fractions
                )

                row[
                    f"{prefix}_{group}_supported_samples"
                ] = supported_count

                row[
                    f"{prefix}_{group}_IF_range"
                ] = (
                    ranges[group]
                    if ranges[group] is not None
                    else "NA"
                )

            supported = all(supported_groups)

            delta = numeric(
                transcript.get("delta_IF")
            )

            consistent = (
                delta is not None
                and all(
                    value is not None
                    for value in ranges.values()
                )
                and abs(delta) > max(
                    ranges.values()
                )
            )

            row[
                f"{prefix}_total_unique"
            ] = sum(count_map.values())

            row[
                f"{prefix}_read_support"
            ] = supported

            row[
                f"{prefix}_effect_exceeds_group_ranges"
            ] = consistent

            pair_support.append(supported)
            pair_consistency.append(consistent)

        row["pair_support_pass"] = all(
            pair_support
        )

        row["both_isoforms_consistent"] = all(
            pair_consistency
        )

        # Descriptive QC only, not a statistical test.
        row["pair_qc_pass"] = (
            row["pair_support_pass"]
            and row["both_isoforms_consistent"]
        )

        dtu_result = dtu_by_gene.get(gene)

        row["DTU_tested"] = (
            dtu_result is not None
        )

        row["DTU_pvalue"] = (
            dtu_result.get("pvalue", "")
            if dtu_result else ""
        )

        row["DTU_adj_pvalue"] = (
            dtu_result.get("adj_pvalue", "")
            if dtu_result else ""
        )

        adj = numeric(row["DTU_adj_pvalue"])

        row["DTU_FDR_0_05"] = (
            adj is not None and adj < 0.05
        )

        results.append(row)

    results.sort(
        key=lambda row: (
            -row["pair_effect"],
            row["gene_id"]
        )
    )

    # ------------------------------------------
    # Output
    # ------------------------------------------

    static_columns = [
        "gene_id",
        "gene_name",
        "major_transcript",
        "partner_transcript",
        "major_delta_IF",
        "partner_delta_IF",
        "pair_effect",
        "CDS_relation",
        "exon_chain_relation",
        "utr5_relation",
        "utr3_relation",
        "major_annotation_status",
        "partner_annotation_status",
        "major_biotype",
        "partner_biotype",
        "major_utr5",
        "partner_utr5",
        "major_utr3",
        "partner_utr3",
    ]

    dynamic_columns = []

    for prefix in ("major", "partner"):

        for sample in samples:
            dynamic_columns.extend([
                f"{prefix}_{sample}_count",
                f"{prefix}_{sample}_IF"
            ])

        for group in group_names:
            dynamic_columns.extend([
                f"{prefix}_{group}_supported_samples",
                f"{prefix}_{group}_IF_range"
            ])

        dynamic_columns.extend([
            f"{prefix}_total_unique",
            f"{prefix}_read_support",
            f"{prefix}_effect_exceeds_group_ranges",
        ])

    columns = (
        static_columns
        + dynamic_columns
        + [
            "pair_support_pass",
            "both_isoforms_consistent",
            "pair_qc_pass",
            "DTU_tested",
            "DTU_pvalue",
            "DTU_adj_pvalue",
            "DTU_FDR_0_05",
        ]
    )

    passed = [
        row
        for row in results
        if row["pair_qc_pass"]
    ]

    write_table(
        outdir / "all_isoform_pairs.tsv",
        results,
        columns
    )

    write_table(
        outdir / "pair_qc_pass_pairs.tsv",
        passed,
        columns
    )

    statuses = Counter(
        row["CDS_relation"]
        for row in results
    )

    metrics = [
        ("comparison_groups", ";".join(group_names)),
        ("total_pairs", len(results)),
        ("genes_with_opposite_pairs", len(results)),
        (
            "pairs_with_read_support",
            sum(
                row["pair_support_pass"]
                for row in results
            )
        ),
        (
            "pairs_with_descriptive_consistency",
            sum(
                row["both_isoforms_consistent"]
                for row in results
            )
        ),
        ("pair_qc_pass_pairs", len(passed)),
        (
            "DTU_tested_pairs",
            sum(row["DTU_tested"] for row in results)
        ),
        (
            "DTU_FDR_0_05_pairs",
            sum(row["DTU_FDR_0_05"] for row in results)
        ),
        (
            "missing_structure_annotations_for_pair_members",
            missing_annotations
        ),
    ]

    for group in group_names:
        metrics.append(
            (f"{group}_sample_count", len(groups[group]))
        )

    for status, count in sorted(statuses.items()):
        metrics.append(
            (f"CDS_{status}", count)
        )

    write_table(
        outdir / "functional_qc_summary.tsv",
        [
            {"metric": metric, "value": value}
            for metric, value in metrics
        ],
        ["metric", "value"]
    )

    print("===== FUNCTIONAL CONSEQUENCES QC =====")

    for metric, value in metrics:
        print(f"{metric}\t{value}")

    print(
        "Interpretation: CDS relations compare "
        "reference genomic coordinates and phase. "
        "They do not establish protein sequence "
        "identity or functional effects."
    )


if __name__ == "__main__":
    main()
