#!/usr/bin/env python3
"""Krill Renew 08: batch Ensembl 104 structural annotation (no ORF prediction)."""
import argparse, csv, gzip, re
from collections import defaultdict
from pathlib import Path

PATTERN = re.compile(r'([A-Za-z_][A-Za-z0-9_]*) "([^"]*)"')

def reader(path):
    op = gzip.open if str(path).endswith('.gz') else open
    return op(path, 'rt', encoding='utf-8', errors='replace')

def gtf_attrs(s):
    return dict(PATTERN.findall(s))

def load_gtf(path, wanted):
    result = defaultdict(lambda: {'exon': [], 'CDS': [], 'five_prime_utr': [], 'three_prime_utr': [], 'stop_codon': []})
    with reader(path) as handle:
        for line in handle:
            if not line or line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9:
                continue
            kind = f[2]
            if kind not in ('transcript', 'exon', 'CDS', 'five_prime_utr', 'three_prime_utr', 'stop_codon'):
                continue
            attrs = gtf_attrs(f[8])
            tid = attrs.get('transcript_id')
            if tid not in wanted:
                continue
            d = result[tid]
            for key in ('gene_id', 'gene_name', 'transcript_biotype', 'protein_id'):
                if attrs.get(key):
                    d[key] = attrs[key]
            d['strand'] = f[6]
            d['chrom'] = f[0]
            if kind != 'transcript':
                d[kind].append((int(f[3]), int(f[4]), f[7] if kind == 'CDS' else '.'))
    return result

def signature(features, include_phase=False):
    if include_phase:
        return tuple(sorted(features))
    return tuple(sorted((a, b) for a, b, _ in features))

def feature_len(features):
    return sum(b - a + 1 for a, b, _ in features)

def save_tsv(path, cols, rows):
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, delimiter='\t', fieldnames=cols, extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--candidates', required=True, help='isoform_candidate_summary_v3.tsv')
    p.add_argument('--usage', required=True, help='isoform_usage_supported.tsv')
    p.add_argument('--reference-gtf', required=True)
    p.add_argument('--bambu-gtf', required=True)
    p.add_argument('--outdir', required=True)
    a = p.parse_args()
    outdir = Path(a.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    with open(a.candidates, newline='') as f:
        candidates = list(csv.DictReader(f, delimiter='\t'))
    if not candidates or not {'gene_id', 'top_transcript'}.issubset(candidates[0]):
        raise ValueError('Candidates must include gene_id and top_transcript')
    genes = {row['gene_id'] for row in candidates}
    usage = defaultdict(set)
    with open(a.usage, newline='') as f:
        r = csv.DictReader(f, delimiter='\t')
        if not {'gene_id', 'feature_id'}.issubset(r.fieldnames or []):
            raise ValueError('Usage must include gene_id and feature_id')
        for row in r:
            if row['gene_id'] in genes:
                usage[row['gene_id']].add(row['feature_id'])
    wanted = set().union(*usage.values()) if usage else set()
    wanted.update(row['top_transcript'] for row in candidates)
    print('Candidate genes:', len(genes), 'Transcripts:', len(wanted), flush=True)
    print('Reading Ensembl reference ...', flush=True)
    ref = load_gtf(a.reference_gtf, wanted)
    print('Reading Bambu annotation ...', flush=True)
    bambu = load_gtf(a.bambu_gtf, wanted)

    columns = ['gene_id','transcript_id','gene_name','biotype','annotation_status','bambu_ref_exon_match',
               'chrom','strand','bambu_exons','reference_exons','reference_spliced_length',
               'cds_segments','cds_length_excluding_stop','stop_codon_length','utr5_length','utr3_length',
               'protein_id','comparison_to_top','same_annotated_cds_as_top']
    rows = []
    mismatch = 0
    for c in candidates:
        gid = c['gene_id']
        top = c['top_transcript']
        top_ref = ref.get(top)
        for tx in sorted(usage[gid] | {top}):
            r = ref.get(tx)
            b = bambu.get(tx)
            b_exons = signature(b['exon']) if b else ()
            r_exons = signature(r['exon']) if r else ()
            exon_match = bool(b_exons and r_exons and b_exons == r_exons)
            if b_exons and r_exons and not exon_match:
                mismatch += 1
            ref_ok = bool(r and r['exon'] and exon_match)
            cds = r['CDS'] if ref_ok else []
            top_ok = bool(top_ref and top_ref['exon'] and bambu.get(top) and
                          signature(top_ref['exon']) == signature(bambu[top]['exon']))
            same_cds = 'NA'
            comparison = 'not_comparable'
            if tx == top:
                comparison = 'top_transcript'
            elif ref_ok and top_ok:
                comparison = 'same_exon_chain' if signature(r['exon']) == signature(top_ref['exon']) else 'different_exon_chain'
            if ref_ok and top_ok and cds and top_ref['CDS']:
                same_cds = 'TRUE' if signature(cds, True) == signature(top_ref['CDS'], True) else 'FALSE'
            if not r:
                status = 'not_in_ensembl104'
            elif not b or not b['exon']:
                status = 'no_bambu_exons'
            elif not exon_match:
                status = 'exon_chain_mismatch'
            elif not cds:
                status = 'matched_no_CDS'
            else:
                status = 'matched_with_CDS'
            rows.append(dict(
                gene_id=gid, transcript_id=tx, gene_name=(r or b or {}).get('gene_name', ''),
                biotype=(r or {}).get('transcript_biotype', ''), annotation_status=status,
                bambu_ref_exon_match='TRUE' if exon_match else 'FALSE',
                chrom=(b or r or {}).get('chrom', ''), strand=(b or r or {}).get('strand', ''),
                bambu_exons=len(b['exon']) if b else 0, reference_exons=len(r['exon']) if r else 0,
                reference_spliced_length=feature_len(r['exon']) if ref_ok else 'NA',
                cds_segments=len(cds), cds_length_excluding_stop=feature_len(cds) if cds else 'NA',
                stop_codon_length=feature_len(r['stop_codon']) if ref_ok and r['stop_codon'] else 'NA',
                utr5_length=feature_len(r['five_prime_utr']) if ref_ok and r['five_prime_utr'] else 'NA',
                utr3_length=feature_len(r['three_prime_utr']) if ref_ok and r['three_prime_utr'] else 'NA',
                protein_id=(r or {}).get('protein_id', '') if ref_ok else '',
                comparison_to_top=comparison, same_annotated_cds_as_top=same_cds))
    save_tsv(outdir/'isoform_structure_annotations.tsv', columns, rows)
    counter = defaultdict(int)
    for row in rows:
        counter[row['annotation_status']] += 1
    summary = [{'metric':'candidate_genes', 'value':len(genes)},
               {'metric':'candidate_transcripts', 'value':len(rows)},
               {'metric':'reference_bambu_exon_mismatches', 'value':mismatch}]
    summary += [{'metric':k, 'value':v} for k,v in sorted(counter.items())]
    save_tsv(outdir/'isoform_structure_qc.tsv', ['metric','value'], summary)
    print('Saved:', outdir/'isoform_structure_annotations.tsv', flush=True)
    print('Statuses:', dict(counter), flush=True)
    print('Note: CDS coordinates do not prove identical translated protein sequences.', flush=True)

if __name__ == '__main__':
    main()
