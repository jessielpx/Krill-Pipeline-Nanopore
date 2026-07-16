# Pipeline_Nanopore

A modular Nextflow DSL2 pipeline for Oxford Nanopore long-read RNA sequencing analysis.

## Features

* Run-level quality control (pycoQC)
* FASTQ concatenation (FASTCAT)
* Genome alignment (Minimap2 + Samtools)
* Transcript assembly (StringTie)
* Transcript quantification (Salmon)
* Gene and transcript count matrices
* Differential gene expression (edgeR)
* Differential transcript usage (DRIMSeq / DEXSeq)
* Transcript annotation (gffcompare)

## Workflow

```text
FASTQ
 │
 ▼
FASTCAT
 │
 ▼
Genome alignment
 │
 ▼
StringTie assembly
 │
 ▼
Merged transcriptome
 │
 ▼
Transcriptome alignment
 │
 ▼
Salmon quantification
 │
 ▼
Count matrices
 │
 ├── Differential expression
 ├── Differential transcript usage
 └── gffcompare annotation
```

## Requirements

* Nextflow (DSL2)
* Java
* Apptainer
* Slurm (Rorqual profile)

## Run

```bash
nextflow run . \
    -profile rorqual \
    --summary sequencing_summary.txt \
    --fastq_dir batch1-8 \
    --ref_genome Homo_sapiens.GRCh38.fa \
    --ref_annotation genes.filtered_to_fasta.gtf \
    --run_id Batch1 \
    --outdir results \
    --de_analysis true
```

## Main outputs

```text
results/
├── pycoqc/
├── alignment/
├── bam/
├── transcriptome_alignment/
├── quantification/
├── differential_expression/
└── gffcompare/
```

## Repository structure

```text
modules/
subworkflows/
bin/
assets/
tests/
main.nf
nextflow.config
```

## Current version

**v0.1.0**

Tested on the Alliance Canada Rorqual cluster using Slurm and Apptainer.
