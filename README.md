# Krill-Pipeline-Nanopore

**Pipeline that can be digested well by "Rorqual"**

A modular **Nextflow DSL2** pipeline for Oxford Nanopore long-read RNA sequencing analysis.

Krill-Pipeline-Nanopore provides an automated workflow for sequencing quality control, read preprocessing, fusion detection, transcriptome analysis, and downstream differential and isoform analyses.

---

# Workflow

```text
                 sequencing_summary.txt
                           │
                           ▼
                        pycoQC
                           │
                           ▼
                        FASTCAT
                    ┌──────────────┐
                    │              │
                    ▼              ▼
                 JAFFAL      Genome Alignment
                                 │
                                 ▼
                        Transcript Assembly
                                 │
                                 ▼
                      Transcript Quantification
                                 │
                                 ▼
                     Differential Expression
                                 │
                                 ▼
                    IsoformSwitchAnalyzeR
```

---

# Features

- Oxford Nanopore RNA sequencing workflow
- Nextflow DSL2 modular design
- pycoQC sequencing QC
- FASTCAT read merging
- JAFFAL fusion detection
- Genome alignment
- Transcript assembly
- Transcript quantification
- Differential expression analysis
- Isoform switching analysis
- Easily extensible with additional modules

---

# Requirements

The pipeline requires:

- Nextflow (>=25)
- Java 21
- Apptainer (Singularity)
- Slurm (recommended)

Example module loading (Alliance Rorqual):

```bash
module load StdEnv/2023
module load java/21.0.1
module load nextflow/25.10.2
module load apptainer/1.4.5
```

---

# Installation

Clone the repository.

```bash
git clone git@github.com:jessielpx/Krill-Pipeline-Nanopore.git

cd Krill-Pipeline-Nanopore
```

---

# Quick Start

The pipeline can be started in **five steps**.

## Step 1. Prepare your sequencing data

Your FASTQ directory should look like

```text
batch1-8/

├── barcode01/
│      *.fastq.gz
│
├── barcode02/
│      *.fastq.gz
│
├── barcode03/
│      *.fastq.gz
│
└── ...
```

Each barcode folder contains FASTQ files from one barcode.

You should also have the sequencing summary file produced during basecalling.

Example

```text
sequencing_summary_run1.txt
```

---

## Step 2. Prepare the sample sheet

Edit

```text
assets/samples.csv
```

Example

| barcode | sample | alias | condition |
|---------|---------|---------|-----------|
| barcode01 | sample20 | sample20 | LR |
| barcode02 | sample27 | sample27 | LR |
| barcode03 | sample74 | sample74 | HR |

### Notes

- `barcode` must exactly match the FASTQ folder name.
- `sample` should be unique.
- `condition` is used for downstream differential analysis.

---

## Step 3. Download the JAFFAL reference

Reference files are **not included** in this repository because they are several gigabytes in size.

Run the setup script once:

```bash
bash scripts/download_jaffal_reference.sh \
    /path/to/reference_directory
```

The script will

- download the official JAFFAL hg38 + GENCODE49 reference
- extract the archive
- check all required files
- build Bowtie2 indexes if necessary
- print the reference directory

This only needs to be done once.

---

## Step 4. Configure the reference

Open

```text
nextflow.config
```

Update

```groovy
params {

    jaffal_ref_dir = "/path/to/JAFFAL_reference"

}
```

---

## Step 5. Run the pipeline

Example

```bash
nextflow run . \
    -profile rorqual \
    --summary /path/to/sequencing_summary_run1.txt \
    --fastq_dir /path/to/batch1-8
```

This executes

```text
pycoQC

↓

FASTCAT

↓

JAFFAL
```

---

# Running the complete transcriptome workflow

The complete workflow additionally requires

- reference genome FASTA
- annotation GTF

Example

```bash
nextflow run . \
    -profile rorqual \
    --summary /path/to/sequencing_summary_run1.txt \
    --fastq_dir /path/to/batch1-8 \
    --ref_genome genome.fa \
    --ref_annotation genes.gtf
```

The workflow becomes

```text
pycoQC

↓

FASTCAT

↓

Genome alignment

↓

Transcript assembly

↓

Transcript quantification

↓

Differential expression

↓

IsoformSwitchAnalyzeR
```

---

# Input

Required

| Input | Description |
|--------|-------------|
| sequencing_summary.txt | Nanopore sequencing summary |
| FASTQ directory | Barcode folders containing FASTQ files |
| samples.csv | Sample metadata |

Additional for transcriptome workflow

| Input | Description |
|--------|-------------|
| Genome FASTA | Reference genome |
| Annotation GTF | Gene annotation |

---

# Output

Typical output structure

```text
results/

├── pycoqc/
│
├── fastcat/
│
├── jaffal/
│
├── alignment/
│
├── transcriptome/
│
├── differential_expression/
│
└── isoform_switch/
```

---

# JAFFAL output

The primary fusion result is

```text
sampleX.jaffa_results.csv
```

Important columns

| Column | Description |
|---------|-------------|
| Fusion | Candidate fusion genes |
| Classification | HighConfidence / LowConfidence / PotentialTransSplicing |
| Spanning Reads | Number of supporting reads |
| Rearrangement | Genomic rearrangement evidence |
| Known COSMIC | COSMIC annotation |
| GTEx Samples | Number of normal GTEx samples |

---

# Troubleshooting

## Missing sequencing summary

Provide

```text
--summary sequencing_summary_run1.txt
```

---

## Missing FASTQ directory

Provide

```text
--fastq_dir path/to/batch1-8
```

---

## JAFFAL cannot find the reference

Check

```groovy
params.jaffal_ref_dir
```

in

```text
nextflow.config
```

---

## Exit status 137

The process exceeded available memory.

Increase the memory allocated to the corresponding Nextflow process.

---

# Repository structure

```text
Krill-Pipeline-Nanopore/

├── assets/
│
├── docs/
│
├── modules/
│
├── scripts/
│      download_jaffal_reference.sh
│
├── tests/
│
├── main.nf
│
├── nextflow.config
│
└── README.md
```

---

# Citation

If you use JAFFAL, please cite

Davidson NM et al.

**JAFFAL: Detecting fusion genes with long-read transcriptome sequencing.**

Genome Biology (2022)

Please also cite the original software used in this workflow, including

- Nextflow
- pycoQC
- FASTCAT
- JAFFAL
- minimap2
- StringTie
- Salmon
- GFFCompare
- IsoformSwitchAnalyzeR

---

# Roadmap

Planned future features

- [ ] Automatic reference download
- [ ] MultiQC report
- [ ] Differential transcript usage
- [ ] Gene fusion visualization
- [ ] nf-core style parameter validation
- [ ] Docker support
- [ ] CI testing with GitHub Actions
