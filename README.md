# Krill-Pipeline-Nanopore

Pipeline that can be digested well by "Rorqual"

**Krill-Pipeline-Nanopore** is a modular Nextflow DSL2 pipeline for Oxford Nanopore long-read RNA sequencing analysis.

Current workflow includes:

- pycoQC sequencing QC
- FASTCAT read merging
- JAFFAL fusion detection
- Genome alignment (minimap2)
- Transcript assembly (StringTie)
- Transcript quantification (Salmon)
- Differential expression analysis
- GFFCompare transcript annotation
- IsoformSwitchAnalyzeR

---

# Workflow

```
                sequencing_summary.txt
                          │
                          ▼
                       pycoQC
                          │
                          ▼
                       FASTCAT
                    ┌───────────┐
                    │           │
                    ▼           ▼
                 JAFFAL     Genome Alignment
                    │              │
                    │              ▼
                    │      Transcript Assembly
                    │              │
                    │              ▼
                    │      Transcript Quantification
                    │              │
                    │              ▼
                    │   Differential Expression
                    │              │
                    │              ▼
                    └──── IsoformSwitchAnalyzeR
```

---

# Requirements

This pipeline has been tested on the Alliance Canada Rorqual cluster.

Required software

- Nextflow >= 25
- Java 21
- Apptainer
- Slurm

Load required modules

```bash
module load StdEnv/2023
module load java/21.0.1
module load nextflow/25.10.2
module load apptainer/1.4.5
```

---

# Installation

Clone the repository

```bash
git clone git@github.com:jessielpx/Krill-Pipeline-Nanopore.git

cd Krill-Pipeline-Nanopore
```

---

# Before running the pipeline

You need to prepare:

✔ sequencing_summary.txt

✔ FASTQ barcode folders

✔ samples.csv

✔ JAFFAL reference

✔ (optional) reference genome

✔ (optional) annotation GTF

---

# Step 1. Prepare your FASTQ files

The FASTQ directory should look like this

```
batch1-8/

    barcode01/

        *.fastq.gz

    barcode02/

        *.fastq.gz

    barcode03/

        *.fastq.gz

    ...
```

Each barcode folder contains one sequencing barcode.

---

# Step 2. Prepare samples.csv

Edit

```
assets/samples.csv
```

Example

| barcode | sample | alias | condition |
|----------|--------|-------|-----------|
| barcode01 | sample20 | sample20 | LR |
| barcode02 | sample27 | sample27 | LR |
| barcode03 | sample74 | sample74 | HR |

Requirements

- barcode must exactly match the FASTQ folder name.
- sample should be unique.
- condition is used for downstream differential analysis.

---

# Step 3. Download the JAFFAL reference (one-time setup)

Download the official JAFFAL reference

```bash
mkdir -p ~/references

cd ~/references

curl -L \
https://api.figshare.com/v2/file/download/61624573 \
-o JAFFA_REFERENCE_FILES_hg38_gencode49.tar.gz

tar -xzf JAFFA_REFERENCE_FILES_hg38_gencode49.tar.gz
```

After extraction, the directory should contain

```
JAFFA_reference_hg38_gencode49/

    hg38.fa

    hg38_gencode49.fa

    hg38_gencode49.tab

    hg38_gencode49.bed

    hg38_gencode49.*.bt2

    Masked_hg38.*.bt2
```

If the Bowtie2 index files (*.bt2) are missing, build them before running the pipeline.

---

# Step 4. Configure nextflow.config

Update the JAFFAL reference directory

```
params {

    jaffal_ref_dir="/path/to/JAFFA_reference_hg38_gencode49"

}
```

Normally no other settings need to be changed.

---

# Step 5. Run JAFFAL

Run

```bash
nextflow run . \
    -profile rorqual \
    --summary /path/to/sequencing_summary_run1.txt \
    --fastq_dir /path/to/batch1-8
```

This workflow performs

```
pycoQC

↓

FASTCAT

↓

JAFFAL
```

Output

```
results/

    pycoqc/

    fastcat/

    jaffal/

        sample20/

            sample20.jaffa_results.csv

            sample20.jaffa_results.fasta

        sample27/

            sample27.jaffa_results.csv
```

---

# Step 6. Run the complete transcriptome workflow

The complete workflow additionally requires

- reference genome FASTA
- annotation GTF

Run

```bash
nextflow run . \
    -profile rorqual \
    --summary /path/to/sequencing_summary.txt \
    --fastq_dir /path/to/batch1-8 \
    --ref_genome /path/to/genome.fa \
    --ref_annotation /path/to/genes.gtf
```

The complete workflow performs

```
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

# Output

```
results/

    pycoqc/

    fastcat/

    jaffal/

    genome_alignment/

    transcriptome/

    differential_expression/

    annotation/

    isoform_switch/
```

---

# Important JAFFAL output

The primary fusion result is

```
sampleX.jaffa_results.csv
```

Important columns

| Column | Description |
|----------|------------|
| fusion genes | Candidate fusion |
| classification | HighConfidence / LowConfidence / PotentialTransSplicing |
| spanning reads | Number of supporting long reads |
| rearrangement | Genomic rearrangement evidence |
| known cosmic | COSMIC annotation |
| gtex samples | Number of GTEx normal samples |

---

# Troubleshooting

### "Please provide --summary"

Specify

```
--summary sequencing_summary_run1.txt
```

---

### "Please provide --fastq_dir"

Specify the directory containing the barcode folders.

---

### Exit status 137

The job ran out of memory.

Increase the memory allocation for the JAFFAL process in `nextflow.config`.

---

### JAFFAL cannot find reference files

Check

```
jaffal_ref_dir
```

in

```
nextflow.config
```

---

# Citation

If you use JAFFAL, please cite

Davidson NM et al.

**JAFFAL: detecting fusion genes with long-read transcriptome sequencing**

Genome Biology (2022)

Please also cite the original software used in this workflow, including pycoQC, FASTCAT, minimap2, StringTie, Salmon, GFFCompare, and IsoformSwitchAnalyzeR.
