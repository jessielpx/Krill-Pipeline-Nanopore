# Krill-Pipeline-Nanopore

Hi colleagues, this is the pipeline I made for Nanopore long-read RNA sequencing data. The pipeline performs quality control, genome alignment, transcriptome analysis, differential expression analysis, and gene fusion detection using JAFFAL.

The pipeline was designed for running on Rorqual, that's why I named it "Krill".

# Overview
[Pipeline workflow](https://github.com/jessielpx/Krill-Pipeline-Nanopore/blob/main/README.md#pipeline-workflow)

[Before you start](https://github.com/jessielpx/Krill-Pipeline-Nanopore/blob/main/README.md#before-you-start)

[Input files](https://github.com/jessielpx/Krill-Pipeline-Nanopore/blob/main/README.md#input-files)

[Prepare the pipeline](https://github.com/jessielpx/Krill-Pipeline-Nanopore/blob/main/README.md#prepare-the-pipeline)

[Running the pipeline](https://github.com/jessielpx/Krill-Pipeline-Nanopore/blob/main/README.md#running-the-pipeline)

[Outputs](https://github.com/jessielpx/Krill-Pipeline-Nanopore/blob/main/README.md#outputs)


# Pipeline workflow

```
Raw FASTQ files
        │
        ▼
FASTCAT
        │
        ▼
Genome alignment
        │
        ├──────────────► JAFFAL
        │                     
        ▼
Transcriptome analysis
        │
        ▼
Annotation analysis
        │
        ▼
Differential expression analysis
```

---

# Before you start

Before running the pipeline, prepare the following:

* Raw FASTQ directory
* Sequencing summary file
* Sample sheet
* Reference genome (FASTA)
* Reference annotation (GTF)
* JAFFAL reference directory
* A working directory for the sequencing run

The following sections describe each item in detail.

---

# Input files

## 1. Raw FASTQ directory

Provide the directory containing barcode folders produced after Nanopore basecalling.

Example:

```
batch1-8/

barcode01/
barcode02/
barcode03/
...
barcode96/
```

Pipeline parameter:

```
--fastq_dir
```

---

## 2. Sequencing summary

Provide the sequencing summary generated during basecalling.

Example:

```
sequencing_summary_run2_5.txt
```

Pipeline parameter:

```
--summary
```

---

## 3. Sample sheet

The sample sheet links sequencing barcodes to biological samples and experimental groups.

Pipeline parameter:

```
--samplesheet
```

Required columns:

| Column    | Description                 |
| --------- | --------------------------- |
| barcode   | Nanopore barcode directory  |
| sample    | Sample identifier           |
| alias     | Sample name used in reports |
| condition | Experimental condition      |

Example:

| barcode   | sample   | alias      | condition |
| --------- | -------- | ---------- | --------- |
| barcode13 | sample16 | Patient001 | HR        |
| barcode14 | sample72 | Patient002 | control   |

---

## 4. Reference genome

Genome FASTA used for alignment.

Pipeline parameter:

```
--ref_genome
```

Example:

```
Homo_sapiens.GRCh38.fa
```

---

## 5. Reference annotation

Genome annotation in GTF format.

Pipeline parameter:

```
--ref_annotation
```

Example:

```
genes.gtf
```

---

## 6. JAFFAL reference

JAFFAL requires a pre-built reference directory.

Pipeline parameter:

```
--jaffal_ref_dir
```

Example:

```
JAFFA_reference_hg38_gencode49/

hg38.fa
hg38_gencode49.fa
hg38_gencode49.bed
hg38_gencode49.tab
...
```

For our laboratory, this directory is available in the shared folder.

---

# Prepare the pipeline

## 1. Choose a project directory

Move to the directory where you want to store the pipeline:

```bash
cd /home/$USER/links/projects/rrg-lefranco/$USER/Nanopore
```

Create the directory if it does not already exist:

```bash
mkdir -p /home/$USER/links/projects/rrg-lefranco/$USER/Nanopore
cd /home/$USER/links/projects/rrg-lefranco/$USER/Nanopore
```

---

## 2. Clone the GitHub repository

Clone the pipeline:

```bash
git clone https://github.com/jessielpx/Krill-Pipeline-Nanopore.git
```

Enter the repository:

```bash
cd Krill-Pipeline-Nanopore
```


The repository should contain files and directories similar to:

```text
Krill-Pipeline-Nanopore/
├── main.nf
├── nextflow.config
├── README.md
├── modules/
├── subworkflows/
├── scripts/
└── run_batch2.sh
```

---

## 4. Update an existing clone

If the repository was cloned previously, do not clone it again.

Enter the existing repository:

```bash
cd /home/$USER/links/projects/rrg-lefranco/$USER/Nanopore/Krill-Pipeline-Nanopore
```

Download the latest version:

```bash
git pull origin main
```

Check the current commit:

```bash
git log --oneline -n 3
```

---

## 5. Check the pipeline configuration

Make sure all the paths to files/folders and the username are correct.

---


# Running the pipeline


A complete example SLURM submission script is [run_batch2.sh](run_batch2.sh):

Submit with:

```bash
sbatch run_batch2.sh
```

---

# Outputs

After completion, the run directory contains:

```
alignment/
```

Genome alignment files.

```
bam/
```

Final BAM files.

```
fastcat/
```

FASTCAT reports and merged FASTQ files.

```
jaffal/
```

Fusion detection results.

Each sample has its own directory containing:

```
jaffa_results.csv
jaffa_results.fasta
```

The pipeline also automatically generates:

```
merged_jaffal_results.csv
```

which combines all sample-level fusion calls into a single table.

```
transcripts/
```

Transcript assemblies.

```
quantification/
```

Expression quantification results.

```
differential_expression/
```

Differential expression analysis results.

```
logs/
```

Nextflow log files.

```
nextflow_work/
```

Intermediate Nextflow working files.

