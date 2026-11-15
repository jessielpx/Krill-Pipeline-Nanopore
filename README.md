# Krill Renew

**A modular Nextflow DSL2 workflow for Nanopore long-read cDNA RNA sequencing.**

Krill Renew integrates long-read preprocessing, genome alignment, transcript discovery, expression analysis, transcript usage, structural annotation, alternative splicing, and fusion detection.

The pipeline was developed for basal cell carcinoma (BCC) research, with particular emphasis on high-risk (HR) versus control comparisons. Its workflow architecture supports reusable sample inputs and independent analysis stages.

## 1. Overview

The complete workflow includes:

1. FASTQ merging and basic read QC
2. Genome alignment with minimap2
3. BAM indexing and alignment QC with samtools
4. Joint transcript discovery and quantification with bambu
5. SQANTI3 structural transcript QC
6. Unified sample QC
7. Differential gene expression (DESeq2)
8. Differential transcript usage (DRIMSeq)
9. Isoform usage and structural annotation
10. Functional consequence annotation
11. Gene set enrichment analysis (GSEA)
12. Junction analysis and differential junction usage
13. Long-read fusion detection with JAFFAL

Modules are connected through Nextflow channels. Independent stages can use previously generated results without repeating upstream sequencing analysis.

## 2. Repository structure

```text
Krill_Renew/
├── main.nf
├── nextflow.config
├── subworkflows/
│   ├── core.nf
│   ├── preprocessing.nf
│   ├── genome_alignment.nf
│   ├── raw_read_stages.nf
│   ├── transcriptome_stages.nf
│   ├── analysis_stages.nf
│   └── remaining_stages.nf
├── modules/
│   └── local/
├── bin/
├── assets/
└── examples/
    └── params.yaml
```

`main.nf` is the single supported workflow entrypoint.

## 3. Requirements

### Workflow environment

- Nextflow DSL2 version 25.10.2 is the currently validated version. Compatibility with Nextflow 26.x has not yet been established.
- Java compatible with the installed Nextflow version
- Linux environment
- SLURM for the provided `rorqual` profile

### Analysis software

The pipeline uses:

- minimap2
- samtools
- R and Bioconductor packages, including bambu, DESeq2 and DRIMSeq
- Python 3
- SQANTI3
- JAFFAL/JAFFA
- Additional R/Python dependencies required by individual modules

The existing configuration targets the Rorqual HPC environment using native module installations. Docker, Singularity and Apptainer are not enabled.

External installations, R libraries and reference paths must be configured by the user.

The workflow is not currently distributed as a containerized, automatically provisioned software environment.

## 4. Configuration

Copy the example configuration:

```bash
cp examples/params.yaml params.local.yaml
```

Edit `params.local.yaml` to provide your local paths, reference files, software installations and SLURM account.

Keep this local configuration private if it contains personal or institution-specific paths.

Key parameters include:

| Parameter | Purpose |
|---|---|
| `stage` | Workflow stage to execute |
| `samplesheet` | Sample metadata CSV |
| `raw_fastq_dir` | Root directory containing barcode FASTQ directories |
| `merged_fastq_dir` | Existing merged FASTQ output root |
| `bam_dir` | Existing sorted BAM output root |
| `bambu_dir` | Existing bambu results directory |
| `ref_genome` | Reference genome FASTA |
| `ref_annotation` | Reference annotation GTF |
| `outdir` | Output root |
| `slurm_account` | SLURM account |
| `r_libs_user` | R package library |
| `bioconductor_config_file` | Local Bioconductor configuration |
| `sqanti_dir` | SQANTI3 installation |
| `sqanti_env` | SQANTI3 Python environment |
| `jaffal_home` | JAFFA installation |
| `jaffal_ref_dir` | JAFFA reference files |

See `examples/params.yaml` and `nextflow.config` for additional parameters and analysis thresholds.


### Optional resource configurations

The `examples/` directory also contains:

- `bambu_resources.config`: optional increased Bambu resources
- `fusion_resources.config`: optional JAFFAL concurrency limit

These are examples, not automatically enabled. Adapt them to the
available computing resources and use `-c` when needed.

### MSigDB gene-set resources

GSEA requires three GMT files:

- `hallmark.gmt`
- `reactome.gmt`
- `go_bp.gmt`

These files are not distributed with the public workflow source.
Obtain the appropriate MSigDB collections from
https://www.gsea-msigdb.org/gsea/msigdb and prepare GMT files
containing Ensembl gene identifiers compatible with the DGE results.

Set `gsea_gmt_dir` to the directory containing these three files.
The default location is `assets/msigdb_human_ensembl/`.

The file `examples/msigdb_manifest.tsv` records the provenance
of the original project's GMT inputs, derived using msigdbr 26.1.1
and MSigDB 2026.1.Hs. It does not replace the GMT files.

Users must comply with MSigDB and underlying source licensing
and attribution requirements. Gene-set source and release versions
should be recorded for reproducibility.


## 5. Samplesheet format

The samplesheet is a CSV containing sample metadata.

Example:

```csv
barcode,sample,alias,condition,run,subtype,fastq_path
barcode1,sample01,sample01,control,1,Nodular,/data/run1/barcode1
barcode2,sample02,sample02,HR,1,Infiltrative,/data/run1/barcode2
barcode1,sample03,sample03,HR,2,Micronodular,/data/run2/barcode1
```

The raw-read parser uses `barcode`, `sample` and an optional `alias` (which defaults to the sample ID). The `condition` column is used by downstream comparisons.

Additional metadata may include `run`, `subtype` and `patient`.

Two input directory modes are supported:

**Mode A — Shared FASTQ root**

Omit `fastq_path` and supply `raw_fastq_dir`. The workflow looks for:

```text
raw_fastq_dir/
├── barcode1/
├── barcode2/
└── barcode3/
```

**Mode B — Explicit FASTQ directories**

Supply `fastq_path` for each sample. This permits barcode reuse across different sequencing runs.

Sample IDs must be unique, and two samples cannot use the same FASTQ input directory.

Currently, most differential analysis scripts are designed for the `HR` and `control` conditions. They should not be assumed to support arbitrary comparison groups without adaptation and validation.

## 6. Available stages

| Stage | Function |
|---|---|
| `all` | CORE plus all standard downstream analyses |
| `core` | FASTQ preparation, alignment, bambu, SQANTI3 and unified QC; optional downstream modules can be enabled |
| `alignment` | FASTQ preparation and genome alignment |
| `bambu` | Joint bambu analysis using existing BAM/BAI files |
| `post_bambu` | SQANTI3 and gene expression export using existing bambu outputs |
| `statistics` | DGE and DTU using existing bambu outputs |
| `isoforms` | Isoform usage and all-transcript structural annotation |
| `functional` | Functional consequence analysis from existing DTU, usage and structural annotation results |
| `gsea` | Enrichment analysis using existing DGE results |
| `junction` | Junction quantification and differential junction usage using existing BAMs |
| `fusion` | JAFFAL detection and cohort fusion aggregation using existing merged FASTQs |
| `help` | Print supported stages |

Some analyses have additional optional execution flags when running `core`. Consult `nextflow.config` for their names and defaults.

The `all` stage activates the standard downstream branches, including DGE, DTU, isoform analysis, Functional Consequences, GSEA, Junction and Fusion.

The legacy `ISOFORM_STRUCTURE` analysis is separate from the all-transcript `ISOFORM_STRUCTURE_ALL` branch.

## 7. Running the pipeline

### Inspect the pipeline

```bash
nextflow run main.nf --stage help
```

### Full workflow

After configuring all input paths and external dependencies:

```bash
nextflow run main.nf \
    -profile rorqual \
    -params-file params.local.yaml \
    --stage all \
    -w /path/to/work \
    -with-report \
    -with-trace
```

**Warning:** This command launches the computational workflow. Do not run it on an existing cohort unless a new analysis is intended.

### Independent DGE and DTU

```bash
nextflow run main.nf \
    -profile rorqual \
    -params-file params.local.yaml \
    --stage statistics \
    --bambu_dir /path/to/results/02_transcriptome/bambu_results \
    --samplesheet /path/to/samplesheet.csv \
    --outdir /path/to/new_results
```

### Independent isoform analysis

```bash
nextflow run main.nf \
    -profile rorqual \
    -params-file params.local.yaml \
    --stage isoforms \
    --bambu_dir /path/to/results/02_transcriptome/bambu_results \
    --samplesheet /path/to/samplesheet.csv \
    --ref_annotation /path/to/reference.gtf \
    --outdir /path/to/new_results
```

### DAG preview without executing analysis processes

```bash
nextflow run main.nf \
    -params-file params.local.yaml \
    --stage all \
    -preview \
    -with-dag workflow_dag.html
```

A DAG preview verifies workflow construction but does not validate runtime execution or biological results.

## 8. Output organization

Results are published under the configured `outdir`.

Important directories include:

```text
01_preprocessing/
02_transcriptome/
02_gene_expression/
03_sqanti3/
05_dge/
06_dtu/
07_isoform_usage/
08_isoform_structure_all/
09_functional_consequences/
10_gsea/
11_alternative_splicing/
12_fusion_detection/
```

Each analysis module writes its own results and QC files.

Examples include:

- `counts_gene.txt` and `counts_transcript.txt`
- `dge_HR_vs_control.tsv`
- `dtu_gene_results.tsv`
- `isoform_usage_supported.tsv`
- `isoform_structure_annotations.tsv`
- `junction_counts.tsv`
- `fusion_qc.tsv`

A subset of output directories is generated for independently selected stages.

## 9. Statistical interpretation

DGE uses DESeq2. DTU and differential junction usage use DRIMSeq-based analyses.

The supplied analysis code primarily targets HR-versus-control comparisons and includes configurable covariates and filtering thresholds.

Differential results depend on sample composition, count filtering, design matrix specification and the underlying software environment.

Fusion predictions are candidate events and require independent validation before biological or clinical interpretation.

## 10. Validation status

The modularized pipeline has undergone:

- Nextflow configuration and stage-dispatch validation
- Static process-output interface checks
- Static full-workflow dependency checks
- `all` stage DAG preview
- Samplesheet parser tests, including barcode reuse
- Existing output schema checks for downstream analysis
- Synthetic junction extraction and GTF annotation testing
- Zero-fusion aggregation testing
- Python syntax checks
- HPC software module compatibility checks

These checks do not constitute a new end-to-end execution of the entire modularized pipeline on the full cohort.

## 11. Limitations

- The provided execution profile is specific to Rorqual and SLURM.
- Software dependencies are not automatically installed.
- Downstream differential models are currently oriented toward HR/control analyses.
- Full-workflow runtime behavior should be validated when deploying in a new computing environment.
- Reference assets and third-party tools may require separate installation or licensing arrangements.

## 12. Project status

Krill Renew is a research workflow under active development.

The current repository provides a unified Nextflow entrypoint, reusable subworkflows and independent analysis stages.

