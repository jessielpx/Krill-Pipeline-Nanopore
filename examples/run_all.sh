#!/bin/bash
#SBATCH --account=def-lefranco
#SBATCH --job-name=krill_91_all
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=48:00:00
#SBATCH --output=krill_91_%j.out
#SBATCH --error=krill_91_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=peixi.liu@mail.mcgill.ca

set -euo pipefail

module --force purge
module load StdEnv/2023
module load java/21.0.1
module load nextflow/25.10.2

unset PYTHONPATH PYTHONHOME PYTHONUSERBASE

cd "$SLURM_SUBMIT_DIR"

export NXF_DISABLE_CHECK_LATEST=true
export NXF_ANSI_LOG=false

mkdir -p /scratch/peixiliu/BCC_cohort91/{work,logs}

nextflow run main.nf \
  -profile rorqual \
  -params-file params.local.yaml \
  --stage all \
  -work-dir /scratch/peixiliu/BCC_cohort91/work \
  -with-report \
  -with-trace \
  -resume
