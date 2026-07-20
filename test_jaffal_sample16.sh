#!/bin/bash
#SBATCH --account=def-lefranco
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --job-name=jaffal_s16
#SBATCH --output=jaffal_s16_%j.out
#SBATCH --error=jaffal_s16_%j.err

set -euo pipefail

module purge
module load StdEnv/2023
module load apptainer/1.4.5

JAFFA_SIF="/project/def-lefranco/peixiliu/NXF_SINGULARITY_CACHEDIR/jaffa-2.5.sif"

JAFFAL_REF_DIR="/home/peixiliu/links/projects/rrg-lefranco/shared/JAFFA_reference_hg38_gencode49"

FASTQ="/lustre09/project/6070433/peixiliu/Nanopore/Batch2/fastcat/sample16/merged_fastq/sample16.fastq.gz"

OUTDIR="/lustre09/project/6070433/peixiliu/Nanopore/Batch2/jaffal_test/sample16"

mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

apptainer exec \
    --cleanenv \
    --bind "${JAFFAL_REF_DIR}:/ref" \
    --bind /lustre09:/lustre09 \
    "${JAFFA_SIF}" \
    bash --noprofile --norc -c "
        set -euo pipefail

        /JAFFA/tools/bin/bpipe run \
            -p threads=${SLURM_CPUS_PER_TASK} \
            /JAFFA/JAFFAL.groovy \
            '${FASTQ}'
    "
