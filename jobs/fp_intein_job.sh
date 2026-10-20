#!/bin/bash
# =============================================================================
# fp_intein_job.sh -- FragPipe 24.0 semi-tryptic diaPASEF search (Dataset 3)
# Thesis: Methods 3.6. UCT ICTS HPC, account cbio, partition ada.
#
# Usage:  sbatch fp_intein_job.sh dry        # check setup only (--dry-run)
#         sbatch fp_intein_job.sh            # real search
# A small dry run starts faster:  sbatch --cpus-per-task=2 --mem=8G fp_intein_job.sh dry
#
# Reconstructed from the session record of 14 Oct 2026, including the later
# change of --config-tools-folder. Compare with the copy on the cluster:
#   diff fp_intein_job.sh $PROJ/fragpipe/fp_intein_job.sh
# =============================================================================
#SBATCH --account=cbio
#SBATCH --partition=ada
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=120G
#SBATCH --time=48:00:00
#SBATCH --job-name=fp_intein
#SBATCH --output=fp_intein_%j.out

# Set PROJ to the directory holding the FragPipe run, the search
# database and the results, for example:
#   export PROJ=/scratch/$USER/TB_rif_DIA
: "${PROJ:?set PROJ to your project directory}"
set -euo pipefail                      # any failure stops the job and is reported
cd $PROJ/fragpipe
EXTRA=""
if [ "${1:-run}" = "dry" ]; then EXTRA="--dry-run"; fi
echo "Start: $(date)  Mode: ${1:-run}  Node: $(hostname)"
# The tools folder must be FragPipe's own tools/ (it holds speclib/ scripts);
# MSFragger 4.4.1, IonQuant 1.11.20 and diaTracer 2.2.1 were copied into it.
~/software/fragpipe-24.0/bin/fragpipe --headless \
  --workflow intein_semi_diaPASEF.workflow \
  --manifest intein_semi.manifest \
  --workdir $PROJ/fragpipe/results \
  --ram 110 --threads 20 \
  --config-tools-folder ~/software/fragpipe-24.0/tools \
  --config-diann ~/software/fragpipe-24.0/tools/diann/1.8.2_beta_8/linux/diann-1.8.1.8 \
  --config-python ~/envs/fragpipe/bin/python3 \
  $EXTRA
echo "End: $(date)"
