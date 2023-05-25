#!/bin/bash
#SBATCH -p short
#SBATCH --mem=8g
#SBATCH -o hem_int_des.log
source /software/pyrosetta3/setup.sh
sed -n ${SLURM_ARRAY_TASK_ID}p jobs.file.hem_int_des | bash
