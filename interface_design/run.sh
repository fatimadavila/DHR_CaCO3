#!/bin/bash
#SBATCH -p medium
#SBATCH --mem=2g
#SBATCH -o log

# get line number ${SLURM_ARRAY_TASK_ID} from tasks file
CMD=$(sed -n "${SLURM_ARRAY_TASK_ID}p" task.list)
# tell bash to run $CMD
echo "${CMD}" | bash
#sbatch -a 1-$(cat task.list|wc -l) run.sh

