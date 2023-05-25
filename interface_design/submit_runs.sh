sbatch -a 1-$(cat jobs.file.hem_int_des | wc -l) sbatch_array_job.sh
