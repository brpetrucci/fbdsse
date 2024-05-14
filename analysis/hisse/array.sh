#!/bin/bash

#SBATCH --nodes=1 # 1 node per job
#SBATCH --time=01:00:00 # just an hour--calling the array itself doesn't take long 
#SBATCH --array=1-49 # 114 parameter combinations, 50 at a time

#SBATCH --output=output/hisse/jobs/array/job_%A_%a.out
#SBATCH --error=output/hisse/jobs/array/job_%A_%a.err

#SBATCH --job-name="array_hisse"

#SBATCH --mail-user=petrucci@iastate.edu   # my e-mail
#SBATCH --mail-type=BEGIN # get notifications for all job cases
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

mkdir output/hisse/jobs/job_${SLURM_ARRAY_TASK_ID}
sbatch --time=5-00:00:00 --array=1-100 --output="output/hisse/jobs/job_${SLURM_ARRAY_TASK_ID}/hisse_${SLURM_ARRAY_TASK_ID}_%A_%a.out" --error="output/hisse/jobs/job_${SLURM_ARRAY_TASK_ID}/hisse_${SLURM_ARRAY_TASK_ID}_%A_%a.err" --mail-user=petrucci@iastate.edu --mail-type=BEGIN --mail-type=END --mail-type=FAIL --job-name="${SLURM_ARRAY_TASK_ID}_array" --wrap="sh analysis/hisse/hisse.sh ${SLURM_ARRAY_TASK_ID}"
