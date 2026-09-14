#!/bin/bash
#SBATCH --job-name=Simulation-SeqNMF-vs-Bregman-Iteration
#SBATCH -t 72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G

#SBATCH -A acharl15
#SBATCH --partition=parallel

#SBATCH -o Compare_MUR_SBI-%A_%a.out

#SBATCH --array=3-5
ml load matlab

# The client's --mem sets SLURM_MEM_PER_NODE, which sbatch would otherwise
# export into the worker-pool job. There it collides with the pool's
# --mem-per-cpu (SLURM_MEM_PER_CPU), making the workers' inner srun fail with:
#   "SLURM_MEM_PER_CPU, SLURM_MEM_PER_GPU, and SLURM_MEM_PER_NODE are mutually exclusive"
# Unset it here so the worker job sees only one memory spec. This does NOT
# reduce the client's memory (the cgroup limit is already applied).
unset SLURM_MEM_PER_NODE SLURM_MEM_PER_CPU

matlab -nodisplay -nosplash -nodesktop -r "run('simulation_compare_algorithms.m');exit;"
