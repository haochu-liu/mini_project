#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=30
#SBATCH --mem-per-cpu=3700
#SBATCH --time=08:00:00

module purge
module load GCC/13.2.0 OpenMPI/4.1.6 R/4.3.3

R CMD BATCH mini_project/parallel_time.R
