#!/usr/bin/env bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=10
#SBATCH --mem 24G
#SBATCH -t 0-12:00:00
#SBATCH -p standard
#SBATCH --account berglandlab

module purge
module load gcc
module load openmpi
module load R/3.5.3

nFounders=${1}

nGenerations=${2}

Rscript forward_simulator.Rscript \
-bed recombination.bed \
-prefix dgrp_F50_outbred \
-n0 10000 \
-rate 1.0 \
-sex dioecious \
-nfounders ${nFounders} \
-ngenerations ${nGenerations} \
-lineIDs lines.txt \
-chrx X \
-iter ${SLURM_ARRAY_TASK_ID} \
-recombination femaleOnly \
-dmel TRUE \
-nthreads 10 \
#-nRILs 800 \
#-inbreed_generations 25 \
