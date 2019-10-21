#!/usr/bin/env bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=5
#SBATCH --mem 30G
#SBATCH -t 0-0:30:00
#SBATCH -p standard
#SBATCH --account berglandlab

module load gcc
module load R/3.5.1


foldername=${1}

population_size=$(echo ${foldername} | grep -o -E "[0-9]+$")
population=${foldername%_${population_size}}


Rscript calculate_STITCH_accuracy.R ${population} ${population_size}
