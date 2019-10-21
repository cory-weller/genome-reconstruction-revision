#!/usr/bin/env bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem 8G
#SBATCH -t 0-12:00:00
#SBATCH -p standard
#SBATCH --account berglandlab

module load gcc
module load R/3.5.1

N_individuals="1000"

chromosome="2L"

population="hybrid_swarm_32_F${SLURM_ARRAY_TASK_ID}"

projectDir="/scratch/${USER}/genome-reconstruction-revision/"

outputDir="/scratch/${USER}/genome-reconstruction-revision/04_RABBIT/${population}/readcounts/"

mkdir -p ${outputDir}

referenceGenome="${projectDir}/input_data/dgrp2.reference.fasta"

bamFolder="${projectDir}/02_simulate_reads/${population}/"

variants_file="${projectDir}/input_data/variants_het_2L.vcf"


getMLA.R > ${outputDir}/${ind}.mla

population <- args[1]
population_i <- 1
ind <- args[2]
chromosome <- args[3]
priors.file <- args[4]
