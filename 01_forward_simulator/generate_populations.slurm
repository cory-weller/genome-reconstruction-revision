#!/usr/bin/env bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --mem 16G
#SBATCH -t 0-24:00:00
#SBATCH -p standard
#SBATCH --account berglandlab

# Edit this file to contain the parameters desired for running the forward simulator.
# Exclude -nRILs and -inbreed_generations if not generating RILs. See README.md for more detail.

module purge
module load gcc
module load openmpi
module load R/3.5.3

Rscript forward_simulator.R \
-bed ../input_data/recombination_map.bed \
-prefix hybrid_swarm \
-n0 10000 \
-rate 1.0 \
-sex dioecious \
-nfounders 32 \
-ngenerations 1 \
-lineIDs lines.txt \
-chrx X \
-iter 1 \
-recombination femaleOnly \
-dmel TRUE \
-nthreads 4 \
# -nRILs 800 \
# -inbreed_generations 25 \

Rscript forward_simulator.R \
-bed ../input_data/recombination_map.bed \
-prefix hybrid_swarm \
-n0 10000 \
-rate 1.0 \
-sex dioecious \
-nfounders 32 \
-ngenerations 2 \
-lineIDs lines.txt \
-chrx X \
-iter 1 \
-recombination femaleOnly \
-dmel TRUE \
-nthreads 4 \
# -nRILs 800 \
# -inbreed_generations 25 \

Rscript forward_simulator.R \
-bed ../input_data/recombination_map.bed \
-prefix hybrid_swarm \
-n0 10000 \
-rate 1.0 \
-sex dioecious \
-nfounders 32 \
-ngenerations 5 \
-lineIDs lines.txt \
-chrx X \
-iter 1 \
-recombination femaleOnly \
-dmel TRUE \
-nthreads 4 \
# -nRILs 800 \
# -inbreed_generations 25 \

Rscript forward_simulator.R \
-bed ../input_data/recombination_map.bed \
-prefix hybrid_swarm \
-n0 10000 \
-rate 1.0 \
-sex dioecious \
-nfounders 32 \
-ngenerations 50 \
-lineIDs lines.txt \
-chrx X \
-iter 1 \
-recombination femaleOnly \
-dmel TRUE \
-nthreads 4 \
# -nRILs 800 \
# -inbreed_generations 25 \
