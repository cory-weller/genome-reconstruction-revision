# Readme

## Arguments:

* `-bed` is the (string) path to a file containing recombination rates.
* `-prefix` is a string to identify the population being simulated.
* `-n0` is an integer for the initial population size.
* `-rate` is the per-generation population growth rate. Use `1.0` to remain constant.
* `-sex` can be `dioecious` or `hermaphroditic`.
* `-nfounders` sets the number of unique lines (e.g. columns from VCF file) that found the population.
* `-ngenerations` sets the integer count of generations for recombination, e.g. `1` for F1, `5` for F5.
* `-lineIDs` is the (string) path to input text file containing (one per row) Line IDs.
* `-chrx` is the string used to identify the female sex chromosome.
* `-iter` is an integer identifier for unique replicate simulated populations.
* `-nthreads` is number of parallel threads.

* `-recombination femaleOnly` will prevent recombination in males, use `-recombination both` for recombination in both parents.
* `-dmel TRUE` will cause 2L/2R and 3L/3R to function as a single chromosome for recombination purposes. Set `FALSE` otherwise.
* `-nRILs` sets the number of recombinant inbred lines to generate, following the outbreeding step. If used, also requires `-inbreed_generations`. Exclude if not generating RILs.
* `-inbreed_generations` sets the number of generations of full-sib mating, following the outbreeding steps. If used, also requires `-nRILs`. Exclude if not generating RILs.


## Reproducibility
Each simulation sets a random seed based on the parameters used.


## Running Forward Simulator
Submit a job to SLURM using a command similar to the following heredoc:
```
sbatch --array=1-1%1 --ntasks-per-node=1 --mem=4G --time=0-0:05:00 --partition=dev --account=berglandlab --nodes=1 <<EOF
#!/usr/bin/env bash

module purge
module load gcc
module load openmpi
module load R/3.5.3

Rscript forward_simulator.Rscript \
-bed recombination.bed \
-prefix sbatch_test \
-n0 100 \
-rate 1.0 \
-sex dioecious \
-nfounders 32 \
-ngenerations 2 \
-lineIDs lines.txt \
-chrx X \
-iter \${SLURM_ARRAY_TASK_ID} \
-recombination femaleOnly \
-dmel TRUE \
-nthreads 1 \
#-nRILs 800 # uncomment start of line if generating RILs \
#-inbreed_generations 25 # uncomment start of line if generating RILs \

EOF
```
